!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2018-2023 Philipp Pracht, Stefan Grimme
!
! crest is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! crest is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with crest.  If not, see <https://www.gnu.org/licenses/>.
!================================================================================!

!=========================================================================================!
! estimate conformational flexibility from geom. and WBO
! the value of flex is normalized i.e. between 0(rigid) and 1(long alkane)
!=========================================================================================!
subroutine flexi(mol,rednat,includeAtom,flex,nfix,fixat)
  use iso_fortran_env,only:wp => real64
  use crest_parameters,only:stdout,autoaa
  use strucrd
  use zdata,only:readwbo
  use miscdata, only: rcov
  use crest_cn_module
  implicit none
  !> INPUT
  type(coord),intent(in) :: mol                  !> input molecule
  integer,intent(in)     :: rednat               !> should equal sum(includeAtom(:))
  integer,intent(in)     :: includeAtom(mol%nat) !> 0 if not considered, 1 if considered
  !> OPTIONAL INPUT: bonds the user declared explicitly (--cinp / TOML
  !> "constraints").  They are inserted into the bond graph EVEN IF THEY ARE
  !> LONGER THAN THE COVALENT RADIUS SUM, so that a bond which is still forming
  !> (or already breaking) in a transition state still closes its ring.  The
  !> only guard is CREST_FLEXI_FORCEFAC (default 1.5), see below.  NOTE that a
  !> "atoms: 1,81,82,86,87" block registers ALL C(5,2) = 10 distances, not 5,
  !> and most of them are not bonds - the gate is what filters those out.
  integer,intent(in),optional :: nfix            !> number of declared bonds
  integer,intent(in),optional :: fixat(:,:)      !> (2,nfix) 1-based FULL indices
  !> OUTPUT
  real(wp),intent(out) :: flex
  !> LOCAL
  real(wp),allocatable::xyz(:,:),cn(:),cring(:,:),wbo(:,:),wbofull(:,:)
  integer,allocatable::at(:),map(:),b(:,:),sring(:),map2(:)
!>--- LOCAL EXTENSION (2026-09-29): reverse map (reduced index -> FULL 1-based
!> index) so the diagnostic table can quote the numbers the user's own xyz file
!> uses.  See the construction below, right after map2.
  integer,allocatable::irev(:)
  real(wp) :: dx,dy,dz,r,r2,rco,val,ringf,doublef,branch,effectivNat,hybf
  real(wp) :: av2
  integer :: i,j,k,l,m,n,rn
  logical :: ex

  !>--- diagnostic dump (CREST_FLEXI_DEBUG=1), see the bottom of this routine
  !> NOTE: the size parameter must NOT be called NTOP - Fortran is
  !> case-insensitive, so NTOP and the running counter `ntop` below are the
  !> same symbol (that collision does not compile).
  integer,parameter :: NTOPLIST = 20   !> how many bonds to list
  integer :: dbg,rd2,kk,ntop,nrb,nab,srh,itmp,jtmp,npart
  character(len=32) :: envstr
  integer :: ti(NTOPLIST),tj(NTOPLIST)
  real(wp) :: tv(NTOPLIST),vtmp,ssum,sr2,sa2
  integer :: srhn(0:14)

  !>--- continuous ring membership (CREST_FLEXI_RINGMODE=1, the default)
  integer :: ringmode,iu,iv,it,nhop
  logical :: inring
  real(wp) :: ww,nd,cstr,ccho
  !> winsorising cap on val (CREST_FLEXI_VALCAP) and power-mean order
  !> (CREST_FLEXI_PMEAN).  vcap = 1.0 disables the cap, pmean = 2.0 restores
  !> the original root-mean-square.
  real(wp) :: vcap,pmean
  real(wp),allocatable :: wt(:,:),rho(:,:),dst(:)
  integer,allocatable :: nc(:,:),hps(:)
  logical,allocatable :: vis(:)
  !> user declared (--cinp) bonds, remapped onto the reduced index space
  integer :: nfr,lf,ii,jj
  integer,allocatable :: fcat(:,:)
  logical,allocatable :: frc(:,:)
  !> CREST_FLEXI_FORCEFAC: sanity gate for the declared bonds.  A declared bond
  !> is only inserted if  dist <= FORCEFAC * (rcov_i+rcov_j).  The default of
  !> 1.3 accepts every forming/breaking bond one ever sees at a transition
  !> state (1.3 x (rcov_C+rcov_C) = 2.6 Ang, 1.3 x (rcov_Rh+rcov_C) = 3.3 Ang)
  !> while it drops the pairs that a  "atoms: 1,81,82,86,87"  block produces ON
  !> TOP of the real bonds: that block registers ALL C(5,2) = 10 distances, and
  !> e.g. 1-81 at 5.2 Ang or 86-87 at 2.9 Ang would otherwise become bonds and
  !> invent rings that do not exist.  FORCEFAC <= 0 disables the gate, i.e.
  !> every single declared distance becomes a bond, whatever its length.
  real(wp) :: ffac
  integer :: nrej
  !> steepness of the pair bond order sigmoid (identical to calculate_CN)
  real(wp),parameter :: WSIG = 16.0_wp
  !> Cycle coherence of a fully locked ring, i.e. the geometric mean bond order
  !> around it.  A normal covalent ring sits at ~0.99 and therefore saturates at
  !> rho = 1, which reproduces the OLD ringf exactly (no recalibration needed).
  !> Strained or half-open rings fall below it and get a partial penalty instead
  !> of the full one, which is what removes the cliff.
  real(wp),parameter :: RHOFULL = 0.98_wp
  !> ring size damping, unchanged from the original implementation
  real(wp),parameter :: RINGK = 0.06_wp
  !> "unreachable" cost for the path search
  real(wp),parameter :: BIGC = 1.0e30_wp

  !> CN less than this is considerd as terminating atom (H, F, ...)
  real(wp),parameter :: thr = 1.2_wp

  flex = 0.0_wp

!>--- read the diagnostic switch BY VALUE (a present-but-unparsable value must
!>    fall back to 0, and a failed internal read leaves dbg undefined)
  dbg = 0
  call get_environment_variable('CREST_FLEXI_DEBUG',value=envstr,status=rd2)
  if (rd2 == 0) then
    if (len_trim(envstr) > 0) read (envstr,*,iostat=rd2) dbg
    if (rd2 /= 0) dbg = 0
  end if
!>--- read the ring-membership switch BY VALUE (same convention as above).
!>    ringmode = 0 restores the original discrete sring() based ringf.
  ringmode = 1
  call get_environment_variable('CREST_FLEXI_RINGMODE',value=envstr,status=rd2)
  if (rd2 == 0) then
    if (len_trim(envstr) > 0) read (envstr,*,iostat=rd2) ringmode
    if (rd2 /= 0) ringmode = 1
  end if
!>--- winsorising cap and power-mean order, both read BY VALUE.
!>    The MTD length is re-tuned by the user via --mdlen, so no attempt is
!>    made here to recalibrate the exp() model in choose_settings.f90.
  vcap = 0.9_wp
  call get_environment_variable('CREST_FLEXI_VALCAP',value=envstr,status=rd2)
  if (rd2 == 0) then
    if (len_trim(envstr) > 0) read (envstr,*,iostat=rd2) vcap
    if (rd2 /= 0) vcap = 0.9_wp
  end if
  if (vcap .le. 0.0_wp) vcap = 1.0_wp   !> <=0 would zero out everything
  pmean = 1.5_wp
  call get_environment_variable('CREST_FLEXI_PMEAN',value=envstr,status=rd2)
  if (rd2 == 0) then
    if (len_trim(envstr) > 0) read (envstr,*,iostat=rd2) pmean
    if (rd2 /= 0) pmean = 1.5_wp
  end if
  if (pmean .lt. 0.1_wp) pmean = 0.1_wp
  ffac = 1.3_wp
  call get_environment_variable('CREST_FLEXI_FORCEFAC',value=envstr,status=rd2)
  if (rd2 == 0) then
    if (len_trim(envstr) > 0) read (envstr,*,iostat=rd2) ffac
    if (rd2 /= 0) ffac = 1.3_wp
  end if
  if (ffac .lt. 0.0_wp) ffac = 0.0_wp     !> <= 0 : force EVERY declared pair
  ntop = 0
  ti = 0
  tj = 0
  tv = 0.0_wp
  nrb = 0
  nab = 0
  npart = 0
  sr2 = 0.0_wp
  sa2 = 0.0_wp

  n = mol%nat !> local modification
  rn = rednat !> Number of atoms selected by atomlist+/-
  !> rn=n if atomlist+/- is not used

  allocate (xyz(3,rn),at(rn),map(rn),cn(rn),b(rn,rn))
  allocate (sring(rn),cring(12,rn),wbo(rn,rn),wbofull(n,n),map2(n))

!>--- map given structure to considered substructure
  if (rn .ne. n) then
    call rdcoord_reduced('coord',n,rn,xyz,at,includeAtom)
  else
    !call rdcoord('coord',rn,at,xyz)
    at = mol%at
    xyz = mol%xyz
  end if
  call calculate_CN(rn,at,xyz,cn)

!>--- map the new (reduced) coordinate order to the original
  j = 0
  map2 = 0   !> every EXcluded atom should get a zero
  do i = 1,n
    if (includeAtom(i) .gt. 0) then
      j = j+1
      map2(i) = j
    end if
  end do

!>--- LOCAL EXTENSION (2026-09-29): the INVERSE map, reduced -> original.
!> Why it is needed: every array this routine works with (at, xyz, cn, wbo, b,
!> sring, ...) is indexed in the REDUCED space, so when the diagnostic table
!> below reports `i` and `j` it is quoting reduced indices.  With nat=97 and
!> rednat=92 (a --cinp constraint locks 5 atoms) the printed number is therefore
!> NOT the number the user sees in his own xyz file - phosphorus sits at 50/51
!> in the input but the table called it 47/48, and the offset GROWS with the
!> position (C 60 -> 64, C 71 -> 75, C 82 -> 86), which is exactly why it is
!> easy to misread as "almost right".
!> fixat, by contrast, is documented as 1-based FULL indices, so the table was
!> the only place that leaked the reduced space to the user.  Print
!> irev(i) instead; if the reduction is disabled (rn == n) it is the identity.
!> NOTE the counter must not be called `i` - that is the enclosing loop's
!> variable and reusing it here would corrupt the caller's iteration.
  allocate (irev(rn))
  do i = 1,rn
    irev(i) = i          !> identity when nothing was excluded
  end do
  k = 0
  do i = 1,n
    if (includeAtom(i) .gt. 0) then
      k = k+1
      irev(k) = i
    end if
  end do

!>--- read the WBO file and map it to the (reduced) coordinates
!>--- if it does not exist, single and double bonds are not distinguished
  wbo = 0
  wbofull = 0
  inquire (file='wbo',exist=ex)
  if (ex) then
    call readwbo('wbo',n,wbofull)

!>--- map wbofull to wbo (only if atoms are excluded from the rmsd)
    if (rn .ne. n) then
      do i = 1,n-1
        if (includeAtom(i) .gt. 0) then
          k = map2(i)
          do j = i+1,n
            if (includeAtom(j) .gt. 0) then
              l = map2(j)
              wbo(k,l) = wbofull(i,j)
              wbo(l,k) = wbo(k,l)
            else
              cycle !> j cycle
            end if
          end do
        else
          cycle !> i cycle
        end if
      end do
    else !> rn == n
      wbo = wbofull
    end if !>-- if(rn.ne.n) end
  end if !>-- if(ex) end

!>--- remap the user declared bonds onto the (reduced) coordinate order.
!>    This has to happen here, while map2 is still alive; pairs that touch an
!>    excluded atom are dropped.
  nfr = 0
  nrej = 0
  if (present(nfix)) then
    if (nfix .gt. 0.and.present(fixat)) then
      allocate (fcat(2,nfix),source=0)
      do lf = 1,nfix
        ii = fixat(1,lf)
        jj = fixat(2,lf)
        if (ii .lt. 1.or.ii .gt. n) cycle
        if (jj .lt. 1.or.jj .gt. n) cycle
        if (rn .ne. n) then
          if (map2(ii) .eq. 0.or.map2(jj) .eq. 0) cycle
          ii = map2(ii)
          jj = map2(jj)
        end if
        if (ffac .gt. 0.0_wp) then
          r = sqrt((xyz(1,ii)-xyz(1,jj))**2+(xyz(2,ii)-xyz(2,jj))**2+ &
          &        (xyz(3,ii)-xyz(3,jj))**2)
          rco = rcov(at(ii))+rcov(at(jj))
          if (r .gt. ffac*rco) then
            nrej = nrej+1
            if (dbg .gt. 0) then
              write (stdout,'(5x,a,i0,a,i0,a,f7.3,a,f6.2,a)') 'skipped ', &
              &  fixat(1,lf),' - ',fixat(2,lf),'  r = ',r*autoaa, &
              &  ' Ang  (r/rcov = ',r/max(1.0e-6_wp,rco),' > FORCEFAC)'
            end if
            cycle
          end if
        end if
        nfr = nfr+1
        fcat(1,nfr) = ii
        fcat(2,nfr) = jj
      end do
    end if
  end if
  deallocate (map2,wbofull)

!+++++++++++++++++++++++++++++++++++++++!
  n = rn !<---- so I don't have to change the following original code
!+++++++++++++++++++++++++++++++++++++++!

!>--- adjacency matrix b setup
  map = 0
  b(:,:) = 0
  do i = 1,n
    if (cn(i) .lt. thr) cycle
    do j = 1,n
      if (i .eq. j) cycle
      dx = xyz(1,j)-xyz(1,i)
      dy = xyz(2,j)-xyz(2,i)
      dz = xyz(3,j)-xyz(3,i)
      r2 = dx*dx+dy*dy+dz*dz
      r = sqrt(r2)
      rco = rcov(at(i))+rcov(at(j))
      if (r .lt. rco) then
        if (cn(j) .gt. thr) then
          b(j,i) = 1
          map(i) = map(i)+1
        end if
      end if
    end do
  end do

!>--- insert the user declared bonds.  They are accepted regardless of the
!>    actual distance: in a transition state the forming/breaking bonds are
!>    exactly the ones that are too long (or too short) for the geometric
!>    criterion above, and they are the ones that must close the ring.
  allocate (frc(n,n),source=.false.)
  do lf = 1,nfr
    i = fcat(1,lf)
    j = fcat(2,lf)
    if (i .lt. 1.or.i .gt. n) cycle
    if (j .lt. 1.or.j .gt. n) cycle
    if (b(j,i) .eq. 0) then
      b(j,i) = 1
      b(i,j) = 1
      map(i) = map(i)+1
      map(j) = map(j)+1
    end if
    frc(i,j) = .true.
    frc(j,i) = .true.
  end do

!>--- setup ring data for all atoms, i.e., in which ring and how large
  call minringsizes(n,at,xyz,sring)

!>=========================================================================!
!> continuous ring membership rho(i,j) and ring size nc(i,j)
!>
!> The original implementation asks the topology module "is this bond in a
!> ring? yes/no" and then jumps between ringf = 1.0 and ringf = 0.151.  That
!> binary answer comes from simpletopo(), which rounds the coordination number
!> (nint(cn2)) before it picks neighbours, so a tiny geometry change can flip
!> it - and with it the whole flexibility measure.
!>
!> Instead of a yes/no answer we measure HOW STRONGLY the bond is locked:
!> for every bond i-j we search the strongest path from i back to j that does
!> NOT use the i-j bond itself (that path closes the cycle).  The strength of
!> a path is the product of the smooth pair bond orders
!>      w = 1/(1+exp(-16*(rco/r-1)))      (the sigmoid of calculate_CN)
!> so it decays continuously when a contact is stretched, instead of switching
!> on and off.  rho is the geometric mean of w along the cycle, normalised by
!> RHOFULL.  A well formed covalent ring gives rho = 1 (identical to the old
!> ringf), a distorted or half-open ring gives something in between, and an
!> acyclic bond gives rho = 0 (free rotation, as before).
!>=========================================================================!
  allocate (rho(n,n),nc(n,n))
  rho = 0.0_wp
  nc = 0
!>--- the path search is O(m*n^2).  For very large systems fall back to the
!>    original discrete ring term rather than slowing the setup down.
  if (ringmode .gt. 0 .and. n .gt. 1000) ringmode = 0
  if (ringmode .gt. 0) then
    allocate (wt(n,n),dst(n),hps(n),vis(n))
    wt = 0.0_wp
    do i = 1,n
      do j = 1,i-1
        if (b(j,i) .le. 0) cycle
        r = sqrt((xyz(1,i)-xyz(1,j))**2+(xyz(2,i)-xyz(2,j))**2+(xyz(3,i)-xyz(3,j))**2)
        rco = rcov(at(i))+rcov(at(j))
        ww = 1.0_wp/(1.0_wp+exp(-WSIG*(rco/r-1.0_wp)))
!>------ a declared bond is a bond: give it full weight, otherwise the long
!>       forming bond would drag the cycle coherence down and the ring would
!>       still not be recognised
        if (frc(i,j)) ww = 1.0_wp
        if (ww .lt. 1.0e-8_wp) ww = 1.0e-8_wp
        wt(i,j) = ww
        wt(j,i) = ww
      end do
    end do
    do i = 1,n
      do j = 1,i-1
        if (b(j,i) .le. 0) cycle
        if (map(i) .eq. 1.or.map(j) .eq. 1) cycle
!>------ Dijkstra on -log(w): maximising a product == minimising a sum
        dst = BIGC
        hps = 0
        vis = .false.
        dst(i) = 0.0_wp
        do it = 1,n
          iu = 0
          nd = BIGC
          do iv = 1,n
            if (.not.vis(iv)) then
              if (dst(iv) .lt. nd) then
                nd = dst(iv)
                iu = iv
              end if
            end if
          end do
          if (iu .eq. 0) exit
          vis(iu) = .true.
          if (iu .eq. j) exit
          do iv = 1,n
            if (iv .eq. iu) cycle
            if (wt(iu,iv) .le. 0.0_wp) cycle
            if (iu .eq. i.and.iv .eq. j) cycle   !> the bond itself
            if (iu .eq. j.and.iv .eq. i) cycle
            nd = dst(iu)-log(wt(iu,iv))
            if (nd .lt. dst(iv)) then
              dst(iv) = nd
              hps(iv) = hps(iu)+1
            end if
          end do
        end do
        if (dst(j) .lt. BIGC*0.5_wp) then
          cstr = exp(-dst(j))                    !> product of w along the cycle
          nhop = hps(j)+1                        !> + the i-j bond = ring size
          ccho = cstr**(1.0_wp/dble(max(1,nhop))) !> geometric mean, size free
          rho(i,j) = min(1.0_wp,ccho/RHOFULL)
          nc(i,j) = nhop
        end if
        rho(j,i) = rho(i,j)
        nc(j,i) = nc(i,j)
      end do
    end do
    deallocate (wt,dst,hps,vis)
  end if

!>--- the actual flexibility measure calculation
  m = 0
  av2 = 0.0d0
  do i = 1,n
    do j = 1,i-1
      !> count all "bonds"
      if (b(j,i) .gt. 0) m = m+1

      !> no branch on terminating bonds (e.g. Me)
      if (map(i) .eq. 1.or.map(j) .eq. 1) cycle

      !> adjacent pair ij
      if (b(j,i) .gt. 0) then

        hybf = 1.0

        !> sp2 C are less flexible
        if (at(i) .eq. 6.and.cn(i) .lt. 3.3) hybf = hybf*0.5_wp
        if (at(j) .eq. 6.and.cn(j) .lt. 3.3) hybf = hybf*0.5_wp

        !> double bond term
        doublef = 1.0-exp(-4.0*(wbo(j,i)-2.0)**6)

        !> branching measure
        branch = 2.0/sqrt(dble(map(i))*dble(map(j)))

        !> flex. reduction for atoms in rings
        if (ringmode .gt. 0) then
!>-------- continuous ring membership: 0 = free rotation, 1 = locked in a ring
          inring = nc(i,j) .gt. 0
          ringf = 1.0_wp
          if (inring) then
            ringf = 1.0_wp-(1.0_wp-0.5_wp*(1.0_wp-exp(-RINGK*dble(nc(i,j)))))*rho(i,j)
          end if
        else
!>-------- original (discrete) version, kept for A/B comparison
          ringf = 1.0
          k = min(sring(i),sring(j))
          inring = k .gt. 0
          !> a ring is (even at infinite size) a factor of 2 less flexible, factors are empirical
          if (inring) ringf = 0.5*(1.0-exp(-RINGK*float(k)))
          !> (adjusted such that c-C20 is converged)
        end if

        !>  put it together
        val = branch*ringf*doublef*hybf
!>------ winsorise: no single bond may dominate the average
        if (val .gt. vcap) val = vcap
        !> power mean of order pmean (pmean = 2 would be the original RMS)
        av2 = av2+val**pmean

!>--- diagnostic: keep the NTOPLIST largest contributions (insertion sort)
        if (dbg > 0) then
          if (inring) then
            nrb = nrb+1
            sr2 = sr2+val**pmean
            if (ringmode .gt. 0) then
              if (rho(i,j) .gt. 1.0e-6_wp.and.rho(i,j) .lt. 1.0_wp-1.0e-6_wp) then
                npart = npart+1
              end if
            end if
          else
            nab = nab+1
            sa2 = sa2+val**pmean
          end if
          if (ntop .lt. NTOPLIST) then
            ntop = ntop+1
            ti(ntop) = i
            tj(ntop) = j
            tv(ntop) = val
          else if (val .gt. tv(ntop)) then
            ti(ntop) = i
            tj(ntop) = j
            tv(ntop) = val
          end if
          do kk = ntop,2,-1
            if (tv(kk) .le. tv(kk-1)) exit
            itmp = ti(kk)
            ti(kk) = ti(kk-1)
            ti(kk-1) = itmp
            jtmp = tj(kk)
            tj(kk) = tj(kk-1)
            tj(kk-1) = jtmp
            vtmp = tv(kk)
            tv(kk) = tv(kk-1)
            tv(kk-1) = vtmp
          end do
        end if
      end if
    end do
  end do
  if (m .gt. 0) then
    av2 = (av2/dble(m))**(1.0_wp/pmean)
  end if

  flex = av2

!>=========================================================================!
!> optional diagnostic dump: show how this flexibility number came about.
!> Everything is still allocated here, so the per-bond factors can simply be
!> recomputed for the listed bonds instead of being stored.
!>=========================================================================!
  if (dbg > 0) then
    ssum = flex**pmean*dble(m)
    srhn = 0
    do i = 1,n
      srh = sring(i)
      if (srh .gt. 14) srh = 14
      srhn(srh) = srhn(srh)+1
    end do
    write (stdout,'(1x,a)') '--- covalent flexibility: how this number was obtained ---'
    write (stdout,'(3x,a,i0,a,i0)') 'nat / reduced nat            : ',mol%nat,' / ',n
    write (stdout,'(3x,a,l1)') 'wbo file found and read      : ',ex
    write (stdout,'(3x,a,i0)') 'user declared bonds (--cinp) : ',nfr
    if (nfr .gt. 0) then
      do lf = 1,nfr
        ii = fcat(1,lf)
        jj = fcat(2,lf)
        r = sqrt((xyz(1,ii)-xyz(1,jj))**2+(xyz(2,ii)-xyz(2,jj))**2+ &
        &        (xyz(3,ii)-xyz(3,jj))**2)
        rco = rcov(at(ii))+rcov(at(jj))
        if (r .gt. rco) then
          write (stdout,'(5x,a,i0,a,i0,3x,a,f7.3,a,f6.3,a)') 'forced ',ii, &
          &  ' - ',jj,' r = ',r*autoaa,' Ang   r/rcov = ',r/rco, &
          &  '   (long: forming/breaking)'
        else
          write (stdout,'(5x,a,i0,a,i0,3x,a,f7.3,a,f6.3)') 'forced ',ii, &
          &  ' - ',jj,' r = ',r*autoaa,' Ang   r/rcov = ',r/rco
        end if
      end do
    end if
    if (nrej .gt. 0) then
      write (stdout,'(3x,a,i0,a,f6.2,a)') 'declared pairs SKIPPED       : ', &
      &  nrej,'   (r > ',ffac,' x (rcov_i+rcov_j), see CREST_FLEXI_FORCEFAC)'
    end if
    write (stdout,'(3x,a,i0)') 'bonds m (both ends CN>1.2)   : ',m
    write (stdout,'(3x,a,i0,a)') 'bonds entering the average   : ',nrb+nab, &
    &                            '   (terminal bonds skipped)'
    write (stdout,'(3x,a,i0,a,i0)') '  of which flagged as ring   : ',nrb, &
    &                               '   acyclic                   : ',nab
    if (ringmode .gt. 0) then
      write (stdout,'(3x,a,i0)') 'partially locked (0<rho<1)   : ',npart
    end if
    write (stdout,'(3x,a,f6.2)') 'winsorising cap on val       : ',vcap
    write (stdout,'(3x,a,f6.2)') 'power mean order p           : ',pmean
    write (stdout,'(3x,a,f12.6)') 'sum of val^p                 : ',ssum
    write (stdout,'(3x,a,f12.6)') '  contributed by ring bonds  : ',sr2
    write (stdout,'(3x,a,f12.6)') '  contributed by acyclic     : ',sa2
    write (stdout,'(3x,a,f9.6)') 'flex = (sum/m)^(1/p)         : ',flex
    write (stdout,'(3x,a)') 'atoms per smallest-ring size (topology module view, 0 = not in a ring):'
    do kk = 0,14
      if (srhn(kk) > 0) write (stdout,'(5x,a,i0,a,i0)') 'size ',kk,' : ',srhn(kk)
    end do
    write (stdout,'(3x,a)') 'largest contributions to sum(val^p):'
    write (stdout,'(5x,a)') &
    &  'bond       r/Ang  branch  ringf     rho   hybf  doublef     wbo  valraw  valcap  share%  fix'
    do kk = 1,ntop
      i = ti(kk)
      j = tj(kk)
      r = sqrt((xyz(1,i)-xyz(1,j))**2+(xyz(2,i)-xyz(2,j))**2+(xyz(3,i)-xyz(3,j))**2)
      hybf = 1.0_wp
      if (at(i) .eq. 6.and.cn(i) .lt. 3.3) hybf = hybf*0.5_wp
      if (at(j) .eq. 6.and.cn(j) .lt. 3.3) hybf = hybf*0.5_wp
      doublef = 1.0_wp-exp(-4.0_wp*(wbo(j,i)-2.0_wp)**6)
      branch = 2.0_wp/sqrt(dble(map(i))*dble(map(j)))
      if (ringmode .gt. 0) then
        ringf = 1.0_wp
        if (nc(i,j) .gt. 0) then
          ringf = 1.0_wp-(1.0_wp-0.5_wp*(1.0_wp-exp(-RINGK*dble(nc(i,j)))))*rho(i,j)
        end if
        ccho = rho(i,j)
      else
        k = min(sring(i),sring(j))
        ringf = 1.0_wp
        if (k .gt. 0) ringf = 0.5_wp*(1.0_wp-exp(-RINGK*float(k)))
        ccho = 0.0_wp
      end if
      vtmp = branch*ringf*doublef*hybf        !> before the cap
!>--- LOCAL EXTENSION (2026-09-29): report FULL (input-file) atom numbers, not
!> the reduced-order indices.  i/j index the reduced arrays; irev() maps them
!> back to the numbering the user sees in his own xyz file (identity when no
!> atom was excluded).  Printing i/j here is what made phosphorus at 50/51 show
!> up as 47/48.
      write (stdout,'(5x,a2,i0,"-",a2,i0,2x,f7.3,3f8.3,f7.2,f9.4,f8.3,2f8.3,f8.1,2x,l1)') &
      &  asym(at(i)),irev(i),asym(at(j)),irev(j),r*autoaa,branch,ringf,ccho,hybf,doublef,wbo(j,i), &
      &  vtmp,tv(kk),100.0_wp*tv(kk)**pmean/max(1.0e-30_wp,ssum),frc(i,j)
    end do
    write (stdout,'(1x,a)') '--------------------------------------------------------------'
  end if

  effectivNat = av2*dble(n)
  deallocate (xyz,at,map,cn,b,sring,cring,wbo)
  deallocate (irev)   !> LOCAL EXTENSION (2026-09-29): reverse map
  deallocate (rho,nc)
  if (allocated(frc)) deallocate (frc)
  if (allocated(fcat)) deallocate (fcat)

  return
end subroutine flexi
subroutine minringsizes(nat,at,xyz,sring)
  use crest_parameters
  use zdata
  implicit none
  integer :: nat
  integer :: at(nat)
  real(wp) :: xyz(3,nat)
  integer :: sring(nat)
  type(zmolecule) :: zmol
  integer :: i,j

  call simpletopo(nat,at,xyz,zmol,.false.,.true.,'')
  do i = 1,zmol%nat
    sring(i) = 0
    do j = 1,zmol%nri
      if (any(zmol%zri(j)%rlist == i)) then
        if (sring(i) == 0) then
          sring(i) = zmol%zri(j)%rs
        else if (zmol%zri(j)%rs < sring(i)) then
          sring(i) = zmol%zri(j)%rs
        end if
      end if
    end do
  end do
  call zmol%deallocate()
  return
end subroutine minringsizes
!========================================================================================!

subroutine nciflexi_gfnff(mol,flexval)
!*************************************************************
!* new version of the non-covalent flexibility measure.
!* the nci-flexibility is estimated RELATIVE to a bioorganic
!* molecule (crambin) and its hydrogen-bond and dispersion
!* energy per atom.
!*************************************************************
  use crest_parameters
  use crest_calculator
  use strucrd
  use gfnff_api
  implicit none
  !> INPUT
  type(coord),intent(in) :: mol
  !> OUTPUT
  real(wp),intent(out)   :: flexval
  !> LOCAL
  integer :: io
  logical :: ex
  real(wp) :: ehb,edisp
  real(wp) :: energy
  real(wp),allocatable :: grad(:,:)
  type(gfnff_data),allocatable :: ff_dat
!>--- reset
  flexval = 0.0_wp
  ehb = 0.0_wp
  edisp = 0.0_wp

!>--- preprocessor statement (only executed when compiled with GFN-FF support)
#ifdef WITH_GFNFF
  allocate (ff_dat)
  allocate (grad(3,mol%nat),source=0.0_wp)

!>--- set up a GFN-FF calculation
  write (stdout,'(1x,a)',advance='no') 'Calculating NCI flexibility ...'
  flush (stdout)
  call gfnff_api_setup(mol,mol%chrg,ff_dat,io,.false.)
  if (io /= 0) then
    write (stdout,'(a)') ' failed.'
    deallocate (grad,ff_dat)
    return
  end if

!>--- perform a singlepoint calculation
  call gfnff_sp(mol,ff_dat,energy,grad,io)
  if (io /= 0) then
    write (stdout,'(a)') ' failed.'
    deallocate (grad,ff_dat)
    return
  end if
  ehb = ff_dat%res%e_hb
  edisp = ff_dat%res%e_disp
  write (stdout,'(a)') ' done.'

!>--- normalize by number of atoms
  ehb = ehb/mol%nat
  edisp = edisp/mol%nat

!>--- NCI flexi is determined RELATIVE to a reference molecule (Crambin)
  flexval = 0.5_wp*(1.0_wp-(ehb/(-0.00043374_wp)))
  flexval = flexval+0.5_wp*(1.0_wp-(edisp/(-0.00163029_wp)))

  deallocate (grad,ff_dat)
#else

  write (stdout,'(a)') 'Caclulation of NCI flexibility skipped.'
  write (stdout,'(a)') 'Set up compilation with -DWITH_GFNFF=true to enable this feature.'
  flexval = 0.0_wp

#endif

  return
end subroutine nciflexi_gfnff
!========================================================================================!
