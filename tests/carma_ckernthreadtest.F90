!! This test illustrates a thread-safety / uninitialized-memory bug in setupckern.F90.
!!
!! In setupckern, the van der Waals flag use_vw(j1,j2) is a local (stack) array that
!! was only filled for the upper triangle (j2 >= j1):
!!
!!     do j1 = 1, NGROUP
!!       do j2 = j1, NGROUP
!!         use_vw(j1, j2) = is_grp_sulfate(j1) .and. is_grp_sulfate(j2)
!!
!! but the kernel loop that follows visits every (j1,j2) pair with icoag(j1,j2) /= 0,
!! and setupcoag makes icoag symmetric. So for any two groups that coagulate with each
!! other, use_vw(j2,j1) (j2 > j1) is read without having been set. Its value is whatever
!! was left on the calling thread's stack, which depends on what that thread executed
!! before. In CAM this made the answers depend on the number of OpenMP threads (each
!! thread processes a different set of chunks), and an ERC test with threading failed.
!!
!! The test sets up two groups, like the CAM trop_strat models:
!!   group 1 : pure sulfate (is_sulfate = .true.)
!!   group 2 : mixed aerosol (not sulfate)
!! with coagulation 1+1->1, 1+2->2, 2+2->2. Since group 2 is not sulfate, the van der
!! Waals enhancement must not be applied to the 1-2 pair in either order.
!!
!! Many identical columns are processed in an OpenMP parallel loop (as in carma_test.F90).
!! For each column the coagulation kernel is checked for:
!!   1) symmetry : ckernel(k,i1,i2,1,2) must equal ckernel(k,i2,i1,2,1)
!!   2) reproducibility : every column must give the same kernel as the reference
!!      column computed serially before the parallel loop.
!!
!! Because the bad value is whatever happens to be on the stack, simply running the
!! model only exposes the bug some of the time. To make the failure deterministic, the
!! stack is "poisoned" before setupckern is called: a helper routine fills a large
!! region of the stack with a known bit pattern, and then setupckern is called directly
!! from the same routine, so its local arrays (including use_vw) occupy the poisoned
!! memory. Each column's kernel is computed twice, once with the stack filled with
!! zeros (reads as .false.) and once with all bits set (reads as .true. for any
!! logical convention). A correct setupckern gives identical kernels; with the bug the
!! all-ones pass applies the van der Waals factor to the (2,1) pair. Run it with
!! several thread counts, e.g.
!!
!!   setenv CARMA_THREADS 4 ; run-carma.csh CKERNTHREADTEST.exe
!!
!! The test stops with a non-zero exit status (error stop 1) if any mismatch is found.
!!
!! NOTE: The poisoning relies on setupckern's automatic arrays living on the stack,
!! which is the case for ifx/ifort (with -qopenmp or without -heap-arrays). Compilers
!! that put automatic arrays on the heap may still miss the bug. On the master thread
!! the OpenMP runtime overwrites the poisoned word on entry to setupckern, so the
!! parallel loop always uses at least 2 threads, even when CARMA_THREADS=1.
program carma_ckernthreadtest
  implicit none

  write(*,*) "Coagulation Kernel Threading Test (setupckern use_vw)"

  call test_ckern_threads()

  write(*,*) "Done"
end program


subroutine test_ckern_threads()
  use carma_precision_mod
  use carma_constants_mod
  use carma_enums_mod
  use carma_types_mod
  use carmaelement_mod
  use carmagroup_mod
  use carmastate_mod
  use carma_mod
  use atmosphere_mod
  !$ use omp_lib

  implicit none

  integer, parameter    :: NCOL         = 64
  integer, parameter    :: NZ           = 20
  integer, parameter    :: NZP1         = NZ+1
  integer, parameter    :: NELEM        = 2
  integer, parameter    :: NBIN         = 20
  integer, parameter    :: NGROUP       = 2
  integer, parameter    :: NSOLUTE      = 0
  integer, parameter    :: NGAS         = 0
  integer, parameter    :: NWAVE        = 0

  integer, parameter    :: I_SULF       = 1   ! group 1, pure sulfate
  integer, parameter    :: I_MIX        = 2   ! group 2, mixed (non-sulfate) aerosol
  integer, parameter    :: I_H2SO4      = 1   ! element composition

  real(kind=f), parameter   :: dtime  = 1800._f
  real(kind=f), parameter   :: deltaz = 1000._f
  real(kind=f), parameter   :: zmin   = 0._f
  real(kind=f), parameter   :: rtol   = 1.e-12_f

  ! Stack poisoning: number of 8-byte words to fill, and the two fill patterns.
  integer, parameter        :: NPOISON      = 32768        ! 256 KB
  integer(kind=8), parameter :: POISON_ZERO = 0_8
  integer(kind=8), parameter :: POISON_ONES = -1_8         ! all bits set

  type(carma_type), target  :: carma
  type(carma_type), pointer :: carma_ptr
  type(carmastate_type), allocatable :: cstate(:)
  integer                   :: rc = 0

  real(kind=f)              :: zc(NZ), zl(NZP1), p(NZ), pl(NZP1), t(NZ)
  real(kind=f)              :: mmr(NZ)
  real(kind=f)              :: lat, lon, time

  ! Coagulation kernels for the 1-2 and 2-1 group pairs.
  ! (ck*_ref: serial reference; ck*: zero-poisoned stack; cp*: ones-poisoned stack)
  real(kind=f), allocatable :: ck12_ref(:,:,:), ck21_ref(:,:,:)
  real(kind=f), allocatable :: ck12(:,:,:,:), ck21(:,:,:,:)
  real(kind=f), allocatable :: cp12(:,:,:,:), cp21(:,:,:,:)

  integer               :: icol, ithread, nthreads, i1, i2, k
  integer               :: col_thread(NCOL)
  logical               :: col_asym(NCOL), col_diff(NCOL), col_pois(NCOL)
  integer               :: nasym, ndiff, npois

  ! Always use at least 2 threads in the parallel loop. On the master thread, the
  ! OpenMP runtime call that ifx inserts at the start of setupckern (it contains a
  ! CRITICAL section) overwrites part of the poisoned stack, so only columns computed
  ! by worker threads are guaranteed to see the poison pattern.
  nthreads = 1
  !$ nthreads = max(2, omp_get_max_threads())
  write(*,*) "  Number of OpenMP threads : ", nthreads
  write(*,*) "  Number of columns        : ", NCOL

  ! Define the CARMA model
  call CARMA_Create(carma, NBIN, NELEM, NGROUP, NSOLUTE, NGAS, NWAVE, rc, LUNOPRT=6)
  if (rc /= 0) stop "    *** CARMA_Create FAILED ***"
  carma_ptr => carma

  call CARMAGROUP_Create(carma, I_SULF, "pure sulfate", 2.e-8_f, 2._f, I_SPHERE, 1._f, .false., rc, &
                         shortname="PRSUL", is_sulfate=.true.)
  if (rc /= 0) stop "    *** CARMAGROUP_Create FAILED ***"

  call CARMAGROUP_Create(carma, I_MIX, "mixed aerosol", 2.e-8_f, 2._f, I_SPHERE, 1._f, .false., rc, &
                         shortname="MXAER", is_sulfate=.false.)
  if (rc /= 0) stop "    *** CARMAGROUP_Create FAILED ***"

  call CARMAELEMENT_Create(carma, 1, I_SULF, "pure sulfate", 1.923_f, I_INVOLATILE, I_H2SO4, rc, shortname="PRSUL")
  if (rc /= 0) stop "    *** CARMAELEMENT_Create FAILED ***"

  call CARMAELEMENT_Create(carma, 2, I_MIX, "mixed aerosol", 1.923_f, I_INVOLATILE, I_H2SO4, rc, shortname="MXAER")
  if (rc /= 0) stop "    *** CARMAELEMENT_Create FAILED ***"

  call CARMA_AddCoagulation(carma, I_SULF, I_SULF, I_SULF, I_COLLEC_DATA, rc)
  if (rc /= 0) stop "    *** CARMA_AddCoagulation FAILED ***"
  call CARMA_AddCoagulation(carma, I_SULF, I_MIX,  I_MIX,  I_COLLEC_DATA, rc)
  if (rc /= 0) stop "    *** CARMA_AddCoagulation FAILED ***"
  call CARMA_AddCoagulation(carma, I_MIX,  I_MIX,  I_MIX,  I_COLLEC_DATA, rc)
  if (rc /= 0) stop "    *** CARMA_AddCoagulation FAILED ***"

  call CARMA_Initialize(carma, rc, do_coag=.true.)
  if (rc /= 0) stop "    *** CARMA_Initialize FAILED ***"

  ! A simple Cartesian column; every column in the test is identical.
  lat = 40._f
  lon = -105._f
  time = 0._f

  do k = 1, NZ
    zc(k) = zmin + (deltaz * (k - 0.5_f))
  end do
  call GetStandardAtmosphere(zc, p=p, t=t)

  do k = 1, NZP1
    zl(k) = zmin + ((k - 1) * deltaz)
  end do
  call GetStandardAtmosphere(zl, p=pl)

  mmr(:) = 1.e-10_f

  allocate(cstate(nthreads))
  allocate(ck12_ref(NZ,NBIN,NBIN), ck21_ref(NZ,NBIN,NBIN))
  allocate(ck12(NZ,NBIN,NBIN,NCOL), ck21(NZ,NBIN,NBIN,NCOL))
  allocate(cp12(NZ,NBIN,NBIN,NCOL), cp21(NZ,NBIN,NBIN,NCOL))

  ! Reference column, computed serially by the master thread before the parallel region
  ! (with a zero-poisoned stack, i.e. the "benign" garbage).
  call do_column(cstate(1), ck12_ref, ck21_ref, cp12(:,:,:,1), cp21(:,:,:,1))

  write(*,*) ""
  write(*,*) "  Reference column (serial):"
  call report_symmetry(ck12_ref, ck21_ref)

  ! Now process many identical columns in parallel. Each thread reuses its own
  ! carmastate object, as in carma_test.F90 and in CAM.
  col_thread(:) = 0

  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(icol, ithread) SCHEDULE(STATIC) NUM_THREADS(nthreads)
  do icol = 1, NCOL
    ithread = 1
    !$ ithread = omp_get_thread_num() + 1
    col_thread(icol) = ithread

    call do_column(cstate(ithread), ck12(:,:,:,icol), ck21(:,:,:,icol), &
                   cp12(:,:,:,icol), cp21(:,:,:,icol))
  end do
  !$OMP END PARALLEL DO

  ! Check the results.
  nasym = 0
  ndiff = 0
  npois = 0
  do icol = 1, NCOL
    col_asym(icol) = .false.
    col_diff(icol) = .false.
    col_pois(icol) = .false.

    do i2 = 1, NBIN
      do i1 = 1, NBIN
        do k = 1, NZ
          if (abs(ck12(k,i1,i2,icol) - ck21(k,i2,i1,icol)) > rtol * abs(ck12(k,i1,i2,icol))) col_asym(icol) = .true.
          if (abs(cp12(k,i1,i2,icol) - cp21(k,i2,i1,icol)) > rtol * abs(cp12(k,i1,i2,icol))) col_asym(icol) = .true.
          if (cp12(k,i1,i2,icol) /= ck12(k,i1,i2,icol)) col_pois(icol) = .true.
          if (cp21(k,i1,i2,icol) /= ck21(k,i1,i2,icol)) col_pois(icol) = .true.
          if (ck21(k,i1,i2,icol) /= ck21_ref(k,i1,i2)) col_diff(icol) = .true.
          if (ck12(k,i1,i2,icol) /= ck12_ref(k,i1,i2)) col_diff(icol) = .true.
        end do
      end do
    end do

    if (col_asym(icol)) nasym = nasym + 1
    if (col_diff(icol)) ndiff = ndiff + 1
    if (col_pois(icol)) npois = npois + 1
  end do

!!$  write(*,*) ""
!!$  write(*,'(a)') "   column  thread  K(1,2)/=K(2,1)^T  differs-from-reference"
!!$  do icol = 1, NCOL
!!$    write(*,'(i9,i8,l12,l22)') icol, col_thread(icol), col_asym(icol), col_diff(icol)
!!$  end do

  write(*,*) ""
  write(*,*) "  Columns with an asymmetric 1-2 kernel     : ", nasym, " of ", NCOL
  write(*,*) "  Columns that differ from the reference    : ", ndiff, " of ", NCOL
  write(*,*) "  Columns that depend on stack contents     : ", npois, " of ", NCOL

  ! Cleanup
  do ithread = 1, nthreads
    call CARMASTATE_Destroy(cstate(ithread), rc)
  end do
  deallocate(cstate, ck12_ref, ck21_ref, ck12, ck21, cp12, cp21)

  call CARMA_Destroy(carma, rc)
  if (rc /= 0) stop "    *** CARMA_Destroy FAILED ***"

  write(*,*) ""
  if (nasym > 0 .or. ndiff > 0) then
    write(*,*) "  *** FAILED *** The coagulation kernel depends on uninitialized data"
    write(*,*) "                 (use_vw lower triangle in setupckern.F90)."
    error stop 1
  else
    write(*,*) "  PASSED: kernel is symmetric and identical in every column."
  end if

contains

  !! Create a carmastate for one column and step it, then recompute the coagulation
  !! kernel twice by calling setupckern directly, first with the stack filled with
  !! zeros and then with all bits set. Return the 1-2 and 2-1 kernels from both passes.
  subroutine do_column(cs, k12, k21, p12, p21)
    type(carmastate_type), intent(inout) :: cs
    real(kind=f), intent(out)            :: k12(NZ,NBIN,NBIN)   ! zero-poisoned
    real(kind=f), intent(out)            :: k21(NZ,NBIN,NBIN)
    real(kind=f), intent(out)            :: p12(NZ,NBIN,NBIN)   ! ones-poisoned
    real(kind=f), intent(out)            :: p21(NZ,NBIN,NBIN)

    integer :: lrc, ie, ib

    lrc = 0
    call CARMASTATE_Create(cs, carma_ptr, time, dtime, NZ, I_CART, lat, lon, &
                           zc, zl, p, pl, t, lrc)
    if (lrc /= 0) stop "    *** CARMASTATE_Create FAILED ***"

    do ie = 1, NELEM
      do ib = 1, NBIN
        call CARMASTATE_SetBin(cs, ie, ib, mmr, lrc)
      end do
    end do

    ! Step once so that the wet radii, etc. used by setupckern are set up.
    call CARMASTATE_Step(cs, lrc)
    if (lrc /= 0) stop "    *** CARMASTATE_Step FAILED ***"

    ! setupckern is called from this routine right after poison_stack, so its stack
    ! frame (and the automatic array use_vw) lands in the memory just poisoned.
    call poison_stack(POISON_ZERO)
    call setupckern(carma_ptr, cs, lrc)
    if (lrc /= 0) stop "    *** setupckern FAILED ***"
    k12(:,:,:) = cs%f_ckernel(:,:,:,I_SULF,I_MIX)
    k21(:,:,:) = cs%f_ckernel(:,:,:,I_MIX,I_SULF)

    call poison_stack(POISON_ONES)
    call setupckern(carma_ptr, cs, lrc)
    if (lrc /= 0) stop "    *** setupckern FAILED ***"
    p12(:,:,:) = cs%f_ckernel(:,:,:,I_SULF,I_MIX)
    p21(:,:,:) = cs%f_ckernel(:,:,:,I_MIX,I_SULF)
  end subroutine do_column

  !! Fill a large region of the current thread's stack with a bit pattern. The routine
  !! is recursive so that buf is a true stack variable (not SAVEd), and buf is volatile
  !! so the stores can not be optimized away.
  recursive subroutine poison_stack(pattern)
    integer(kind=8), intent(in) :: pattern

    integer(kind=8), volatile   :: buf(NPOISON)
    integer                     :: i

    do i = 1, NPOISON
      buf(i) = pattern
    end do
    if (buf(1) /= pattern) write(*,*) "poison_stack: unexpected value"
  end subroutine poison_stack

  subroutine report_symmetry(k12, k21)
    real(kind=f), intent(in) :: k12(NZ,NBIN,NBIN)
    real(kind=f), intent(in) :: k21(NZ,NBIN,NBIN)

    write(*,'(a,es14.6)') "    K(1,2)(k=1,i1=1,i2=1)  = ", k12(1,1,1)
    write(*,'(a,es14.6)') "    K(2,1)(k=1,i1=1,i2=1)  = ", k21(1,1,1)
    write(*,'(a,es14.6)') "    max |K(1,2)-K(2,1)^T| / K(1,2) = ", &
      maxval(abs(k12 - reshape(k21, shape(k21), order=[1,3,2])) / k12)
  end subroutine report_symmetry

end subroutine test_ckern_threads
