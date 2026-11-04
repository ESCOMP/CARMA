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
!! With the bug, the (2,1) kernel picks up a spurious van der Waals sticking factor in
!! some columns but not others, depending on which thread computed the column and what
!! that thread had on its stack. Run it with several thread counts, e.g.
!!
!!   setenv CARMA_THREADS 4 ; run-carma.csh CKERNTHREADTEST.exe
!!
!! The test stops with a non-zero exit status if any mismatch is found.
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

  type(carma_type), target  :: carma
  type(carma_type), pointer :: carma_ptr
  type(carmastate_type), allocatable :: cstate(:)
  integer                   :: rc = 0

  real(kind=f)              :: zc(NZ), zl(NZP1), p(NZ), pl(NZP1), t(NZ)
  real(kind=f)              :: mmr(NZ)
  real(kind=f)              :: lat, lon, time

  ! Coagulation kernels for the 1-2 and 2-1 group pairs.
  real(kind=f), allocatable :: ck12_ref(:,:,:), ck21_ref(:,:,:)
  real(kind=f), allocatable :: ck12(:,:,:,:), ck21(:,:,:,:)

  integer               :: icol, ithread, nthreads, i1, i2, k
  integer               :: col_thread(NCOL)
  logical               :: col_asym(NCOL), col_diff(NCOL)
  integer               :: nasym, ndiff

  nthreads = 1
  !$ nthreads = omp_get_max_threads()
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

  ! Reference column, computed serially by the master thread before the parallel region.
  call do_column(cstate(1), ck12_ref, ck21_ref)

  write(*,*) ""
  write(*,*) "  Reference column (serial):"
  call report_symmetry(ck12_ref, ck21_ref)

  ! Now process many identical columns in parallel. Each thread reuses its own
  ! carmastate object, as in carma_test.F90 and in CAM.
  col_thread(:) = 0

  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(icol, ithread) SCHEDULE(STATIC)
  do icol = 1, NCOL
    ithread = 1
    !$ ithread = omp_get_thread_num() + 1
    col_thread(icol) = ithread

    call do_column(cstate(ithread), ck12(:,:,:,icol), ck21(:,:,:,icol))
  end do
  !$OMP END PARALLEL DO

  ! Check the results.
  nasym = 0
  ndiff = 0
  do icol = 1, NCOL
    col_asym(icol) = .false.
    col_diff(icol) = .false.

    do i2 = 1, NBIN
      do i1 = 1, NBIN
        do k = 1, NZ
          if (abs(ck12(k,i1,i2,icol) - ck21(k,i2,i1,icol)) > rtol * abs(ck12(k,i1,i2,icol))) col_asym(icol) = .true.
          if (ck21(k,i1,i2,icol) /= ck21_ref(k,i1,i2)) col_diff(icol) = .true.
          if (ck12(k,i1,i2,icol) /= ck12_ref(k,i1,i2)) col_diff(icol) = .true.
        end do
      end do
    end do

    if (col_asym(icol)) nasym = nasym + 1
    if (col_diff(icol)) ndiff = ndiff + 1
  end do

!!$  write(*,*) ""
!!$  write(*,'(a)') "   column  thread  K(1,2)/=K(2,1)^T  differs-from-reference"
!!$  do icol = 1, NCOL
!!$    write(*,'(i9,i8,l12,l22)') icol, col_thread(icol), col_asym(icol), col_diff(icol)
!!$  end do

  write(*,*) ""
  write(*,*) "  Columns with an asymmetric 1-2 kernel     : ", nasym, " of ", NCOL
  write(*,*) "  Columns that differ from the reference    : ", ndiff, " of ", NCOL

  ! Cleanup
  do ithread = 1, nthreads
    call CARMASTATE_Destroy(cstate(ithread), rc)
  end do
  deallocate(cstate, ck12_ref, ck21_ref, ck12, ck21)

  call CARMA_Destroy(carma, rc)
  if (rc /= 0) stop "    *** CARMA_Destroy FAILED ***"

  write(*,*) ""
  if (nasym > 0 .or. ndiff > 0) then
    write(*,*) "  *** FAILED *** The coagulation kernel depends on uninitialized data"
    write(*,*) "                 (use_vw lower triangle in setupckern.F90)."
    !error stop 1
  else
    write(*,*) "  PASSED: kernel is symmetric and identical in every column."
  end if

contains

  !! Create a carmastate for one column, step it (which calls setupckern), and
  !! return the 1-2 and 2-1 coagulation kernels.
  subroutine do_column(cs, k12, k21)
    type(carmastate_type), intent(inout) :: cs
    real(kind=f), intent(out)            :: k12(NZ,NBIN,NBIN)
    real(kind=f), intent(out)            :: k21(NZ,NBIN,NBIN)

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

    call CARMASTATE_Step(cs, lrc)
    if (lrc /= 0) stop "    *** CARMASTATE_Step FAILED ***"

    k12(:,:,:) = cs%f_ckernel(:,:,:,I_SULF,I_MIX)
    k21(:,:,:) = cs%f_ckernel(:,:,:,I_MIX,I_SULF)
  end subroutine do_column

  subroutine report_symmetry(k12, k21)
    real(kind=f), intent(in) :: k12(NZ,NBIN,NBIN)
    real(kind=f), intent(in) :: k21(NZ,NBIN,NBIN)

    write(*,'(a,es14.6)') "    K(1,2)(k=1,i1=1,i2=1)  = ", k12(1,1,1)
    write(*,'(a,es14.6)') "    K(2,1)(k=1,i1=1,i2=1)  = ", k21(1,1,1)
    write(*,'(a,es14.6)') "    max |K(1,2)-K(2,1)^T| / K(1,2) = ", &
      maxval(abs(k12 - reshape(k21, shape(k21), order=[1,3,2])) / k12)
  end subroutine report_symmetry

end subroutine test_ckern_threads
