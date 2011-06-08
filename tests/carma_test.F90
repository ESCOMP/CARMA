!! This is not a complete unit test, but it does try to at least exercise some of the code
!! paths of the CARMA module interface. More complete and rigorous testing of the CARMA module
!! is performed with the SLOD.
!!
!! @author Chuck Bardeen
!! @version Feb-2009
!  use carma_precision_mod
program carma_test
  implicit none

  write(*,*) "Simple CARMA Interface Tester"

  call test_dust_1Dv2()  
end program


subroutine test_dust_1Dv2()
  use carma_precision_mod 
  use carma_constants_mod 
  use carma_enums_mod 
  use carma_types_mod 
  use carmaelement_mod
  use carmagroup_mod
  use carmastate_mod
  use carma_mod
  use atmosphere_mod

  implicit none

  integer, parameter    :: NX           = 15
  integer, parameter    :: NY           = 15
  integer, parameter    :: NZ           = 80
  integer, parameter    :: NZP1         = NZ+1
  integer, parameter    :: NELEM        = 1
  integer, parameter    :: NBIN         = 28
  integer, parameter    :: NGAS         = 0
  integer, parameter    :: NSOLUTE      = 0
  integer, parameter    :: NWAVE        = 0
  integer, parameter    :: nstep        = 30

  real(kind=f), parameter   :: dtime  = 1800._f
  real(kind=f), parameter   :: deltax = 100._f
  real(kind=f), parameter   :: deltay = 100._f
  real(kind=f), parameter   :: deltaz = 1000._f
  real(kind=f), parameter   :: zmin   = 0._f
  
  integer, parameter        :: I_DUST       = 1
  integer, parameter        :: I_ICE        = 2

  type(carma_type), target  :: carma
  type(carma_type), pointer :: carma_ptr
  type(carmastate_type), allocatable     :: cstate(:)
  integer                   :: rc = 0
  
  real(kind=f), allocatable   :: xc(:,:,:)
  real(kind=f), allocatable   :: dx(:,:,:)
  real(kind=f), allocatable   :: yc(:,:,:)
  real(kind=f), allocatable   :: dy(:,:,:)
  real(kind=f), allocatable   :: zc(:,:,:)
  real(kind=f), allocatable   :: zl(:,:,:)
  real(kind=f), allocatable   :: p(:,:,:)
  real(kind=f), allocatable   :: pl(:,:,:)
  real(kind=f), allocatable   :: t(:,:,:)
  
  real(kind=f), allocatable, target  :: mmr(:,:,:,:,:)
  
  real(kind=f), allocatable          :: lat(:,:)
  real(kind=f), allocatable          :: lon(:,:)
  
  integer               :: i
  integer               :: ix
  integer               :: iy
  integer               :: ixy
  integer               :: istep
  integer               :: ielem
  integer               :: ibin
  integer               :: ithread

  real(kind=f)          :: time
  
  integer               :: omp_get_num_threads, omp_get_max_threads, omp_get_thread_num
  

  write(*,*) ""
  write(*,*) "Dust Model, 3D"
  
  ! Allocate the arrays that we need for the model
  allocate(xc(NZ,NY,NX), dx(NZ,NY,NX), yc(NZ,NY,NX), dy(NZ,NY,NX), zc(NZ,NY,NX), &
    zl(NZP1,NY,NX), p(NZ,NY,NX), pl(NZP1,NY,NX), t(NZ,NY,NX))
  allocate(mmr(NZ,NY,NX,NELEM,NBIN))
  allocate(lat(NY,NX), lon(NY,NX))
  
  write(*,*) "  CARMA_Create(carma, ", NBIN, ", 1, 1, 0, 0, rc, 6) ..."
  call CARMA_Create(carma, NBIN, 1, 1, 0, 0, 0, rc, LUNOPRT=6)
  if (rc /=0) write(*, *) "    *** FAILED ***"
	carma_ptr => carma
  
  write(*,*) "  Add Group(s) ..."
  call CARMAGROUP_Create(carma, 1, "meteoric dust", 2e-8_f, 2.0_f, I_SPHERE, 1._f, .FALSE., rc)
  if (rc < 0) stop "    *** FAILED ***"
  
  write(*,*) "  Add Element(s) ..."
  call CARMAELEMENT_Create(carma, 1, 1, "meteoric dust", 2._f, I_INVOLATILE, I_DUST, rc)
  if (rc /=0) stop "    *** FAILED ***"
 
  write(*,*) "  CARMA_AddCoagulation(carma, 1, 1, 1, I_COLLEC_DATA, rc) ..."
  call CARMA_AddCoagulation(carma, 1, 1, 1, I_COLLEC_DATA, rc) 
  if (rc /=0) write(*, *) "    *** FAILED ***, rc=", rc
  
  write(*,*) "  CARMA_Initialize(carma, rc, do_substep=.TRUE., do_vtran=.TRUE., do_coag=.TRUE., vf_const=2._f) ..."
  call CARMA_Initialize(carma, rc, do_vtran=.TRUE., do_coag=.TRUE.)
  if (rc /=0) write(*, *) "    *** FAILED ***, rc=", rc
  
  write(*,*)  ""
  call dumpGroup(carma, rc)
  if (rc /=0) write(*, *) "    *** FAILED ***, rc=", rc
  
  write(*,*) ""
  
  ! For simplicity of setup, do a case with Cartesian coordinates, which are specified
  ! in this interface in meters.
  !
  ! NOTE: For Cartesian coordinates, the first level is the bottom of the model (e.g. z = 0),
  ! while for sigma and hybrid coordinates the first level is the top of the model.
  lat(:,:) = 40.0_f
  lon(:,:) = -105.0_f
  
  ! Horizonal centers
  do ix = 1, NX
    do iy = 1, NY
      dx(:,iy,ix) = deltax
      xc(:,iy,ix) = ix*dx(:,iy,ix) / 2._f
      dy(:,iy,ix) = deltay
      yc(:,iy,ix) = iy*dy(:,iy,ix) / 2._f
    end do
  end do
  
  ! Vertical center
  do i = 1, NZ
    zc(i,:,:) = zmin + (deltaz * (i - 0.5_f))
  end do
  call GetStandardAtmosphere(zc, p=p, t=t)

  write(*,'(a6, 3a12)') "level", "zc", "p", "t"
  write(*,'(a6, 3a12)') "", "(m)", "(Pa)", "(K)"
  do i = 1, NZ
    write(*,'(i6,3f12.3)') i, zc(i,1,1), p(i,1,1), t(i,1,1)
  end do


  ! Vertical edge
  do i = 1, NZP1
    zl(i,:,:) = zmin + ((i - 1) * deltaz)
  end do
  call GetStandardAtmosphere(zl, p=pl)

  write(*,*) ""
  write(*,'(a6, 2a12)') "level", "zl", "pl"
  write(*,'(a6, 2a12)') "", "(m)", "(Pa)"
  do i = 1, NZP1
    write(*,'(i6,2f12.3)') i, zl(i,1,1), pl(i,1,1)
  end do

  			
  ! Put a blob at the top of the model in the first bin.
  mmr(:,:,:,:,:) = 0._f
  mmr(1,:,:,1,1) = 1e-6_f
  mmr(NZ,:,:,1,1) = 1e-6_f

  write(*,*)  ""
  write(*, '(a6, 4a12)') "level", "mmr(1,1,i,1)", "mmr(1,1,i,2)", "mmr(1,2,i,1)", "mmr(1,2,i,2)"
  do i = 1, NZ
	  write(*, '(i6, 4g12.3)') i, mmr(i,1,1,1,1), mmr(i,1,1,1,2), mmr(i,2,1,1,1), mmr(i,2,1,1,2)
  end do
  
  ! Allocate enough carmastate objects for the maximum number of threads.
  allocate(cstate(omp_get_max_threads()))
  
		
  ! Iterate the model over a few time steps.
  write(*,*) ""
  do istep = 1, nstep
  
    ! Calculate the model time.
    time = (istep - 1) * dtime
    write(*,'(a,i6,a,f8.2,a)') "  step ", istep, ", time=", time / 3600._f, " (hr)"
!    write(*,*) ""

		! NOTE: This means that there should not be any looping over NX or NY done
		! in any other CARMA routines. They should only loop over NZ.
		!
		! NOTE: This directive allows each column of the model to be processed in a
		! separate thread. This can allow for faster computation on machines that 
		! allow multiple threads (e.g. have multiple CPUS). This should probably not
		! be used when the the model is embedded in a another model that is already
		! controlling the distribution of the model across multiple threads.
		
		!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ixy,ix,iy,ielem,ibin,ithread)
		do ixy = 1, NX*NY
			ix = ((ixy-1) / NY) + 1
			iy = ixy - (ix-1)*NY
			
			! Use the thread number to determine which member of the cstate pool to use.
			ithread = omp_get_thread_num() + 1
			
!			write(*,"(i2,': (',i5,',',i5,')')") ithread, iy, ix
  
			! Create a CARMASTATE for this column.
	!    write(*,*) "  CARMASTATE_Create"
			call CARMASTATE_Create(cstate(ithread), carma_ptr, time, dtime, NZ, I_CART, I_CART, lat(iy,ix), lon(iy,ix), &
				xc(:,iy,ix), dx(:,iy,ix), yc(:,iy,ix), dy(:,iy,ix), zc(:,iy,ix), zl(:,iy,ix), p(:,iy,ix), &
				pl(:,iy,ix), t(:,iy,ix), rc)
		
			! Send the bin mmrs to CARMA
			do ielem = 1, NELEM
			  do ibin = 1, NBIN
				  call CARMASTATE_SetBin(cstate(ithread), ielem, ibin, mmr(:,iy,ix,ielem,ibin), rc)
				end do
			end do
			
			! Execute the step
	!    write(*,*)  ""
	!    write(*,*) "  CARMA_Step"
	!    write(*,*)  ""
			call CARMASTATE_Step(cstate(ithread), rc)
			
			! Get the updated bin mmr.
			do ielem = 1, NELEM
			  do ibin = 1, NBIN
	!			  write (*, *) "    CARMASTATE_GetBin() ..."
				  call CARMASTATE_GetBin(cstate(ithread), ielem, ibin, mmr(:,iy,ix,ielem,ibin), rc)
				end do
			end do
	
			! Get the updated temperature.
	!    write(*,*) "    CARMASTATE_GetState"
			call CARMASTATE_GetState(cstate(ithread), rc, t=t(:,iy,ix))
		enddo
		!$OMP END PARALLEL DO
		
	end do
	
	! Cleanup the carma state objects
	do i = 1, omp_get_max_threads()
	  call CARMASTATE_Destroy(cstate(i), rc)
	end do
	deallocate(cstate)
	
	
	write(*,*)  ""
	write(*,*)  ""
	write(*, '(a8, 8a14)') "level", "mmr(1,1,i,1)", "mmr(1,1,i,2)", "mmr(1,1,i,27)", "mmr(1,1,i,28)", "mmr(1,2,i,1)", "mmr(1,2,i,2)", "mmr(1,2,i,27)", "mmr(1,2,i,28)"
	do i = 1, NZ
		write(*, '(i8, 8g14.3)') i, mmr(i,1,1,1,1), mmr(i,1,1,1,2), mmr(i,1,1,1,27), mmr(i,1,1,1,28), mmr(i,2,1,1,1), mmr(i,2,1,1,2), mmr(i,2,1,1,27), mmr(i,2,1,1,28)
	end do
		
  write(*,*)  ""
	write(*, '(a8, 2a12)') "level", "t(1,1,:)", "t(1,2,:)"		
	do i = 1, NZ
		write(*, '(i8, 2f12.3)') i, t(i,1,1), t(i,2,1)
	end do

  if (rc /=0) write(*, *) "    *** FAILED ***, rc=", rc

  write(*,*)  ""
  write(*,*) "  CARMA_Destroy() ..."
  call CARMA_Destroy(carma, rc)
  if (rc /=0) write(*, *) "    *** FAILED ***, rc=", rc  
end subroutine

