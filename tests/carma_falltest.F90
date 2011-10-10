!! This code is to demonstrate the CARMA sedimentation routines
!! using a constant fall velocity.
!!
!! Upon execution, a text file (carma_falltest.txt) is generated.
!! The text file can be read with the IDL procedure read_falltest.pro.
!!
!! @author Peter Colarco (based on Chuck Bardeen's code)
!! @version Feb-2009

program carma_falltest
  implicit none

  write(*,*) "Sedimentation Test"

  call test_sedimentation()  
  
  write(*,*) "Done"
end program


subroutine test_sedimentation()
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

  integer, parameter    :: NX           = 1
  integer, parameter    :: NY           = 1
  integer, parameter    :: NZ           = 110
  integer, parameter    :: NZP1         = NZ+1
  integer, parameter    :: NELEM        = 1
  integer, parameter    :: NBIN         = 8
  integer, parameter    :: NGROUP       = 1
  integer, parameter    :: NSOLUTE      = 0
  integer, parameter    :: NGAS         = 0
  integer, parameter    :: NWAVE        = 0
  integer, parameter    :: nstep        = 100*6

  real(kind=f), parameter   :: dtime  = 1000._f
  real(kind=f), parameter   :: deltax = 100._f
  real(kind=f), parameter   :: deltay = 100._f
  real(kind=f), parameter   :: deltaz = 100._f
  real(kind=f), parameter   :: zmin   = 0._f
  
  integer, parameter        :: I_DUST       = 1
  integer, parameter        :: I_ICE        = 2

  type(carma_type), target            :: carma
  type(carma_type), pointer           :: carma_ptr
  type(carmastate_type)               :: cstate
  integer                             :: rc = 0
  
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
  integer, parameter    :: lun = 42

  real(kind=f)          :: time
  real(kind=f)          :: rmin, rmrat, rho
  logical               :: do_explised = .false.
!  logical               :: do_explised = .true.
  real(kind=f)          :: vf_const = 2.0_f
!  real(kind=f)          :: vf_const = 0.0_f
    

!  write(*,*) ""
!  write(*,*) "Sedimentation of Dust Particles"

  ! Open the output text file
  open(unit=lun,file="carma_falltest.txt",status="unknown")
  
  ! Allocate the arrays that we need for the model
  allocate(xc(NZ,NY,NX), dx(NZ,NY,NX), yc(NZ,NY,NX), dy(NZ,NY,NX), &
           zc(NZ,NY,NX), zl(NZP1,NY,NX), p(NZ,NY,NX), pl(NZP1,NY,NX), &
           t(NZ,NY,NX))
  allocate(mmr(NZ,NY,NX,NELEM,NBIN))
  allocate(lat(NY,NX), lon(NY,NX))  

  ! Define the particle-grid extent of the CARMA test
!  write(*,*) "  CARMA_Create(carma, ", NBIN,    ", ", NELEM, ", ", NGROUP, &
!                                 ", ", NSOLUTE, ", ", NGAS, ", rc, 6) ..."
  call CARMA_Create(carma, NBIN, NELEM, NGROUP, NSOLUTE, NGAS, NWAVE, rc, LUNOPRT=6)
  if (rc /=0) write(*, *) "    *** FAILED ***, rc=", rc
	carma_ptr => carma


  ! Define the group
  rho = 2.65_f
  rmrat = 2.0
  rmin = 7.5e-4_f
!  write(*,*) "  Add Group(s) ..."
  call CARMAGROUP_Create(carma, 1, "dust", rmin, rmrat, I_SPHERE, 1._f, .FALSE., rc)
  if (rc /=0) write(*, *) "    *** FAILED ***, rc=", rc

  
  ! Define the element
!  write(*,*) "  Add Element(s) ..."
  call CARMAELEMENT_Create(carma, 1, 1, "dust", rho, I_INVOLATILE, I_DUST, rc)
  if (rc /=0) write(*, *) "    *** FAILED ***, rc=", rc
  
  
!  write(*,*) "  CARMA_AddCoagulation(carma, 1, 1, 1, I_COLLEC_DATA, rc) ..."
!  call CARMA_AddCoagulation(carma, 1, 1, 1, I_COLLEC_DATA, rc) 
!  if (rc /=0) write(*, *) "    *** FAILED ***, rc=", rc
  
! Setup the CARMA processes to exercise
!  do_explised = .true.
!  vf_const = 2.0
!  write(*,*) "  CARMA_Initialize(carma, rc, do_vtran=.TRUE., "// &
!               "vf_const=", vf_const,", do_explised=",do_explised,") ..."
  call CARMA_Initialize(carma, rc, do_vtran=.TRUE., vf_const=vf_const, do_explised=do_explised)
  if (rc /=0) write(*, *) "    *** FAILED ***, rc=", rc
  
  ! Print the Group Information
!  write(*,*)  ""
!  call dumpGroup(carma, rc)
!  if (rc /=0) write(*, *) "    *** FAILED ***, rc=", rc
  
  ! Print the Element Information
!  write(*,*)  ""
!  call dumpElement(carma, rc)
!  if (rc /=0) write(*, *) "    *** FAILED ***, rc=", rc

!  write(*,*) ""
  
  ! For simplicity of setup, do a case with Cartesian coordinates,
  ! which are specified in this interface in meters.
  !
  ! NOTE: For Cartesian coordinates, the first level is the bottom 
  ! of the model (e.g. z = 0), while for sigma and hybrid coordinates
  ! the first level is the top of the model.
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

!  write(*,'(a6, 3a12)') "level", "zc", "p", "t"
!  write(*,'(a6, 3a12)') "", "(m)", "(Pa)", "(K)"
!  do i = 1, NZ
!    write(*,'(i6,3f12.3)') i, zc(i,NY,NX), p(i,NY,NX), t(i,NY,NX)
!  end do


  ! Vertical edge
  do i = 1, NZP1
    zl(i,:,:) = zmin + ((i - 1) * deltaz)
  end do
  call GetStandardAtmosphere(zl, p=pl)

!  write(*,*) ""
!  write(*,'(a6, 2a12)') "level", "zl", "pl"
!  write(*,'(a6, 2a12)') "", "(m)", "(Pa)"
!  do i = 1, NZP1
!    write(*,'(i6,2f12.3)') i, zl(i,NY,NX), pl(i,NY,NX)
!  end do

  			
  ! Put a blob in the model first bin at 8 km
  mmr(:,:,:,:,:) = 0._f
  do i = 1, NZ
   mmr(i,:,:,1,1) = 1e-10_f * exp( - ( ( zc(i,:,:) - 8.e3_f)/3.e3_f) ** 2) / &
                     ( p(i,:,:) / 287._f / t(i,:,:))
  end do

!  write(*,*)  ""
!  write(*, '(a6, 4a12)') "level", "mmr(i,NY,NX,1)", "mmr(i,NY,NX,2)"
!  do i = 1, NZ
!	  write(*, '(i6, 4g12.3)') i, mmr(i,NY,NX,1,1), mmr(i,NY,NX,1,2)
!  end do
  
  ! Write output for the falltest
  write(lun,*) NZ
  do i = 1, NZ
   write(lun,'(i3,2f10.1)') &
    i, zc(i,NY,NX), zl(i+1,NY,NX)-zl(i,NY,NX)
  end do
  
  write(lun,*) 0
  do i = 1, NZ
   write(lun,'(i3,e10.3,e10.3)') &
    i, real(mmr(i,NY,NX,1,1)), real(mmr(i,NY,NX,1,1)*p(i,NY,NX) / 287._f / t(i,NY,NX))
  end do


  ! Iterate the model over a few time steps.
 ! write(*,*) ""
  do istep = 1, nstep
  
    ! Calculate the model time.
    time = (istep - 1) * dtime

    ! NOTE: This means that there should not be any looping over NX or NY done
    ! in any other CARMA routines. They should only loop over NZ.
    do ixy = 1, NX*NY
       ix = ((ixy-1) / NY) + 1
       iy = ixy - (ix-1)*NY
			
       ! Create a CARMASTATE for this column.
       call CARMASTATE_Create(cstate, carma_ptr, time, dtime, NZ, &
                              I_CART, I_CART, lat(iy,ix), lon(iy,ix), &
                              xc(:,iy,ix), dx(:,iy,ix), &
                              yc(:,iy,ix), dy(:,iy,ix), &
                              zc(:,iy,ix), zl(:,iy,ix), p(:,iy,ix), &
                              pl(:,iy,ix), t(:,iy,ix), rc)
		
       ! Send the bin mmrs to CARMA
       do ielem = 1, NELEM
        do ibin = 1, NBIN
         call CARMASTATE_SetBin(cstate, ielem, ibin, &
                                mmr(:,iy,ix,ielem,ibin), rc)
        end do
       end do
			
       ! Execute the step
       call CARMASTATE_Step(cstate, rc)
       ! Get the updated bin mmr.
       do ielem = 1, NELEM
        do ibin = 1, NBIN
         call CARMASTATE_GetBin(cstate, ielem, ibin, &
                                mmr(:,iy,ix,ielem,ibin), rc)
        end do
       end do

       ! Get the updated temperature.
       call CARMASTATE_GetState(cstate, rc, t=t(:,iy,ix))
    enddo

    ! Write output for the falltest
    write(lun,'(f12.0)') istep*dtime
    do i = 1, NZ
      write(lun,'(i3,e10.3,e10.3)') &
        i, real(mmr(i,NY,NX,1,1)), real(mmr(i,NY,NX,1,1)*p(i,NY,NX) / 287._f / t(i,NY,NX))
    end do

  end do   ! time loop
	
  ! Cleanup the carma state objects
  call CARMASTATE_Destroy(cstate, rc)
  if (rc /=0) stop "    *** FAILED ***"

  ! Close the output file
  close(unit=lun)	
	
!  write(*,*)  ""
!  write(*,*)  ""
!  write(*, '(a8, 8a14)') "level", "mmr(i,NY,NX,1)", "mmr(i,NY,NX,2)"
!  do i = 1, NZ
!   write(*, '(i8, 8g14.3)') i, mmr(i,NY,NX,1,1), mmr(i,NY,NX,1,2)
!  end do
		
!  write(*,*)  ""
!  write(*, '(a8, 2a12)') "level", "t(:,1,1)", "t(:,NY,NX)"		
!  do i = 1, NZ
!   write(*, '(i8, 2f12.3)') i, t(i,1,1), t(i,NY,NX)
!  end do

  if (rc /=0) write(*, *) "    *** FAILED ***, rc=", rc

!  write(*,*)  ""
!  write(*,*) "  CARMA_Destroy() ..."
  call CARMA_Destroy(carma, rc)
  if (rc /=0) write(*, *) "    *** FAILED ***, rc=", rc  
end subroutine
