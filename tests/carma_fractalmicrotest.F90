!! This code is to demonstrate the CARMA microphysical treatment
!! of fractal aggregates composed of identical spheres
!!
!! Upon execution, a text file (carma_fractalmicrotest.txt) is generated.
!! The text file can be read with the IDL procedure read_fractaltest.pro.
!!
!! @author Eric Wolf (based on Chuck Bardeen's code)
!! @version Feb-2013

program carma_fractalmicrotest
  implicit none

  write(*,*) "Fractal Test"

  call test_fractalmicro()

  write(*,*) "Done"
end program


subroutine test_fractalmicro()
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
  integer, parameter    :: NBIN         = 10
  integer, parameter    :: NGROUP       = 1
  integer, parameter    :: NSOLUTE      = 0
  integer, parameter    :: NGAS         = 0
  integer, parameter    :: NWAVE        = 0
  integer, parameter    :: nstep        = 1500
  integer, parameter    :: nrun         = 2

  real(kind=f), parameter   :: dtime  = 604800._f
  !real(kind=f), parameter   :: dtime  = 86400._f
  real(kind=f), parameter   :: deltaz = 1000._f
  real(kind=f), parameter   :: zmin   = 0._f

  integer, parameter        :: I_DUST       = 1

  type(carma_type), target            :: carma
  type(carma_type), pointer           :: carma_ptr
  type(carmastate_type)               :: cstate
  integer                             :: rc = 0

  real(kind=f), allocatable   :: zc(:,:,:)
  real(kind=f), allocatable   :: zl(:,:,:)
  real(kind=f), allocatable   :: p(:,:,:)
  real(kind=f), allocatable   :: pl(:,:,:)
  real(kind=f), allocatable   :: t(:,:,:)

  real(kind=f), allocatable, target  :: mmr(:,:,:,:,:)

  real(kind=f), allocatable          :: lat(:,:)
  real(kind=f), allocatable          :: lon(:,:)

  integer		:: z
  integer               :: i
  integer               :: ix
  integer               :: iy
  integer               :: ixy
  integer               :: istep
  integer               :: irun
  integer               :: ielem
  integer               :: ibin
  integer               :: ithread
  integer, parameter    :: lun = 42

  real(kind=f)          :: time
  real(kind=f)          :: rmin, rmrat, rho, rmon_in, rmon
  real(kind=f)		:: falpha_in
  real(kind=f), allocatable   :: df_in(:)
  real(kind=f), allocatable   :: df(:)
  real(kind=f), allocatable   :: nmon(:)
  real(kind=f), allocatable   :: rrat(:)
  real(kind=f), allocatable   :: rprat(:)
  real(kind=f), allocatable   :: r(:)
  real(kind=f), allocatable   :: dr(:)
  real(kind=f), allocatable   :: rmass(:)
  real(kind=f), allocatable   :: mmr_emis(:,:,:,:,:)
  real(kind=f), allocatable   :: mmr_out(:,:,:,:,:,:,:)
  real(kind=f), allocatable   :: rhoa(:,:,:)

  logical               :: do_explised = .false.


!  write(*,*) ""
!  write(*,*) "Dry Fractal Aggregates of Spheres"

  ! Open the output text file
  open(unit=lun,file="carma_fractalmicrotest.txt",status="unknown")

  ! Allocate the arrays that we need for the model
  allocate(zc(NZ,NY,NX), zl(NZP1,NY,NX), p(NZ,NY,NX), pl(NZP1,NY,NX), &
           t(NZ,NY,NX))
  allocate(mmr(NZ,NY,NX,NELEM,NBIN))
  allocate(mmr_out(nrun,nstep,NZ,NY,NX,NELEM,NBIN))
  allocate(mmr_emis(NZ,NY,NX,NELEM,NBIN))
  allocate(lat(NY,NX), lon(NY,NX))
  allocate(df_in(NBIN), df(NBIN), r(NBIN), rrat(NBIN), rprat(NBIN), rmass(NBIN), dr(NBIN), nmon(NBIN))
  allocate(rhoa(NZ,NY,NX))


  do irun = 1, nrun

  ! Define the particle-grid extent of the CARMA test
  !  write(*,*) "  CARMA_Create(carma, ", NBIN,    ", ", NELEM, ", ", NGROUP, &
  !                                 ", ", NSOLUTE, ", ", NGAS, ", rc, 6) ..."
  call CARMA_Create(carma, NBIN, NELEM, NGROUP, NSOLUTE, NGAS, NWAVE, rc, LUNOPRT=6)
  if (rc /=0) write(*, *) "    *** FAILED ***, rc=", rc
	carma_ptr => carma


  ! Define the group
  rho = 2.65_f
  rmrat = 2.5_f
  rmin = 1.0e-5_f
  df_in(:) = 2.0_f
  rmon_in = 1.0e-5_f
  falpha_in = 1.0_f

!  write(*,*) "  Add Group(s) ..."

  if (irun .EQ. 1) then
    ! fractal run
    call CARMAGROUP_Create(carma, 1, "fractal", rmin, rmrat, I_SPHERE, 1._f, .FALSE., rc, &
                         is_fractal=.TRUE., rmon=rmon_in, df=df_in, falpha=falpha_in)


    if (rc /=0) write(*, *) "    *** FAILED ***, rc=", rc
    ! Define the element
    ! write(*,*) "  Add Element(s) ..."
    call CARMAELEMENT_Create(carma, 1, 1, "fractal", rho, I_INVOLATILE, I_DUST, rc)
    if (rc /=0) write(*, *) "    *** FAILED ***, rc=", rc
  end if

  if (irun .EQ. 2) then
    !fractal run
    call CARMAGROUP_Create(carma, 1, "sphere", rmin, rmrat, I_SPHERE, 1._f, .FALSE., rc, &
                         is_fractal=.FALSE. )
    if (rc /=0) write(*, *) "    *** FAILED ***, rc=", rc
    ! Define the element
    ! write(*,*) "  Add Element(s) ..."
    call CARMAELEMENT_Create(carma, 1, 1, "sphere", rho, I_INVOLATILE, I_DUST, rc)
    if (rc /=0) write(*, *) "    *** FAILED ***, rc=", rc
  end if

  call CARMA_AddCoagulation(carma, 1, 1, 1, I_COLLEC_DATA, rc)
  if (rc /=0) write(*, *) "    *** FAILED ***, rc=", rc

! Setup the CARMA processes to exercise
!  do_explised = .true.
!  write(*,*) "  CARMA_Initialize(carma, rc, do_vtran=.TRUE., "// &
!               , do_explised=",do_explised,") ..."
  call CARMA_Initialize(carma, rc, do_vtran=.TRUE., do_coag=.TRUE., do_vdiff=.TRUE. )
  if (rc /=0) write(*, *) "    *** FAILED ***, rc=", rc

  call CARMAGROUP_Get(carma, 1, rc, r=r, rrat=rrat, rprat=rprat, rmon=rmon, df=df, dr=dr, rmass=rmass, nmon=nmon)
  ! Write output for the test

  if (irun .EQ. 1) then
    write(lun,'(i4, "  ", e12.5)') NBIN, RMON_in

    do ibin = 1, NBIN
      write(lun,'(i4,f4.1,e10.3, 2f10.5,e10.3)') ibin, df_in(ibin),r(ibin) , rrat(ibin), rprat(ibin), nmon(ibin)
    end do
  end if

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


  ! Put inital aerosol mass in the model first bin at 60 km
  mmr(:,:,:,:,:) = 0._f
  do i = 1, NZ
   mmr(i,:,:,1,1) = 1e-10_f * exp( - ( ( zc(i,:,:) - 60.e3_f)/3.e3_f) ** 2) / &
                     ( p(i,:,:) / 287._f / t(i,:,:))
  end do

!  write(*,*)  ""
!  write(*, '(a6, 4a12)') "level", "mmr(i,NY,NX,1)", "mmr(i,NY,NX,2)"
!  do i = 1, NZ
!	  write(*, '(i6, 4g12.3)') i, mmr(i,NY,NX,1,1), mmr(i,NY,NX,1,2)
!  end do


  if (irun .EQ. 1) then
    ! Write output for the fractalmicrotest
    write(lun,*) NZ
    do i = 1, NZ
      write(lun,'(i3,2f10.1)') &
      i, zc(i,NY,NX), zl(i+1,NY,NX)-zl(i,NY,NX)
    end do
  end if

  ! Initial particle statistics
  rhoa(:,:,:) = p(:,:,:)/287._f/t(:,:,:)
  !write(lun,*) 0
  !do z = 1, NZ
  ! do i = 1, NBIN
  !
  !  write(lun,'(2i3,e10.3,e10.3,e10.3,e10.3)') &
  !         z, i,real(mmr(z,NY,NX,1,i)*rhoa(1,NY,NX)/rmass(i)*1e-6_f*1e3_f), &
  !              real(mmr(z,NY,NX,1,i)*p(z,NY,NX) / 287._f / t(z,NY,NX)), &
  !              real(mmr(z,NY,NX,1,i)*rhoa(1,NY,NX)/rmass(i)*1e-6_f*1e3_f), &
  !              real(mmr(z,NY,NX,1,i)*p(z,NY,NX) / 287._f / t(z,NY,NX))
  ! end do
  !end do

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
                              I_CART, lat(iy,ix), lon(iy,ix), &
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

    do i=1, NBIN
      mmr_out(irun,istep,:,NY,NX,1,i) = mmr(:,NY,NX,1,i)
    end do

    ! Source aerosl for next timestep
    !do i = 1, NZ
    ! mmr_emis(i,:,:,1,1) = 1e-10_f * exp( - ( ( zc(i,:,:) - 60.e3_f)/3.e3_f) ** 2) / &
    !                  ( p(i,:,:) / 287._f / t(i,:,:))
    !end do
    !mmr = mmr + mmr_emis

  end do   ! time loop

  ! Cleanup the carma state objects
  call CARMASTATE_Destroy(cstate, rc)
  if (rc /=0) stop "    *** FAILED ***"


  end do   ! run loop

   ! Write output for the fractaltest
   do istep = 1, nstep
    write(lun,'(f12.0)') istep*dtime
     do z = 1, NZ
      !do i = 1, NBIN

       write(lun,'(i3,e10.3,e10.3,e10.3,e10.3)') &
           z, SUM(real(mmr_out(1,istep,z,NY,NX,1,:)*rhoa(z,NY,NX)/rmass(:)*1e-6_f*1e3_f)), &
                SUM(real(mmr_out(1,istep,z,NY,NX,1,:)*p(z,NY,NX) / 287._f / t(z,NY,NX))), &
                SUM(real(mmr_out(2,istep,z,NY,NX,1,:)*rhoa(z,NY,NX)/rmass(:)*1e-6_f*1e3_f)), &
                SUM(real(mmr_out(2,istep,z,NY,NX,1,:)*p(z,NY,NX) / 287._f / t(z,NY,NX)))
      !end do
    end do
   end do


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
