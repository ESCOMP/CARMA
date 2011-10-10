!! This code is to test the dry deposition routine by comparing
!! sedimentation with and without dry deposition.
!!
!! Upon execution, a text file (carma_drydeptest.txt) is generated.
!! The text file can be read with the IDL procedure read_drydeptest.pro.
!!
!! @author  Tianyi Fan
!! @version Apr-2011

program carma_drydeptest
  implicit none

  write(*,*) "Dry Deposition Test"

  call test_drydep()  
  
  write(*,*) "Done"
end program

!! Create 2 particle groups, one for particles with dry deposition 
!! using arbitary values of ram(aerodynamic resistance) and fv (friction velocity)
!! the other for particles without dry deposition, to see if they make a difference.
subroutine test_drydep()
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

  integer, parameter        :: NZ           = 150
  integer, parameter        :: NZP1         = NZ+1
  integer, parameter        :: NELEM        = 2
  integer, parameter        :: NBIN         = 16
  integer, parameter        :: NGROUP       = 2
  integer, parameter        :: NSOLUTE      = 0
  integer, parameter        :: NGAS         = 0
  integer, parameter        :: NWAVE        = 0
  integer, parameter        :: LUNOPRT      = 6
  integer, parameter        :: nstep        = 100*6
  
  ! To keep the file processing simpler, only one bin will get written out
  ! to the output file. 
  integer, parameter        :: OUTBIN       = 14


  real(kind=f), parameter   :: dtime  = 1000._f
  real(kind=f), parameter   :: deltax = 100._f
  real(kind=f), parameter   :: deltay = 100._f
  real(kind=f), parameter   :: deltaz = 100._f
  real(kind=f), parameter   :: rhmin  = .4_f
  real(kind=f), parameter   :: rhmax  = 1.05_f
  real(kind=f), parameter   :: zmin   = 0._f
  
  integer, parameter        :: I_PART  = 1

  type(carma_type), target  :: carma
  type(carma_type), pointer :: carma_ptr
  type(carmastate_type)     :: cstate
  integer                   :: rc = 0
  
  real(kind=f), allocatable   :: xc(:)
  real(kind=f), allocatable   :: dx(:)
  real(kind=f), allocatable   :: yc(:)
  real(kind=f), allocatable   :: dy(:)
  real(kind=f), allocatable   :: zc(:)
  real(kind=f), allocatable   :: zl(:)
  real(kind=f), allocatable   :: p(:)
  real(kind=f), allocatable   :: pl(:)
  real(kind=f), allocatable   :: t(:)
  real(kind=f), allocatable   :: relhum(:)
  
  real(kind=f), allocatable, target  :: mmr(:,:,:)
  
  real(kind=f)                :: lat
  real(kind=f)                :: lon
  
  integer               :: i
  integer               :: istep
  integer               :: ielem
  integer               :: ibin
  integer               :: igroup
  integer, parameter    :: lun = 42
  integer, parameter    :: lun1 = 41
  
  real(kind=f)          :: time
  real(kind=f)          :: rmin, rmrat, rho
  real(kind=f)          :: drh
  
  real(kind=f)          :: lndfv    = 1.5_f   ! land friction velocity
  real(kind=f)          :: lndram   = 60._f   ! land aerodynamic resistance
  real(kind=f)          :: lndfrac  = 0.0_f   ! land fraction
  
  real(kind=f)          :: ocnfv    = 2.0_f   ! ocean friction velocity
  real(kind=f)          :: ocnram   = 40._f   ! ocean aerodynamic resistance
  real(kind=f)          :: ocnfrac  = 1.0_f   ! ocean fraction
  
  real(kind=f)          :: icefv    = 2.5_f   ! ice friction velocity
  real(kind=f)          :: iceram   = 20._f   ! ice aerodynamic resistance
  real(kind=f)          :: icefrac  = 0.0_f   ! ice fraction

  real(kind=f)          :: vdry(NBIN, NGROUP)
  

  ! Open the output text file
  open(unit=lun,file="carma_drydeptest.txt",status="unknown") 
  
  ! Allocate the arrays that we need for the model
  allocate(xc(NZ), dx(NZ), yc(NZ), dy(NZ), &
           zc(NZ), zl(NZP1), p(NZ), pl(NZP1), &
           t(NZ))
  allocate(mmr(NZ,NELEM,NBIN))
  allocate(relhum(NZ))

  ! Define the particle-grid extent of the CARMA test
  call CARMA_Create(carma, NBIN, NELEM, NGROUP, NSOLUTE, NGAS, NWAVE, rc, LUNOPRT=LUNOPRT)
  if (rc /=0) stop "    *** CARMA_Create FAILED ***"
  carma_ptr => carma


  ! Define the groups
  rho = 2.65_f
  rmrat = 4.32_f
  rmin  = 1e-6_f
  
  call CARMAGROUP_Create(carma, 1, "DryDep", rmin, rmrat, I_SPHERE, 1._f, .FALSE., rc, &
       do_mie = .FALSE., do_wetdep=.FALSE., do_drydep=.TRUE., do_vtran=.TRUE., shortname="DD") 
  if (rc /=0) stop "    *** CARMAGROUP_Create FAILED ***"

  call CARMAGROUP_Create(carma, 2, "NoDryDep", rmin, rmrat, I_SPHERE, 1._f, .FALSE., rc, &
       do_mie = .FALSE., do_wetdep=.FALSE., do_drydep=.FALSE., do_vtran=.TRUE., shortname="NDD")
  if (rc /=0) stop "    *** CARMAGROUP_Create FAILED ***"
  
  ! Define the elements
  call CARMAELEMENT_Create(carma, 1, 1, "DryDep", rho, I_INVOLATILE, I_PART, rc)
  if (rc /=0) stop "    *** CARMAELEMENT_Create FAILED ***"

  call CARMAELEMENT_Create(carma, 2, 2, "NoDryDep", rho, I_INVOLATILE, I_PART, rc)
  if (rc /=0) stop "    *** CARMAELEMENT_Create FAILED ***"
  
!  write(*,*) "  CARMA_AddCoagulation(carma, 1, 1, 1, I_COLLEC_DATA, rc) ..."
!  call CARMA_AddCoagulation(carma, 1, 1, 1, I_COLLEC_DATA, rc) 
!  if (rc /=0) stop "    *** CARMA_AddCoagulation FAILED ***"
  
  ! Setup the CARMA processes to exercise
  call CARMA_Initialize(carma, rc, do_vtran=.TRUE., do_drydep=.TRUE.)
  if (rc /=0) stop "    *** CARMA_Initialize FAILED ***"
  
  
  ! For simplicity of setup, do a case with Cartesian coordinates,
  ! which are specified in this interface in meters.
  !
  ! NOTE: For Cartesian coordinates, the first level is the bottom 
  ! of the model (e.g. z = 0), while for sigma and hybrid coordinates
  ! the first level is the top of the model.
  lat = -40.0_f
  lon = -105.0_f
  
  ! Horizonal centers
  dx(:) = deltax
  xc(:) = dx(:) / 2._f
  dy(:) = deltay
  yc(:) = dy(:) / 2._f
  
  ! Vertical center
  do i = 1, NZ
    zc(i) = zmin + (deltaz * (i - 0.5_f))
  end do
  call GetStandardAtmosphere(zc, p=p, t=t)

  ! Vertical edge
  do i = 1, NZP1
    zl(i) = zmin + ((i - 1) * deltaz)
  end do
  call GetStandardAtmosphere(zl, p=pl)

  
  ! Set up some arbitrary relative humidities, with the maximum at the bottom and
  ! minimum at the top. Make the RH at the top 0, just to make sure the code can
  ! handle an RH of 0.
  drh    = (rhmax - rhmin) / (NZ-1)
  
  do i = 1, NZ
    relhum(i) = rhmax - ((i - 1) * drh)
  end do
  
  relhum(NZ) = 0._f

        
  ! Put a blob in the model for all elements and bins at 8 km.
  mmr(:,:,:) = 0._f
  do i = 1, NZ
    do ielem = 1, NELEM
      do ibin = 1, NBIN
        mmr(i,ielem,ibin) = 1e-10_f * exp(-((zc(i) - 8.e3_f) / 3.e3_f)**2) / &
                            (p(i) / 287._f / t(i))
      end do
    end do
  end do

  
  ! Write output for the falltest
  write(lun,*) NZ, NELEM
  do i = 1, NZ
   write(lun,'(i3,2f10.1)') i, zc(i), zl(i+1)-zl(i)
  end do
  
  write(lun,*) 0
  do ielem = 1, NELEM
    do i = 1, NZ
      write(lun,'(2i4,e10.3,e10.3)') ielem, i, real(mmr(i,ielem,OUTBIN)), real(mmr(i,ielem,OUTBIN)*p(i) / 287._f / t(i))
    end do
  end do
  
  
  ! Iterate the model over a few time steps.
  do istep = 1, nstep
  
    ! Calculate the model time.
    time = (istep - 1) * dtime

    ! Create a CARMASTATE for this column.
    call CARMASTATE_Create(cstate, carma_ptr, time, dtime, NZ, &
                          I_CART, I_CART, lat, lon, &
                          xc(:), dx(:), &
                          yc(:), dy(:), &
                          zc(:), zl(:), &
                          p(:), pl(:), &
                          t(:), rc, relhum=relhum(:))
    if (rc /=0) stop "    *** CARMASTATE_Create FAILED ***"

    ! Send the bin mmrs to CARMA
    do ielem = 1, NELEM
      do ibin = 1, NBIN
        call CARMASTATE_SetBin(cstate, ielem, ibin, mmr(:,ielem,ibin), rc)
        if (rc /=0) stop "    *** CARMASTATE_SetBin FAILED ***"
      end do
    end do
      
    ! Execute the step
    call CARMASTATE_Step(cstate, rc, &
      lndfv=lndfv,     ocnfv=ocnfv,     icefv=icefv, &
      lndram=lndram,   ocnram=ocnram,   iceram=iceram, &
      lndfrac=lndfrac, ocnfrac=ocnfrac, icefrac=icefrac)
    if (rc /=0) stop "    *** CARMASTATE_StepFAILED ***"
       
    ! Get the updated bin mmr.
    do ielem = 1, NELEM
      do ibin = 1, NBIN
        call CARMASTATE_GetBin(cstate, ielem, ibin, mmr(:,ielem,ibin), rc, vd=vdry(ibin,ielem))
        if (rc /=0) stop "    *** CARMASTATE_GetBin FAILED ***"
      end do
    end do

    ! Get the updated temperature.
    call CARMASTATE_GetState(cstate, rc, t=t(:))
    if (rc /=0) stop "    *** CARMASTATE_GetState FAILED ***"

    ! Write output for the falltest
    write(lun,'(f12.0)') istep*dtime
    do ielem = 1, NELEM
      do i = 1, NZ
        write(lun,'(2i4,e10.3,e10.3)') ielem, i, real(mmr(i,ielem,OUTBIN)), real(mmr(i,ielem,OUTBIN)*p(i) / 287._f / t(i))
      end do
    end do
  end do   ! time loop

  
  ! Close the output file
  close(unit=lun)
  
  ! write the dry deposition velocity
  open(unit=lun1,file="carma_vdry.txt",status="unknown")
  do igroup = 1, NGROUP
    write(lun1,*) igroup, real(vdry(:, igroup))     
  end do  
  close(unit=lun1)  
  
  ! Cleanup the carma state objects
  call CARMASTATE_Destroy(cstate, rc)
  if (rc /=0) stop "    *** CARMASTATE_Destroy FAILED ***"

  call CARMA_Destroy(carma, rc)
  if (rc /=0) stop "    *** CARMA_Destroy FAILED ***"
end subroutine

