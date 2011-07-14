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

  integer, parameter        :: NX           = 1
  integer, parameter        :: NY           = 1
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
  
  real(kind=f), allocatable   :: xc(:,:,:)
  real(kind=f), allocatable   :: dx(:,:,:)
  real(kind=f), allocatable   :: yc(:,:,:)
  real(kind=f), allocatable   :: dy(:,:,:)
  real(kind=f), allocatable   :: zc(:,:,:)
  real(kind=f), allocatable   :: zl(:,:,:)
  real(kind=f), allocatable   :: p(:,:,:)
  real(kind=f), allocatable   :: pl(:,:,:)
  real(kind=f), allocatable   :: t(:,:,:)
  real(kind=f), allocatable   :: relhum(:,:,:)
  
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
  integer               :: igroup
  integer               :: ithread
  integer, parameter    :: lun = 42
  integer, parameter    :: lun1 = 41
  
  real(kind=f)          :: time
  real(kind=f)          :: rmin, rmrat, rho
  real(kind=f)          :: drh
  
  real(kind=f)          :: surfric, ram, landfrac, ocnfrac, icefrac
  real(kind=f),parameter:: surfric_in = 2._f  ! [m/s]
  real(kind=f),parameter:: ram_in = 40._f     ! [s/m]
  real(kind=f)          :: vdry(NBIN, NGROUP)
  
!  write(*,*) ""
!  write(*,*) "Dry deposition of Particles with density of 2.0 g/cm3"

  ! Open the output text file
  open(unit=lun,file="carma_drydeptest.txt",status="unknown") 
  
  ! Allocate the arrays that we need for the model
  allocate(xc(NZ,NY,NX), dx(NZ,NY,NX), yc(NZ,NY,NX), dy(NZ,NY,NX), &
           zc(NZ,NY,NX), zl(NZP1,NY,NX), p(NZ,NY,NX), pl(NZP1,NY,NX), &
           t(NZ,NY,NX))
  allocate(mmr(NZ,NY,NX,NELEM,NBIN))
  allocate(lat(NY,NX), lon(NY,NX))  
  allocate(relhum(NZ,NY,NX))

  ! Define the particle-grid extent of the CARMA test
!  write(*,*) "  CARMA_Create ..."
  call CARMA_Create(carma, NBIN, NELEM, NGROUP, NSOLUTE, NGAS, NWAVE, rc, LUNOPRT=LUNOPRT)
  if (rc /=0) stop "    *** FAILED ***"
  carma_ptr => carma


  ! Define the groups
!  write(*,*) "  Add Group(s) ..."
  rho = 2.65_f
  rmrat = 4.32_f
  rmin  = 1e-6_f
  call CARMAGROUP_Create(carma, 1, "DryDep", rmin, rmrat, I_SPHERE, 1._f, .FALSE., rc, &
       do_mie = .FALSE., do_wetdep=.FALSE., do_drydep=.TRUE., do_vtran=.TRUE., shortname="DD") 
  if (rc /=0) stop "    *** FAILED ***"

  call CARMAGROUP_Create(carma, 2, "NoDryDep", rmin, rmrat, I_SPHERE, 1._f, .FALSE., rc, &
       do_mie = .FALSE., do_wetdep=.FALSE., do_drydep=.FALSE., do_vtran=.TRUE., shortname="NDD")
  if (rc /=0) stop "    *** FAILED ***"
  
  ! Define the elements
!  write(*,*) "  Add Element(s) ..."
  call CARMAELEMENT_Create(carma, 1, 1, "DryDep", rho, I_INVOLATILE, I_PART, rc)
  if (rc /=0) stop "    *** FAILED ***"

  call CARMAELEMENT_Create(carma, 2, 2, "NoDryDep", rho, I_INVOLATILE, I_PART, rc)
  if (rc /=0) stop "    *** FAILED ***"
  
!  write(*,*) "  CARMA_AddCoagulation(carma, 1, 1, 1, I_COLLEC_DATA, rc) ..."
!  call CARMA_AddCoagulation(carma, 1, 1, 1, I_COLLEC_DATA, rc) 
!  if (rc /=0) stop "    *** FAILED ***"
  
  ! Setup the CARMA processes to exercise
!  write(*,*) "  Initialize ..."
  call CARMA_Initialize(carma, rc, do_vtran=.TRUE., do_drydep=.TRUE.)
  if (rc /=0) stop "    *** FAILED ***"

  ! Print the Group Information
!  write(*,*)  ""
!  call dumpGroup(carma, rc)
!  if (rc /=0) stop "    *** FAILED ***"
  
  ! Print the Element Information
!  write(*,*)  ""
!  call dumpElement(carma, rc)
!  if (rc /=0) stop "    *** FAILED ***"

!  write(*,*) ""
  
  
  ! For simplicity of setup, do a case with Cartesian coordinates,
  ! which are specified in this interface in meters.
  !
  ! NOTE: For Cartesian coordinates, the first level is the bottom 
  ! of the model (e.g. z = 0), while for sigma and hybrid coordinates
  ! the first level is the top of the model.
  lat(:,:) = -40.0_f
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
  
  
  ! Set up some arbitrary relative humidities, with the maximum at the bottom and
  ! minimum at the top. Make the RH at the top 0, just to make sure the code can
  ! handle an RH of 0.
  drh    = (rhmax - rhmin) / (NZ-1)
  
  do i = 1, NZ
    relhum(i,:,:) = rhmax - ((i - 1) * drh)
  end do
  
  relhum(NZ,:,:) = 0._f

        
  ! Put a blob in the model for all elements and bins at 8 km.
  mmr(:,:,:,:,:) = 0._f
  do i = 1, NZ
    do ielem = 1, NELEM
      do ibin = 1, NBIN
        mmr(i,:,:,ielem,ibin) = 1e-10_f * exp(-((zc(i,:,:) - 8.e3_f) / 3.e3_f)**2) / &
                                (p(i,:,:) / 287._f / t(i,:,:))
      end do
    end do
  end do

!  write(*,*)  ""
!  write(*, '(a6, 4a12)') "level", "mmr(i,NY,NX,1)", "mmr(i,NY,NX,2)"
!  do i = 1, NZ
!    write(*, '(i6, 4g12.3)') i, mmr(i,NY,NX,1,1), mmr(i,NY,NX,1,2)
!  end do
  
  
  ! Write output for the falltest
  write(lun,*) NZ, NELEM
  do i = 1, NZ
   write(lun,'(i3,2f10.1)') &
    i, zc(i,NY,NX), zl(i+1,NY,NX)-zl(i,NY,NX)
  end do
  
  write(lun,*) 0
  do ielem = 1, NELEM
    do i = 1, NZ
      write(lun,'(2i4,e10.3,e10.3)') &
      ielem, i, real(mmr(i,NY,NX,ielem,OUTBIN)), real(mmr(i,NY,NX,ielem,OUTBIN)*p(i,NY,NX) / 287._f / t(i,NY,NX))
    end do
  end do
  
  
  ! Iterate the model over a few time steps.
!  write(*,*) ""
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
                              pl(:,iy,ix), t(:,iy,ix), rc, relhum=relhum(:,iy,ix))
       !  
       ! Send the bin mmrs to CARMA
       do ielem = 1, NELEM
        do ibin = 1, NBIN
         call CARMASTATE_SetBin(cstate, ielem, ibin, &
                                mmr(:,iy,ix,ielem,ibin), rc)
        end do
       end do
      
       ! Execute the step
       call CARMASTATE_Step(cstate, rc, surfric=surfric_in, ram = ram_in, landfrac=0._f, &
        ocnfrac=1._f, icefrac=0._f)
       
       ! Get the updated bin mmr.
       do ielem = 1, NELEM
        do ibin = 1, NBIN
         call CARMASTATE_GetBin(cstate, ielem, ibin, &
                                mmr(:,iy,ix,ielem,ibin), rc, vd=vdry(ibin,ielem))
        end do
       end do

       ! Get the updated temperature.
       call CARMASTATE_GetState(cstate, rc, t=t(:,iy,ix))
    enddo

    ! Write output for the falltest
    write(lun,*) istep*dtime
    do ielem = 1, NELEM
      do i = 1, NZ
        write(lun,'(2i4,e10.3,e10.3)') &
        ielem, i, real(mmr(i,NY,NX,ielem,OUTBIN)), real(mmr(i,NY,NX,ielem,OUTBIN)*p(i,NY,NX) / 287._f / t(i,NY,NX))
      end do
    end do
  end do   ! time loop
  
  ! Close the output file
  close(unit=lun)
  
  ! write the dry deposition velocity
  open(unit=lun1,file="carma_vdry.txt",status="unknown")
  write(lun1, *) NGROUP
  do igroup = 1, NGROUP
   !write(lun1,'(i4,16e15.5)') &
    write(lun1,*) &
    igroup, real(vdry(:, igroup))
     
  end do  
  close(unit=lun1)  
  
  ! Cleanup the carma state objects
  call CARMASTATE_Destroy(cstate, rc)
  if (rc /=0) stop "    *** FAILED ***"


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
 
!  do igroup = 1, NGROUP
!    write(*,*) &
!    igroup, vdry(:, igroup) 
!  end do  
  
  if (rc /=0) stop "    *** FAILED ***"

!  write(*,*)  ""
!  write(*,*) "  CARMA_Destroy() ..."
  call CARMA_Destroy(carma, rc)
  if (rc /=0) stop "    *** FAILED ***"
end subroutine

