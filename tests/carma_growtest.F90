!! This code is to test the impact of particle swelling from
!! relative humidity on sedimentation. Upon execution, a text
!! file (carma_swelltest.txt) is generated.  The text file can 
!! be read with the IDL procedure read_swelltest.pro.
!!
!! @author  Chuck Bardeen
!! @version May-2009

program carma_growtest
  implicit none

  write(*,*) "Growth Test"

  call test_grow_simple()  
end program

!! Just have one grid box. In that grid box, but an initial concentration
!! of drops at the smallest size, then allow that to grow using a gas. The
!! total mas of drops + gas should be conserved.
subroutine test_grow_simple()
  use carma_precision_mod 
  use carma_constants_mod 
  use carma_enums_mod 
  use carma_types_mod 
  use carmaelement_mod
  use carmagroup_mod
  use carmagas_mod
  use carmastate_mod
  use carma_mod
  use atmosphere_mod

  implicit none

  integer, parameter        :: NX           = 1
  integer, parameter        :: NY           = 1
  integer, parameter        :: NZ           = 1
  integer, parameter        :: NZP1         = NZ+1
  integer, parameter        :: NELEM        = 1
  integer, parameter        :: NBIN         = 40
  integer, parameter        :: NGROUP       = 1
  integer, parameter        :: NSOLUTE      = 0
  integer, parameter        :: NGAS         = 1
  integer, parameter        :: NWAVE        = 0
  integer, parameter        :: LUNOPRT      = 6
  integer, parameter        :: nstep        = 50
  


  ! Different sizes for time steps provide different results
  ! because of the satbility issues.
!  real(kind=f), parameter   :: dtime  = 1._f
!  real(kind=f), parameter   :: dtime  = 5._f
!  real(kind=f), parameter   :: dtime  = 10._f
  real(kind=f), parameter   :: dtime  = 100._f
!  real(kind=f), parameter   :: dtime  = 1000._f
!  real(kind=f), parameter   :: dtime  = 1800._f
!  real(kind=f), parameter   :: dtime  = 5000._f
!  real(kind=f), parameter   :: dtime  = 10000._f
!  real(kind=f), parameter   :: dtime  = 50000._f
  real(kind=f), parameter   :: deltax = 100._f
  real(kind=f), parameter   :: deltay = 100._f
  real(kind=f), parameter   :: deltaz = 100._f
  real(kind=f), parameter   :: rhmin  = .4_f
  real(kind=f), parameter   :: rhmax  = 1.05_f
  real(kind=f), parameter   :: zmin   = 3000._f

  
  integer, parameter        :: I_H2O  = 1

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
  real(kind=f), allocatable   :: relhum(:,:,:)
  real(kind=f), allocatable   :: rho(:,:,:)
  
  real(kind=f), allocatable   :: mmr(:,:,:,:,:)
  real(kind=f), allocatable   :: mmr_gas(:,:,:,:)
  real(kind=f), allocatable   :: satliq(:,:,:,:)
  real(kind=f), allocatable   :: satice(:,:,:,:)
  
  real(kind=f), allocatable   :: r(:)
  real(kind=f), allocatable   :: rmass(:)

  real(kind=f), allocatable          :: lat(:,:)
  real(kind=f), allocatable          :: lon(:,:)
  
  integer               :: i
  integer               :: ix
  integer               :: iy
  integer               :: ixy
  integer               :: istep
  integer               :: igas
  integer               :: igroup
  integer               :: ielem
  integer               :: ibin
  integer               :: ithread
  integer, parameter    :: lun = 42

  real(kind=f)          :: time
  real(kind=f)          :: rmin, rmrat
  real(kind=f)          :: drh
  

  write(*,*) ""
  write(*,*) "Particle Growth - Simple"

  ! Open the output text file
  open(unit=lun,file="carma_growtest.txt",status="unknown")
  
  ! Allocate the arrays that we need for the model
  allocate(xc(NZ,NY,NX), dx(NZ,NY,NX), yc(NZ,NY,NX), dy(NZ,NY,NX), &
           zc(NZ,NY,NX), zl(NZP1,NY,NX), p(NZ,NY,NX), pl(NZP1,NY,NX), &
           t(NZ,NY,NX), rho(NZ,NY,NX)) 
  allocate(mmr(NZ,NY,NX,NELEM,NBIN))
  allocate(mmr_gas(NZ,NY,NX,NGAS))
  allocate(satliq(NZ,NY,NX,NGAS))
  allocate(satice(NZ,NY,NX,NGAS))
  allocate(r(NBIN))
  allocate(rmass(NBIN))
  allocate(lat(NY,NX), lon(NY,NX))  


  ! Define the particle-grid extent of the CARMA test
  write(*,*) "  CARMA_Create ..."
  call CARMA_Create(carma, NBIN, NELEM, NGROUP, NSOLUTE, NGAS, NWAVE, rc, LUNOPRT=LUNOPRT)
  if (rc /=0) stop "    *** FAILED ***"
	carma_ptr => carma


  ! Define the groups
  write(*,*) "  Add Group(s) ..."
  rmrat = 2._f
!  rmin  = 1e-8_f
!  rmin  = 1e-4_f
  rmin  = 5e-4_f
  call CARMAGROUP_Create(carma, 1, "Ice Crystal", rmin, rmrat, I_SPHERE, 1._f, .TRUE., rc)
  if (rc /=0) stop "    *** FAILED ***"
  
  
  ! Define the elements
  write(*,*) "  Add Element(s) ..."
  call CARMAELEMENT_Create(carma, 1, 1, "Ice Crystal", RHO_I, I_VOLATILE, I_H2O, rc)
  if (rc /=0) stop "    *** FAILED ***"
  
  ! Define the gases
  write(*,*) "  Add Gase(s) ..."
  call CARMAGAS_Create(carma, 1, "Water Vapor", WTMOL_H2O, I_VAPRTN_H2O_MURPHY2005, I_GCOMP_H2O, rc)
!  call CARMAGAS_Create(carma, 1, "Water Vapor", WTMOL_H2O, I_VAPRTN_H2O_GOFF1946, I_GCOMP_H2O, rc)
  if (rc /=0) stop "    *** FAILED ***"

  
  ! Setup the CARMA processes to exercise
  call CARMA_AddGrowth(carma, 1, 1, rc)


  write(*,*) "  Initialize ..."
  call CARMA_Initialize(carma, rc, do_grow=.true.)
  if (rc /=0) stop "    *** FAILED ***"
  

  ! Print the Group Information
  write(*,*)  ""
  call dumpGroup(carma, rc)
  if (rc /=0) stop "    *** FAILED ***"
  
  ! Print the Element Information
  write(*,*)  ""
  call dumpElement(carma, rc)
  if (rc /=0) stop "    *** FAILED ***"

  ! Print the Gas Information
  write(*,*)  ""
  call dumpGas(carma, rc)
  if (rc /=0) stop "    *** FAILED ***"

  write(*,*) ""
  
  
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

  ! Vertical edge
  do i = 1, NZP1
    zl(i,:,:) = zmin + ((i - 1) * deltaz)
  end do
  call GetStandardAtmosphere(zl, p=pl)

  
  ! Setup up an arbitray mass mixing ratio of water vapor, so there is someting to
  ! grow the particles.
  mmr_gas(:,:,:,:) = 1e-2_f
  			
  ! Start with some initial water drops in the smallest bin, which can then grow
  ! to larger sizes in the presence of the water vapor.
  mmr(:,:,:,:,1) = 1e-6_f
  
  
  ! Write output for the test
  write(lun,*) NGROUP, NELEM, NBIN, NGAS

  do igroup = 1, NGROUP
    call CARMAGROUP_Get(carma, igroup, rc, r=r, rmass=rmass)
    if (rc /=0) stop "    *** FAILED ***"
    
    do ibin = 1, NBIN
      write(lun,'(2i4,2e10.3)') igroup, ibin, r(ibin) * 1e4_f, rmass(ibin)
    end do
  end do


  ! Try TTL Conditions ...
  !
  ! p, T, z, mmrgas, rmin, particle concentration
  ! 90 hPa, 190 K, 17 km, 3.5e-6 g/g, 5 um, 0.1/cm^3.
  p(1,:,:)         = 90._f * 100._f
  zc(1,:,:)        = 17000._f
  t(1,:,:)         = 190._f
  zl(1,:,:)        = zc(1,:,:) - deltaz
  zl(2,:,:)        = zc(1,:,:) + deltaz
  rho(1,:,:)       = (p(1,:,:) * 10._f) / (R_AIR * t(1,:,:)) * (1e-3_f * 1e6_f)
  pl(1,:,:)        = p(1,:,:) - (zl(1,:,:) - zc(1,:,:)) * rho(1,:,:) * (GRAV / 100._f)
  pl(2,:,:)        = p(1,:,:) - (zl(2,:,:) - zc(1,:,:)) * rho(1,:,:) * (GRAV / 100._f)
  mmr_gas(:,:,:,:) = 3.5e-6_f
  mmr(:,:,:,1,1)   = (0.1_f * rmass(1) * (1e-3_f * 1e6_f)) / rho(:,:,:)
  
  write(*,'(a6, 3a12)') "level", "zc", "p", "t"
  write(*,'(a6, 3a12)') "", "(m)", "(Pa)", "(K)"
  do i = 1, NZ
    write(*,'(i6,3f12.3)') i, zc(i,NY,NX), p(i,NY,NX), t(i,NY,NX)
  end do
  
  write(*,*) ""
  write(*,'(a6, 2a12)') "level", "zl", "pl"
  write(*,'(a6, 2a12)') "", "(m)", "(Pa)"
  do i = 1, NZP1
    write(*,'(i6,2f12.3)') i, zl(i,NY,NX), pl(i,NY,NX)
  end do
  
   
  write(lun,*) 0

  do ielem = 1, NELEM
    do ibin = 1, NBIN
      write(lun,'(2i4,e10.3)') ielem, ibin, real(mmr(1,NY,NX,ielem,ibin))
    end do
  end do

  do igas = 1, NGAS
    write(lun,'(i4,3e10.3)') igas, real(mmr_gas(1,NY,NX,igas)), 0., 0.
  end do

		
  ! Iterate the model over a few time steps.
  write(*,*) ""
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
			
      ! Send the gas mmrs to CARMA
      do igas = 1, NGAS
         call CARMASTATE_SetGas(cstate, igas, &
                                mmr_gas(:,iy,ix,igas), rc)
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

      ! Get the updated gas mmr.
      do igas = 1, NGAS
         call CARMASTATE_GetGas(cstate, igas, &
                                mmr_gas(:,iy,ix,igas), rc, &
                                satliq=satliq(:,iy,ix,igas), &
                                satice=satice(:,iy,ix,igas))
      end do

      ! Get the updated temperature.
      call CARMASTATE_GetState(cstate, rc, t=t(:,iy,ix))
    enddo

    ! Write output for the falltest
    write(lun,*) istep*dtime
    do ielem = 1, NELEM
      do ibin = 1, NBIN
        write(lun,'(2i4,e12.3)') ielem, ibin, real(mmr(1,NY,NX,ielem,ibin))
      end do
    end do
  
    do igas = 1, NGAS
      write(lun,'(i4,3e12.3)') igas, real(mmr_gas(1,NY,NX,igas)), satliq(1,NY,NX,igas), satice(1,NY,NX,igas)
    end do
  end do   ! time loop
	
  ! Cleanup the carma state objects
  call CARMASTATE_Destroy(cstate, rc)
  if (rc /=0) stop "    *** FAILED ***"

  ! Close the output file
  close(unit=lun)	
	
  write(*,*)  ""

  if (rc /=0) stop "    *** FAILED ***"

  write(*,*)  ""
  write(*,*) "  CARMA_Destroy() ..."
  call CARMA_Destroy(carma, rc)
  if (rc /=0) stop "    *** FAILED ***"
end subroutine
