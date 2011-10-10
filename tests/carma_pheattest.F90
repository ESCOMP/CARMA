!! This code is to test the impact of particle heating upon
!! condensational growth.
!!
!! Upon execution, a text file (carma_pheattest.txt) is generated.
!! The text file can be read with the IDL procedure read_pheattest.pro.
!!
!! @author  Chuck Bardeen
!! @version May-2009

program carma_pheattest
  implicit none

  write(*,*) "Particle Heating Test"

  call test_grow_pheat()  
  
  write(*,*) "Done"
end program

!! Just have one grid box. In that grid box, put an initial concentration
!! of drops at the smallest size, then allow that to grow using a gas. The
!! total mas of drops + gas should be conserved.
subroutine test_grow_pheat()
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

  integer, parameter        :: NZ           = 1
  integer, parameter        :: NZP1         = NZ+1
  integer, parameter        :: NELEM        = 1
  integer, parameter        :: NBIN         = 28   ! PMC
!  integer, parameter        :: NBIN         = 18   ! TTL
  integer, parameter        :: NGROUP       = 1
  integer, parameter        :: NSOLUTE      = 0
  integer, parameter        :: NGAS         = 1
  integer, parameter        :: NWAVE        = 4
  integer, parameter        :: LUNOPRT      = 6
  integer, parameter        :: nstep        = 50
  


  ! Different sizes for time steps provide different results
  ! because of the satbility issues.
!  real(kind=f), parameter   :: dtime  = 1._f
!  real(kind=f), parameter   :: dtime  = 5._f
!  real(kind=f), parameter   :: dtime  = 10._f
  real(kind=f), parameter   :: dtime  = 500._f
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
  real(kind=f), allocatable   :: rho(:)
  
  real(kind=f), allocatable   :: mmr(:,:,:)
  real(kind=f), allocatable   :: dtpart(:,:,:)
  real(kind=f), allocatable   :: mmr_gas(:,:)
  real(kind=f), allocatable   :: satliq(:,:)
  real(kind=f), allocatable   :: satice(:,:)
  
  real(kind=f), allocatable   :: r(:)
  real(kind=f), allocatable   :: rmass(:)

  real(kind=f)                :: lat
  real(kind=f)                :: lon
  
  real(kind=f), allocatable          :: wave(:)       ! wavelength centers (cm)
  real(kind=f), allocatable          :: dwave(:)      ! wavelength width (cm)
  real(kind=f), allocatable          :: radint(:,:)   ! radiative intensity (W/m2/sr/cm)
  complex(kind=f), allocatable       :: refidx(:)     ! refractive index
  
  integer               :: i
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
  

!  write(*,*) ""
!  write(*,*) "Particle Growth - Simple"

  ! Open the output text file
  open(unit=lun,file="carma_pheattest.txt",status="unknown")
  
  ! Allocate the arrays that we need for the model
  allocate(xc(NZ), dx(NZ), yc(NZ), dy(NZ), &
           zc(NZ), zl(NZP1), p(NZ), pl(NZP1), &
           t(NZ), rho(NZ)) 
  allocate(mmr(NZ,NELEM,NBIN))
  allocate(dtpart(NZ,NELEM,NBIN))
  allocate(mmr_gas(NZ,NGAS))
  allocate(satliq(NZ,NGAS))
  allocate(satice(NZ,NGAS))
  allocate(r(NBIN))
  allocate(rmass(NBIN))
  allocate(wave(NWAVE), dwave(NWAVE), refidx(NWAVE), radint(NZ,NWAVE))


  ! Define the band centers and widths (in cm) and the refractive indices.
  wave(:)   = (/ 0.26e-4_f, 0.75e-4_f, 3.0e-4_f, 10e-4_f/)
  dwave(:)  = (/ 0.48e-4_f, 0.5e-4_f, 4.0e-4_f, 10e-4_f/)
  refidx(:) = (/ (1.35090_f, 2e-11_f), (1.30590_f, 5.87000e-08_f), (1.03900, 4.38000e-01_f), (1.19260_f, 5.00800e-02_f) /)   


  ! Define the particle-grid extent of the CARMA test
  call CARMA_Create(carma, NBIN, NELEM, NGROUP, NSOLUTE, NGAS, NWAVE, rc, LUNOPRT=LUNOPRT, wave=wave, dwave=dwave)
  if (rc /=0) stop "    *** CARMA_Create FAILED ***"
	carma_ptr => carma


  ! Define the groups

  ! PMC
  rmin  = 2e-8_f
  rmrat = 2.6

  ! TTL
!  rmin  = 1e-4_f
!  rmrat = 2._f
  call CARMAGROUP_Create(carma, 1, "Ice Crystal", rmin, rmrat, I_SPHERE, 1._f, .TRUE., rc, refidx=refidx, do_mie=.true.)
  if (rc /=0) stop "    *** CARMAGROUP_Create FAILED ***"
  
  
  ! Define the elements
  call CARMAELEMENT_Create(carma, 1, 1, "Ice Crystal", RHO_I, I_VOLATILE, I_H2O, rc)
  if (rc /=0) stop "    *** CARMAELEMENT_Create FAILED ***"
  
  ! Define the gases
  call CARMAGAS_Create(carma, 1, "Water Vapor", WTMOL_H2O, I_VAPRTN_H2O_MURPHY2005, I_GCOMP_H2O, rc)
  if (rc /=0) stop "    *** CARMAGAS_Create FAILED ***"

  
  ! Setup the CARMA processes to exercise
  call CARMA_AddGrowth(carma, 1, 1, rc)


  call CARMA_Initialize(carma, rc, do_grow=.true., do_pheat=.true., do_pheatatm=.true., do_thermo=.true.)
!  call CARMA_Initialize(carma, rc, do_grow=.true., do_pheat=.true.)
!  call CARMA_Initialize(carma, rc, do_grow=.true.)
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

  
  ! Write output for the test
  write(lun,*) NGROUP, NELEM, NBIN, NGAS

  do igroup = 1, NGROUP
    call CARMAGROUP_Get(carma, igroup, rc, r=r, rmass=rmass)
    if (rc /=0) stop "    *** CARMAGROUP_Get FAILED ***"
    
    do ibin = 1, NBIN
      write(lun,'(2i4,2e10.3)') igroup, ibin, r(ibin) * 1e4_f, rmass(ibin)
    end do
  end do


  ! Try TTL Conditions ...
  !
  ! p, T, z, mmrgas, rmin, particle concentration
  ! 90 hPa, 190 K, 17 km, 3.5e-6 g/g, 5 um, 0.1/cm^3.
!  p(1)         = 90._f * 100._f
!  zc(1)        = 17000._f
!  t(1)         = 190._f
!  zl(1)        = zc(1) - deltaz
!  zl(2)        = zc(1) + deltaz
!  rho(1)       = (p(1) * 10._f) / (R_AIR * t(1)) * (1e-3_f * 1e6_f)
!  pl(1)        = p(1) - (zl(1) - zc(1)) * rho(1) * (GRAV / 100._f)
!  pl(2)        = p(1) - (zl(2) - zc(1)) * rho(1) * (GRAV / 100._f)
!  mmr_gas(1,1) = 3.5e-6_f
!  mmr(1,1,7)   = (100._f * rmass(7) * (1e-3_f * 1e6_f)) / rho(1)
!  mmr(1,1,8)   = (100._f * rmass(8) * (1e-3_f * 1e6_f)) / rho(1)
  

  ! Try Mesopause Conditions ...
  !
  ! p, T, z, mmrgas, rmin, particle concentration
  ! 0.005 hPa, 150 K, 82.5 km, 5e-6 g/g, 0.1 nm, 100/cm^3.
  p(1)         = 0.005_f * 100._f
  zc(1)        = 82500._f
  t(1)         = 140._f
  zl(1)        = zc(1) - deltaz
  zl(2)        = zc(1) + deltaz
  rho(1)       = (p(1) * 10._f) / (R_AIR * t(1)) * (1e-3_f * 1e6_f)
  pl(1)        = p(1) - (zl(1) - zc(1)) * rho(1) * (GRAV / 100._f)
  pl(2)        = p(1) - (zl(2) - zc(1)) * rho(1) * (GRAV / 100._f)
  mmr_gas(1,1) = 4e-6_f
  mmr(1,1,7)   = (100._f * rmass(7) * (1e-3_f * 1e6_f)) / rho(1)
  mmr(1,1,8)   = (100._f * rmass(8) * (1e-3_f * 1e6_f)) / rho(1)
  
  ! A crude estimate of band radiative intensity per band (in W/m2/sr/cm). Note the band fluxes
  ! need to be scaled by pi to convert to intensity and need to be scaled by
  ! the band width.
!  radint(1,:)  = (/ 171._f, 171._f, 100._f, 300._f/)    ! TTL
  radint(1,:)  = (/ 171._f, 171._f, 60._f, 180._f/)     ! Mesopause
  
  radint(1,:)  = radint(1, :) / 2._f / PI / dwave(:)
  
  
  write(lun,*) 0

  do ielem = 1, NELEM
    do ibin = 1, NBIN
      write(lun,'(2i4,2e10.3)') ielem, ibin, real(mmr(1,ielem,ibin)), real(dtpart(1,ielem,ibin))
    end do
  end do

  do igas = 1, NGAS
    write(lun,'(i4,3e10.3)') igas, real(mmr_gas(1,igas)), 0., 0.
  end do

  write(lun,'(1e12.3)') real(t(1))
		
  ! Iterate the model over a few time steps.
  do istep = 1, nstep
  
    ! Calculate the model time.
    time = (istep - 1) * dtime

    ! Create a CARMASTATE for this column.
    call CARMASTATE_Create(cstate, carma_ptr, time, dtime, NZ, &
                        I_CART, I_CART, lat, lon, &
                        xc(:), dx(:), &
                        yc(:), dy(:), &
                        zc(:), zl(:), p(:), &
                        pl(:), t(:), rc, radint=radint)
  
    ! Send the bin mmrs to CARMA
    do ielem = 1, NELEM
      do ibin = 1, NBIN
       call CARMASTATE_SetBin(cstate, ielem, ibin, &
                              mmr(:,ielem,ibin), rc)
      end do
    end do
    
    ! Send the gas mmrs to CARMA
    do igas = 1, NGAS
       call CARMASTATE_SetGas(cstate, igas, &
                              mmr_gas(:,igas), rc)
    end do

    ! Execute the step
    call CARMASTATE_Step(cstate, rc)
     
    ! Get the updated bin mmr.
    do ielem = 1, NELEM
      do ibin = 1, NBIN
       call CARMASTATE_GetBin(cstate, ielem, ibin, &
                              mmr(:,ielem,ibin), rc, dtpart=dtpart(:,ielem,ibin))
      end do
    end do

    ! Get the updated gas mmr.
    do igas = 1, NGAS
       call CARMASTATE_GetGas(cstate, igas, &
                              mmr_gas(:,igas), rc, &
                              satliq=satliq(:,igas), &
                              satice=satice(:,igas))
    end do

    ! Get the updated temperature.
    call CARMASTATE_GetState(cstate, rc, t=t(:))


    ! Write output for the falltest
    write(lun,'(f12.0)') istep*dtime
    do ielem = 1, NELEM
      do ibin = 1, NBIN
        write(lun,'(2i4,2e12.3)') ielem, ibin, real(mmr(1,ielem,ibin)), real(dtpart(1,ielem,ibin))
      end do
    end do

    do igas = 1, NGAS
      write(lun,'(i4,3e12.3)') igas, real(mmr_gas(1,igas)), satliq(1,igas), satice(1,igas)
    end do
    
    write(lun,'(1e12.3)') real(t(1))
  end do   ! time loop
	
  ! Close the output file
  close(unit=lun)	
	
  if (rc /=0) stop "    *** Stepping FAILED ***"

  ! Cleanup the carma state objects
  call CARMASTATE_Destroy(cstate, rc)
  if (rc /=0) stop "    *** CARMASTATE_Destroy FAILED ***"

  call CARMA_Destroy(carma, rc)
  if (rc /=0) stop "    *** CARMA_Destroy FAILED ***"
end subroutine
