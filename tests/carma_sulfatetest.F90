!! This code is to test the impact of particle swelling from
!! relative humidity on sedimentation. Upon execution, a text
!! file (carma_swelltest.txt) is generated.  The text file can 
!! be read with the IDL procedure read_swelltest.pro.
!!
!! @author  Chuck Bardeen
!! @version May-2009

program carma_sulfatetest
  implicit none

  write(*,*) "Sulfate Test"
  
  call test_sulfate_simple()  
end program

!! Just have one grid box. In that grid box, but an initial concentration
!! of drops at the smallest size, then allow that to grow using a gas. The
!! total mas of drops + gas should be conserved.
subroutine test_sulfate_simple()
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
  integer, parameter        :: NBIN         = 38
  integer, parameter        :: NGROUP       = 1
  integer, parameter        :: NSOLUTE      = 0
  integer, parameter        :: NGAS         = 2
  integer, parameter        :: NWAVE        = 0
  integer, parameter        :: LUNOPRT      = 6
  


  ! Different sizes for time steps provide different results
  ! because of the satbility issues.
!  real(kind=f), parameter   :: dtime  = 1._f
!  real(kind=f), parameter   :: dtime  = 5._f
!   real(kind=f), parameter   :: dtime  = 10._f
!  real(kind=f), parameter   :: dtime  = 50._f
!  real(kind=f), parameter   :: dtime  = 100._f
!  real(kind=f), parameter   :: dtime  = 1000._f
  real(kind=f), parameter   :: dtime  = 1800._f
!  real(kind=f), parameter   :: dtime  = 5000._f
!  real(kind=f), parameter   :: dtime  = 10000._f
!  real(kind=f), parameter   :: dtime  = 50000._f
  real(kind=f), parameter   :: deltax = 100000._f
  real(kind=f), parameter   :: deltay = 100000._f
  real(kind=f), parameter   :: deltaz = 10000._f
   real(kind=f), parameter   :: zmin   = 145000._f

  integer, parameter        :: nstep  = 180000 / int(dtime)

  
  integer, parameter        :: I_H2SO4  = 1

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
  real(kind=f), allocatable   :: mmr_gas(:,:)
  real(kind=f), allocatable   :: new_gas(:,:)
  real(kind=f), allocatable   :: satliq(:,:)
  real(kind=f), allocatable   :: satice(:,:)
  
  real(kind=f), allocatable   :: r(:)
  real(kind=f), allocatable   :: rmass(:)

  real(kind=f)          :: lat
  real(kind=f)          :: lon
  
  integer               :: i
  integer               :: istep
  integer               :: igas
  integer               :: igroup
  integer               :: ielem
  integer               :: ibin
  integer, parameter    :: lun = 42
  integer               :: nsubsteps
  integer               :: lastsub = 0
  
  real(kind=f)          :: nretries
  real(kind=f)          :: lastret = 0._f

  real(kind=f)          :: time
  real(kind=f)          :: rmin, rmrat
  real(kind=f)          :: RHO_SULFATE   
  real(kind=f)          :: drh
  real(kind=f)          :: t_orig
  
  ! Open the output text file
  open(unit=lun,file="carma_sulfatetest.txt",status="unknown")
  
  ! Allocate the arrays that we need for the model
  allocate(xc(NZ), dx(NZ), yc(NZ), dy(NZ), &
           zc(NZ), zl(NZP1), p(NZ), pl(NZP1), &
           t(NZ),rho(NZ)) 
  allocate(mmr(NZ,NELEM,NBIN))
  allocate(mmr_gas(NZ,NGAS))
  allocate(new_gas(NZ,NGAS))
  allocate(satliq(NZ,NGAS))
  allocate(satice(NZ,NGAS))
  allocate(r(NBIN))
  allocate(rmass(NBIN))

  ! Define the particle-grid extent of the CARMA test
  call CARMA_Create(carma, NBIN, NELEM, NGROUP, NSOLUTE, NGAS, NWAVE, rc, &
      LUNOPRT=LUNOPRT)
  if (rc /=0) stop "    *** CARMA_Create FAILED ***"
  
  carma_ptr => carma


  ! Define the groups
  rmrat = 2._f
! rmin  = 1e-8_f
! rmin  = 1e-4_f
  rmin  = 2.e-8_f
  RHO_SULFATE = 1.923_f  ! dry density of sulfate particles (g/cm3)
  
  call CARMAGROUP_Create(carma, 1, "sulfate", rmin, rmrat, I_SPHERE, 1._f, .false., &
                        rc, irhswell=I_WTPCT_H2SO4, do_drydep=.true., &
                        shortname="SULF", is_sulfate=.true.)
  if (rc /=0) stop "    *** CARMAGROUP_Create FAILED ***"
  
  ! Define the elements
  call CARMAELEMENT_Create(carma, 1, 1, "Sulfate", RHO_SULFATE, I_VOLATILE, I_H2SO4, rc, shortname="SULF")    
  if (rc /=0) stop "    *** CARMAELEMENT_Create FAILED ***"
  
  ! Define the gases
  call CARMAGAS_Create(carma, 1, "Water Vapor", WTMOL_H2O, I_VAPRTN_H2O_MURPHY2005, &
    I_GCOMP_H2O, rc, shortname = "Q", dgc_threshold=0.1_f, ds_threshold=0.1_f)    
  if (rc /=0) stop "    *** CARMAGAS_Create FAILED ***"

  call CARMAGAS_Create(carma, 2, "Sulpheric Acid", 98.078479_f, I_VAPRTN_H2SO4_AYERS1980, &
    I_GCOMP_H2SO4, rc, shortname = "H2SO4", dgc_threshold=0.1_f, ds_threshold=0.1_f)
  if (rc /=0) stop "    *** CARMAGAS_Create FAILED ***"
    

  
  ! Setup the CARMA processes to exercise
  call CARMA_AddGrowth(carma, 1, 2, rc)   ! set H2SO4 to be the condensing gas
  if (rc /=0) stop "    *** CARMA_AddGrowth FAILED ***"

   call CARMA_AddNucleation(carma, 1, 1, I_HOMNUC, 0._f, rc, igas=2)
  if (rc /=0) stop "    *** CARMA_AddNucleation FAILED ***"

   call CARMA_AddCoagulation(carma, 1, 1, 1, I_COLLEC_FUCHS, rc)
  if (rc /=0) stop "    *** CARMA_AddCoagulation FAILED ***"


  call CARMA_Initialize(carma, rc, do_grow=.true., do_coag=.true., do_substep=.true., &
          do_thermo=.true., maxretries=16, maxsubsteps=32, dt_threshold=1._f)
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


  ! Try WACCM model top conditions
!  p(1)         = 5.960299999999999e-06_f * 100._f
!  zc(1)        = 145000._f
!  t(1)         = 872.3763285535849_f
!  zl(1)        = zc(1) - deltaz
!  zl(2)        = zc(1) + deltaz
!  rho(1)       = (p(1) * 10._f) / (R_AIR * t(1)) * (1e-3_f * 1e6_f)
!  pl(1)        = p(1) - (zl(1) - zc(1)) * rho(1) * (GRAV / 100._f)
!  pl(2)        = p(1) - (zl(2) - zc(1)) * rho(1) * (GRAV / 100._f)

  ! Initial H2O and H2SO4 concentrations
!  mmr_gas(:,1)  = 8.588236513537504e-09_f     ! H2O
!  mmr_gas(:,2)  = 2.435825934528716e-11_f     ! H2SO4
  


  ! Try TTL Conditions ...
  !
  ! p, T, z, mmrgas, rmin, particle concentration
  ! 90 hPa, 190 K, 17 km, H2O mmr 3.5e-6 g/g, H2SO4 mmr 100 ppb
  p(1)         = 90._f * 100._f
  zc(1)        = 17000._f
  t(1)         = 250._f
  zl(1)        = zc(1) - deltaz
  zl(2)        = zc(1) + deltaz
  rho(1)       = (p(1) * 10._f) / (R_AIR * t(1)) * (1e-3_f * 1e6_f)
  pl(1)        = p(1) - (zl(1) - zc(1)) * rho(1) * (GRAV / 100._f)
  pl(2)        = p(1) - (zl(2) - zc(1)) * rho(1) * (GRAV / 100._f)

  ! Initial H2O and H2SO4 concentrations
!!  mmr_gas(:,1)  = 3.5e-6_f     ! H2O
  mmr_gas(:,1)  = 100.e-6_f     ! H2O
!!  mmr_gas(:,2)  = 100.e-9_f    ! H2SO4
!!  mmr_gas(:,2)  = 30.e-9_f    ! H2SO4
  mmr_gas(:,2)  = 0.1e-9_f * (98._f / 29._f)    ! H2SO4

  satliq(:,:)   = -1._f
  satice(:,:)   = -1._f
  
  ! Initial sulfate concentration
  mmr(:,:,:)   = 0._f
  
  t_orig = t(1)

  
  write(lun,*) 0

  write(lun,'(2i6,g16.5)') 0, 0, 0._f

  do ielem = 1, NELEM
    do ibin = 1, NBIN
      write(lun,'(2i4,e10.3)') ielem, ibin, real(mmr(1,ielem,ibin))
    end do
  end do

  do igas = 1, NGAS
    write(lun,'(i4,3e10.3)') igas, real(mmr_gas(1,igas)), 0., 0.
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
                        p(:),  pl(:), &
                        t(:), rc, &
                        told=t(:))
    if (rc /=0) stop "    *** CARMASTATE_Create FAILED ***"

    ! Send the bin mmrs to CARMA
    do ielem = 1, NELEM
      do ibin = 1, NBIN
        call CARMASTATE_SetBin(cstate, ielem, ibin, mmr(:,ielem,ibin), rc)
        if (rc /=0) stop "    *** CARMASTATE_SetBin FAILED ***"
      end do
    end do
      
    ! Send the gas mmrs to CARMA
    !
    ! For substepping to do anything, during a step, the old an current
    ! gas mmrs or temperatures need to be different.
    
    ! If you want to add some H2SO4, you can do it here using one or the other
    ! of theses lines.
    new_gas = mmr_gas(:,:)
    
    if (istep == 1) then
      new_gas(:,2) = new_gas(:,2) + .05 * new_gas(:,2)  ! H2SO4 - add a source of H2SO4, 5% of the initial value
    end if

!    mmr_gas(:,2) = 100.e-9_f                       ! H2SO4 - reset to the initial condition (i.e. is constant)


    do igas = 1, NGAS
      call CARMASTATE_SetGas(cstate, igas, new_gas(:,igas), rc, &
              mmr_old=mmr_gas(:,igas),&
              satice_old=satice(:,igas), &
              satliq_old=satliq(:,igas))
      if (rc /=0) stop "    *** CARMASTATE_SetGas FAILED ***"
    end do

    ! Execute the step
    call CARMASTATE_Step(cstate, rc)
    if (rc /=0) stop "    *** CARMASTATE_Step FAILED ***"
     

    ! Get the retry stats and the updated temperature.
    call CARMASTATE_Get(cstate, rc, nsubstep=nsubsteps, nretry=nretries)
    if (rc /=0) stop "    *** CARMASTATE_Get FAILED ***"

    call CARMASTATE_GetState(cstate, rc, t=t(:))
    if (rc /=0) stop "    *** CARMASTATE_GetState FAILED ***"

    ! Get the updated bin mmr.
    do ielem = 1, NELEM
      do ibin = 1, NBIN
        call CARMASTATE_GetBin(cstate, ielem, ibin, mmr(:,ielem,ibin), rc)
        if (rc /=0) stop "    *** CARMASTATE_GetBin FAILED ***"
      end do
    end do

    ! Get the updated gas mmr.
    do igas = 1, NGAS
      call CARMASTATE_GetGas(cstate, igas, &
                             mmr_gas(:,igas), rc, &
                             satliq=satliq(:,igas), &
                             satice=satice(:,igas))
      if (rc /=0) stop "    *** CARMASTATE_GetGas FAILED ***"
    end do
      
   
    ! Write output for the sulfatetest
     write(lun,'(f12.0)') istep*dtime

    write(lun,'(2i6,g16.5)') nsubsteps - lastsub, int(nretries - lastret), t(1) - t_orig
    lastsub = nsubsteps
    lastret = nretries

    do ielem = 1, NELEM
      do ibin = 1, NBIN
        write(lun,'(2i4,e12.3)') ielem, ibin, real(mmr(1,ielem,ibin))
      end do
    end do
  
    do igas = 1, NGAS
      write(lun,'(i4,3e12.3)') igas, real(mmr_gas(1,igas)), satliq(1,igas), satice(1,igas)
    end do

  end do   ! time loop

  ! Close the output file
  close(unit=lun) 
  
  ! Cleanup the carma state objects
  call CARMASTATE_Destroy(cstate, rc)
  if (rc /=0) stop "    *** CARMASTATE_Destroy FAILED ***"

  call CARMA_Destroy(carma, rc)
  if (rc /=0) stop "    *** CARMA_Destroy FAILED ***"
end subroutine
