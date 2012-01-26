!! This code is to test condensational growth with substepping using the
!! clear sky microphysics.
!!
!! Upon execution, a text file (carma_growintest.txt) is generated.
!! The text file can be read with the IDL procedure read_growintest.pro.
!!
!! @author  Chuck Bardeen
!! @version Jan 2012

program carma_growclrtest
  implicit none

  write(*,*) "Clear Sky Growth Test with Substepping"

  call test_growclr_sub()  
  
  write(*,*) "Done"
end program

!! Just have one grid box. In that grid box, put an initial concentration
!! of drops at the smallest size, then allow that to grow using a gas. The
!! total mass of drops + gas should be conserved.
!!
!! Use substepping, which means that the model should start at ice saturation
!! and then get a perturbation in the first timestep to the specified gas
!! concentration.
subroutine test_growclr_sub()
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
  integer, parameter        :: NELEM        = 2
  integer, parameter        :: NBIN         = 18
  integer, parameter        :: NGROUP       = 2
  integer, parameter        :: NSOLUTE      = 0
  integer, parameter        :: NGAS         = 1
  integer, parameter        :: NWAVE        = 0
  integer, parameter        :: LUNOPRT      = 6
  


  ! Different sizes for time steps provide different results
  ! because of the satbility issues.
  integer                    :: maxsubsteps = 32
  integer                    :: maxretries  = 20
  
!  real(kind=f), parameter   :: dtime  = 1._f
!  real(kind=f), parameter   :: dtime  = 5._f
!  real(kind=f), parameter   :: dtime  = 10._f
!  real(kind=f), parameter   :: dtime  = 60._f
!  real(kind=f), parameter   :: dtime  = 100._f
  real(kind=f), parameter   :: dtime  = 1000._f
!  real(kind=f), parameter   :: dtime  = 1800._f
!  real(kind=f), parameter   :: dtime  = 5000._f
!  real(kind=f), parameter   :: dtime  = 10000._f
!  real(kind=f), parameter   :: dtime  = 50000._f

  real(kind=f), parameter   :: deltax = 100._f
  real(kind=f), parameter   :: deltay = 100._f
  real(kind=f), parameter   :: deltaz = 100._f

  real(kind=f), parameter   :: zmin   = 3000._f

  integer, parameter        :: nstep        = 5000 / dtime
!  integer, parameter        :: nstep        = 50

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
  real(kind=f), allocatable   :: cldfrc(:)
  real(kind=f), allocatable   :: rhcrit(:)
  
  
  real(kind=f), allocatable   :: mmr(:,:,:)
  real(kind=f), allocatable   :: mmr_gas(:,:)
  real(kind=f), allocatable   :: mmr_gas_old(:,:)
  real(kind=f), allocatable   :: satliq(:,:)
  real(kind=f), allocatable   :: satice(:,:)
  real(kind=f), allocatable   :: pvapice(:,:)
  
  real(kind=f), allocatable   :: r(:)
  real(kind=f), allocatable   :: rmass(:)

  real(kind=f)                :: lat
  real(kind=f)                :: lon
  
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
  real(kind=f)          :: drh
  real(kind=f)          :: t_orig
  

  ! Open the output text file
  open(unit=lun,file="carma_growclrtest.txt",status="unknown")
  
  ! Allocate the arrays that we need for the model
  allocate(xc(NZ), dx(NZ), yc(NZ), dy(NZ), &
           zc(NZ), zl(NZP1), p(NZ), pl(NZP1), &
           t(NZ), rho(NZ), cldfrc(NZ), rhcrit(NZ)) 
  allocate(mmr(NZ,NELEM,NBIN))
  allocate(mmr_gas(NZ,NGAS))
  allocate(mmr_gas_old(NZ,NGAS))
  allocate(satliq(NZ,NGAS))
  allocate(satice(NZ,NGAS))
  allocate(pvapice(NZ,NGAS))
  allocate(r(NBIN))
  allocate(rmass(NBIN))


  ! Define the particle-grid extent of the CARMA test
  call CARMA_Create(carma, NBIN, NELEM, NGROUP, NSOLUTE, NGAS, NWAVE, rc, &
      LUNOPRT=LUNOPRT)
  if (rc /=0) stop "    *** CARMA_Create FAILED ***"
  carma_ptr => carma


  ! Define the groups
  rmrat = 2._f
!  rmin  = 1e-8_f
!  rmin  = 1e-4_f
  rmin  = 1e-4_f
  call CARMAGROUP_Create(carma, 1, "In-cloud Ice Crystal", rmin, rmrat, I_SPHERE, 1._f, .TRUE., rc, is_cloud=.true.)
  if (rc /=0) stop "    *** CARMAGROUP_Create FAILED ***"

  call CARMAGROUP_Create(carma, 2, "Clear Sky Ice Crystal", rmin, rmrat, I_SPHERE, 1._f, .TRUE., rc, is_cloud=.false.)
  if (rc /=0) stop "    *** CARMAGROUP_Create FAILED ***"
  
  
  ! Define the elements
  call CARMAELEMENT_Create(carma, 1, 1, "In-cloud Ice Crystal", RHO_I, I_VOLATILE, I_H2O, rc)
  if (rc /=0) stop "    *** CARMAELEMENT_Create FAILED ***"
  
  call CARMAELEMENT_Create(carma, 2, 2, "Clear Sky Ice Crystal", RHO_I, I_VOLATILE, I_H2O, rc)
  if (rc /=0) stop "    *** CARMAELEMENT_Create FAILED ***"

  ! Define the gases
  call CARMAGAS_Create(carma, 1, "Water Vapor", WTMOL_H2O, I_VAPRTN_H2O_MURPHY2005, &
      I_GCOMP_H2O, rc, ds_threshold=0.2_f)
!  call CARMAGAS_Create(carma, 1, "Water Vapor", WTMOL_H2O, I_VAPRTN_H2O_GOFF1946, I_GCOMP_H2O, rc)
  if (rc /=0) stop "    *** CARMAGAS_Create FAILED ***"

  
  ! Setup the CARMA processes to exercise
  call CARMA_AddGrowth(carma, 1, 1, rc)
  if (rc /=0) stop "    ***  CARMA_AddGrowth FAILED ***"

  call CARMA_AddGrowth(carma, 2, 1, rc)
  if (rc /=0) stop "    ***  CARMA_AddGrowth FAILED ***"


  call CARMA_Initialize(carma, rc, do_grow=.true., do_substep=.true., do_thermo=.true., &
     minsubsteps=1, maxsubsteps=maxsubsteps, maxretries=maxretries, dt_threshold=2._f, &
     do_incloud=.true., do_clearsky=.true.)
  if (rc /=0) stop "    ***  CARMA_Initialize FAILED ***"
  
  
  
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
  p(1)         = 90._f * 100._f
  zc(1)        = 17000._f
  t(1)         = 190._f
  zl(1)        = zc(1) - deltaz
  zl(2)        = zc(1) + deltaz
  rho(1)       = (p(1) * 10._f) / (R_AIR * t(1)) * (1e-3_f * 1e6_f)
  pl(1)        = p(1) - (zl(1) - zc(1)) * rho(1) * (GRAV / 100._f)
  pl(2)        = p(1) - (zl(2) - zc(1)) * rho(1) * (GRAV / 100._f)
  
  mmr(:,1,1)   = (0.1_f * rmass(1) * (1e-3_f * 1e6_f)) / rho(:) / 2._f
  mmr(:,2,1)   = (0.1_f * rmass(1) * (1e-3_f * 1e6_f)) / rho(:) / 2._f

  ! Calculate saturation using Murphy and Koop [2005].
  pvapice(:,1)      = exp(9.550426_f - (5723.265_f / t(:)) + (3.53068_f * log(t(:))) - (0.00728332_f * t(:)))
  mmr_gas_old(:,1)  = pvapice(:,1) * 18.0_f / p(:) * 28.9_f
  mmr_gas(:,1)      = 3.5e-6_f
  
  t_orig = t(1)
  
   
  write(lun,*) 0

  write(lun,'(2i6,e12.3)') 0, 0, real(t(1) - t_orig)

  do ielem = 1, NELEM
    do ibin = 1, NBIN
      write(lun,'(2i4,e10.3)') ielem, ibin, real(mmr(1,ielem,ibin))
    end do
  end do

  do igas = 1, NGAS
    write(lun,'(i4,3e10.3)') igas, real(mmr_gas(1,igas)), 0._f, 0._f
  end do

    
  ! Iterate the model over a few time steps.
  do istep = 1, nstep
  
    ! Calculate the model time.
    time = (istep - 1) * dtime

    ! NOTE: This means that there should not be any looping over NX or NY done
    ! in any other CARMA routines. They should only loop over NZ.
    ! Create a CARMASTATE for this column.
    call CARMASTATE_Create(cstate, carma_ptr, time, dtime, NZ, &
                          I_CART, I_CART, lat, lon, &
                          xc(:), dx(:), &
                          yc(:), dy(:), &
                          zc(:), zl(:), &
                          p(:),  pl(:), &
                          t(:),  rc, told=t(:))
    if (rc /=0) stop "    *** CARMASTATE_Create FAILED ***"
    
    ! Send the bin mmrs to CARMA
    do ielem = 1, NELEM
      do ibin = 1, NBIN
       call CARMASTATE_SetBin(cstate, ielem, ibin, mmr(:,ielem,ibin), rc)
       if (rc /=0) stop "    *** CARMASTATE_SetBin FAILED ***"
      end do
    end do
      
    ! Send the gas mmrs to CARMA
    do igas = 1, NGAS
      
      ! For the first time step, set mmr_old to the original value and
      ! then set the current mmr to the new value. This will create a
      ! change in gas during the timestep that substepping can use to
      ! control the stability of the solution.
      if (istep == 1) then
        satice(:, igas) = -1._f
        satliq(:, igas) = -1._f
        
        call CARMASTATE_SetGas(cstate, igas, mmr_gas(:,igas), rc, &
          mmr_old=mmr_gas_old(:,igas), satice_old=satice(:,igas), satliq_old=satliq(:,igas))
        if (rc /=0) stop "    *** CARMASTATE_SetGas FAILED ***"
      else
      
        call CARMASTATE_SetGas(cstate, igas, mmr_gas(:,igas), rc, &
          mmr_old=mmr_gas(:,igas), satice_old=satice(:,igas), satliq_old=satliq(:,igas))
        if (rc /=0) stop "    *** CARMASTATE_SetGas FAILED ***"
      end if
    end do

    ! Execute the step
    cldfrc(:) = 0.5_f
    rhcrit(:) = 0.8_f
    
    call CARMASTATE_Step(cstate, rc, cldfrc, rhcrit)
    if (rc /=0) stop "    *** CARMASTATE_Step FAILED ***"
     
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
                              mmr_gas(:,igas), &
                              rc, & 
                              satliq=satliq(:,igas), &
                              satice=satice(:,igas))
      if (rc /=0) stop "    *** CARMASTATE_GetGas FAILED ***"
    end do

    ! Get the updated temperature.
    call CARMASTATE_Get(cstate, rc, nsubstep=nsubsteps, nretry=nretries)
    if (rc /=0) stop "    *** CARMASTATE_Get FAILED ***"
    call CARMASTATE_GetState(cstate, rc, t=t(:))
    if (rc /=0) stop "    *** CARMASTATE_GetState FAILED ***"
    

    ! Write output for the test.
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

  ! Cleanup the carma state objects
  call CARMASTATE_Destroy(cstate, rc)
  if (rc /=0) stop "    *** CARMASTATE_Destroy FAILED ***"

  ! Close the output file
  close(unit=lun) 
  
  call CARMA_Destroy(carma, rc)
  if (rc /=0) stop "    *** CARMA_Destroy FAILED ***"
end subroutine
