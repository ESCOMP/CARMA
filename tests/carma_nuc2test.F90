!! This code is to test the aerosol freezing routines for ice.
!!
!! Upon execution, a text file (carma_nuctest.txt) is generated.
!! The text file can be read with the IDL procedure read_nuctest.pro.
!!
!! @author  Chuck Bardeen
!! @version July-2009

program carma_nuctest
  implicit none

  write(*,*) "Nucleation & Growth Test"

  call test_nuc_ttl()  
  
  write(*,*) "Done"
end program

!! Just have one grid box. In that grid box, but an initial concentration
!! of sulfate drops at the smallest size, then allow that to nucleate ice
!! and then the ice can grow using a gas. The total mass of ice + gas should
!! be conserved.
subroutine test_nuc_ttl()
  use carma_precision_mod 
  use carma_constants_mod 
  use carma_enums_mod 
  use carma_types_mod 
  use carmaelement_mod
  use carmagroup_mod
  use carmagas_mod
  use carmasolute_mod
  use carmastate_mod
  use carma_mod
  use atmosphere_mod

  implicit none

  integer, parameter        :: NZ           = 1
  integer, parameter        :: NZP1         = NZ+1
  integer, parameter        :: NELEM        = 3
  integer, parameter        :: NBIN         = 16
  integer, parameter        :: NGROUP       = 2
  integer, parameter        :: NSOLUTE      = 1
  integer, parameter        :: NGAS         = 1
  integer, parameter        :: NWAVE        = 0
  integer, parameter        :: LUNOPRT      = 6
  integer, parameter        :: nstep        = 1000
  


  real(kind=f), parameter   :: dtime  = 1._f
!  real(kind=f), parameter   :: dtime  = 5._f
!  real(kind=f), parameter   :: dtime  = 10._f
!  real(kind=f), parameter   :: dtime  = 20._f
!  real(kind=f), parameter   :: dtime  = 100._f
!  real(kind=f), parameter   :: dtime  = 1000._f
!  real(kind=f), parameter   :: dtime  = 1800._f
!  real(kind=f), parameter   :: dtime  = 5000._f
!  real(kind=f), parameter   :: dtime  = 10000._f
!  real(kind=f), parameter   :: dtime  = 50000._f
  real(kind=f), parameter   :: deltax = 100._f
  real(kind=f), parameter   :: deltay = 100._f
  real(kind=f), parameter   :: deltaz = 100._f
  real(kind=f), parameter   :: zmin   = 3000._f

  real(kind=f), parameter             :: n    = 100._f     !! concentration (cm-3) 
  real(kind=f), parameter             :: r0   = 2.5e-6_f   !! mean radius (cm)
  real(kind=f), parameter             :: rsig = 1.5_f      !! distribution width
  
  real(kind=f)                        :: rhop
  real(kind=f)                        :: rhoa

  integer, parameter                  :: I_H2SO4   = 1               !! sulfate aerosol composition
  integer, parameter                  :: I_ICE     = 2               !! ice
  integer, parameter                  :: I_WATER   = 3               !! water

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
  real(kind=f), allocatable   :: satliq(:,:)
  real(kind=f), allocatable   :: satice(:,:)
  
  real(kind=f), allocatable   :: r(:)
  real(kind=f), allocatable   :: dr(:)
  real(kind=f), allocatable   :: rmass(:)

  real(kind=f)                :: lat
  real(kind=f)                :: lon
  
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
  
  real(kind=f)          :: RHO_CN = 1.78_f
  

  ! Open the output text file
  open(unit=lun,file="carma_nuc2test.txt",status="unknown")
  
  ! Allocate the arrays that we need for the model
  allocate(xc(NZ), dx(NZ), yc(NZ), dy(NZ), &
           zc(NZ), zl(NZP1), p(NZ), pl(NZP1), &
           t(NZ), rho(NZ)) 
  allocate(mmr(NZ,NELEM,NBIN))
  allocate(mmr_gas(NZ,NGAS))
  allocate(satliq(NZ,NGAS))
  allocate(satice(NZ,NGAS))
  allocate(r(NBIN))
  allocate(dr(NBIN))
  allocate(rmass(NBIN))


  ! Define the particle-grid extent of the CARMA test
  call CARMA_Create(carma, NBIN, NELEM, NGROUP, NSOLUTE, NGAS, NWAVE, rc, LUNOPRT=LUNOPRT)
  if (rc /=0) stop "    *** FAILED ***"
	carma_ptr => carma


  ! Define the groups
  call CARMAGROUP_Create(carma, 1, "Sulfate IN", 1.e-7_f, 4._f, I_SPHERE, 1._f, .false., &
                         rc, do_wetdep=.true., do_drydep=.false., solfac=0.3_f, &
                         scavcoef=0.1_f, shortname="CRIN", do_mie=.false.)
  if (rc /=0) stop "    *** FAILED ***"

  call CARMAGROUP_Create(carma, 2, "Ice Crystal", 5.e-5_f, 4.0_f, I_SPHERE, 3._f, .true., &
                         rc, do_wetdep=.true., do_drydep=.false., solfac=0.3_f, &
                         scavcoef=0.1_f, shortname="CRICE", do_mie=.false.)
  if (rc /=0) stop "    *** FAILED ***"
  
  
  ! Define the elements
  call CARMAELEMENT_Create(carma, 1, 1, "Sulfate IN", RHO_CN, I_INVOLATILE, I_H2SO4, rc, shortname="CRIN", isolute=1)
  if (rc /=0) stop "    *** FAILED ***"

  call CARMAELEMENT_Create(carma, 2, 2, "Ice Crystal", RHO_I, I_VOLATILE, I_ICE, rc, shortname="CRICE")
  if (rc /=0) stop "    *** FAILED ***"
  
  call CARMAELEMENT_Create(carma, 3, 2, "Core Mass", RHO_CN, I_COREMASS, I_H2SO4, rc, shortname="CRCORE", isolute=1)
  if (rc /=0) stop "    *** FAILED ***"

  
  ! Define the Solutes
  call CARMASOLUTE_Create(carma, 1, "Sulfuric Acid", 2, 98._f, 1.38_f, rc)
  if (rc /=0) stop "    *** FAILED ***"


  ! Define the gases
  call CARMAGAS_Create(carma, 1, "Water Vapor", WTMOL_H2O, I_VAPRTN_H2O_MURPHY2005, I_GCOMP_H2O, rc, shortname='Q')
  if (rc /=0) stop "    *** FAILED ***"

  
  ! Setup the CARMA processes to exercise growth and nucleation.
  call CARMA_AddGrowth(carma, 2, 1, rc)
  if (rc /=0) stop "    *** FAILED ***"

!  call CARMA_AddNucleation(carma, 1, 3, I_AERFREEZE + I_AF_KOOP_2000, 0._f, rc, igas=1, ievp2elem=1)
!  call CARMA_AddNucleation(carma, 1, 3, I_AERFREEZE + I_AF_TABAZADEH_2000, 0._f, rc, igas=1, ievp2elem=1)
  call CARMA_AddNucleation(carma, 1, 3, I_AERFREEZE + I_AF_MOHLER_2010, 0._f, rc, igas=1, ievp2elem=1)
!  call CARMA_AddNucleation(carma, 1, 3, I_AERFREEZE + I_AF_MURRAY_2010, 0._f, rc, igas=1, ievp2elem=1)

!  call CARMA_AddNucleation(carma, 1, 3, I_AERFREEZE + I_AF_KOOP_2000 + I_AF_MURRAY_2010, 0._f, rc, igas=1, ievp2elem=1)
  if (rc /=0) stop "    *** FAILED ***"

!  write(*,*) "  Initialize ..."
  call CARMA_Initialize(carma, rc, do_grow=.true.)
  if (rc /=0) stop "    *** FAILED ***"
  

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
    call CARMAGROUP_Get(carma, igroup, rc, r=r, dr=dr, rmass=rmass)
    if (rc /=0) stop "    *** FAILED ***"
    
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
!  t(1)         = 190._f
!  t(1)         = 205._f
  t(1)         = 220._f
  zl(1)        = zc(1) - deltaz
  zl(2)        = zc(1) + deltaz
  rho(1)       = (p(1) * 10._f) / (R_AIR * t(1)) * (1e-3_f * 1e6_f)
  pl(1)        = p(1) - (zl(1) - zc(1)) * rho(1) * (GRAV / 100._f)
  pl(2)        = p(1) - (zl(2) - zc(1)) * rho(1) * (GRAV / 100._f)
!  mmr_gas(:,:) = 3.5e-6_f
  mmr_gas(:,:) = 4e-5_f
  
  ! Put in an intial distribution of sulfates.
  call CARMAGROUP_Get(carma, 1, rc, r=r, dr=dr, rmass=rmass)
  if (rc /=0) stop "    *** FAILED ***"
  mmr(:,:,:) = 0._f
  do ibin = 1, NBIN
    rhop = (100._f * dr(ibin) / (sqrt(2._f*PI) * r(ibin) * log(rsig))) * exp(-((log(r(ibin)) - log(r0))**2) / (2._f*(log(rsig))**2)) * rmass(ibin)

    ! We don't know rhoa for the initial condition, but assume something typical of
    ! the conditions at 100 mb and 200K. (mb -> dynes, since R_AIR in cgs)
    rhoa = 100._f * 1000._f / (R_AIR) / (200._f)

    mmr(:,1,ibin)   = rhop / rhoa
  end do
  
  
  write(lun,*) 0

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
                        zc(:), zl(:), p(:), &
                        pl(:), t(:), rc)
    if (rc /=0) stop "    *** FAILED ***"
  
    ! Send the bin mmrs to CARMA
    do ielem = 1, NELEM
      do ibin = 1, NBIN
        call CARMASTATE_SetBin(cstate, ielem, ibin, mmr(:,ielem,ibin), rc)
        if (rc /=0) stop "    *** FAILED ***"
      end do
    end do
    
    ! Send the gas mmrs to CARMA
    do igas = 1, NGAS
      call CARMASTATE_SetGas(cstate, igas, mmr_gas(:,igas), rc)
      if (rc /=0) stop "    *** FAILED ***"
    end do

    ! Execute the step
    call CARMASTATE_Step(cstate, rc)
    if (rc /=0) stop "    *** FAILED ***"
     
    ! Get the updated bin mmr.
    do ielem = 1, NELEM
      do ibin = 1, NBIN
        call CARMASTATE_GetBin(cstate, ielem, ibin, mmr(:,ielem,ibin), rc)
        if (rc /=0) stop "    *** FAILED ***"
      end do
    end do

    ! Get the updated gas mmr.
    do igas = 1, NGAS
      call CARMASTATE_GetGas(cstate, igas, &
                              mmr_gas(:,igas), rc, &
                              satliq=satliq(:,igas), &
                              satice=satice(:,igas))
      if (rc /=0) stop "    *** FAILED ***"
    end do

    ! Get the updated temperature.
    call CARMASTATE_GetState(cstate, rc, t=t(:))
    if (rc /=0) stop "    *** FAILED ***"
    
    ! Cool it down ...
    t(1) = t(1) - .05_f

    if (mod(istep, 10) .eq. 0) then
    
      ! Write output for the falltest
      write(lun,'(f12.0)') istep*dtime
      do ielem = 1, NELEM
        do ibin = 1, NBIN
          write(lun,'(2i4,e12.3)') ielem, ibin, real(mmr(1,ielem,ibin))
        end do
      end do
    
      do igas = 1, NGAS
        write(lun,'(i4,3e12.3)') igas, real(mmr_gas(1,igas)), satliq(1,igas), satice(1,igas)
      end do
      
    end if
  end do   ! time loop
	
  ! Cleanup the carma state objects
  call CARMASTATE_Destroy(cstate, rc)
  if (rc /=0) stop "    *** FAILED ***"

  ! Close the output file
  close(unit=lun)	

  call CARMA_Destroy(carma, rc)
  if (rc /=0) stop "    *** FAILED ***"
end subroutine
