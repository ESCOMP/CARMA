!! This code is to test condensational growth.
!!
!! Upon execution, a text file (carma_growtest.txt) is generated.
!! The text file can be read with the IDL procedure read_growtest.pro.
!!
!! @author  Chuck Bardeen
!! @version May-2009

program carma_kappawetrtest
  implicit none

  write(*,*) "Kappawetr Test"

  call test_kappawetr_simple()

  write(*,*) "Done"
end program

!! Just have one grid box. In that grid box, put an initial concentration
!! of drops at the smallest size, then allow that to grow using a gas. The
!! total mass of drops + gas should be conserved.
subroutine test_kappawetr_simple()
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
  integer, parameter        :: NELEM        = 6
  integer, parameter        :: NBIN         = 20
  integer, parameter        :: NGROUP       = 2
  integer, parameter        :: NSOLUTE      = 0
  integer, parameter        :: NGAS         = 2
  integer, parameter        :: NWAVE        = 0
  integer, parameter        :: LUNOPRT      = 6

  ! Different sizes for time steps provide different results
  ! because of the satbility issues.
!  real(kind=f), parameter   :: dtime  = .1_f
!  real(kind=f), parameter   :: dtime  = 1._f
!  real(kind=f), parameter   :: dtime  = 5._f
!  real(kind=f), parameter   :: dtime  = 10._f
!  real(kind=f), parameter   :: dtime  = 100._f
!  real(kind=f), parameter   :: dtime  = 1000._f
  real(kind=f), parameter   :: dtime  = 1800._f

  real(kind=f), parameter   :: deltax = 100._f
  real(kind=f), parameter   :: deltay = 100._f
  real(kind=f), parameter   :: deltaz = 100._f
  real(kind=f), parameter   :: zmin   = 3000._f

  integer, parameter        :: nstep        = 54000 / dtime

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
  real(kind=f), allocatable   :: rlheat(:)

  real(kind=f), allocatable   :: mmr(:,:,:)
  real(kind=f), allocatable   :: mmr_gas(:,:)
  real(kind=f), allocatable   :: satliq(:,:)
  real(kind=f), allocatable   :: satice(:,:)
  real(kind=f), allocatable   :: kappa(:)
  real(kind=f), allocatable   :: r_wet(:)
  real(kind=f), allocatable   :: rhop_wet(:)

  real(kind=f), allocatable   :: r(:)
  real(kind=f), allocatable   :: rmass(:)
  real(kind=f), allocatable   :: dr(:)

  real(kind=f)          :: lat
  real(kind=f)          :: lon

  integer               :: i
  integer               :: istep
  integer               :: igas
  integer               :: igroup
  integer               :: ielem
  integer               :: ibin
  integer, parameter    :: lun = 42

  real(kind=f)          :: time
  real(kind=f)          :: rmin_PRSUL, rmrat_PRSUL,rmin_MXAER, rmrat_MXAER
  real(kind=f)          :: drh
  real(kind=f)          :: t_orig
  real(kind=f)          :: RHO_obc, RHO_SULFATE, RHO_DUST, RHO_SALT

  integer, parameter      :: I_H2SO4          = 1       !! H2SO4 coposition
  integer, parameter      :: I_OC             = 2       !! OC composition
  integer, parameter      :: I_BC             = 3       !! BC composition
  integer, parameter      :: I_DUST           = 4       !! dust composition
  integer, parameter      :: I_SALT           = 5       !! sea salt composition

  integer, parameter              :: I_GRP_PRSUL     = 1       !! sulfate aerosol
  integer, parameter              :: I_GRP_MXAER     = 2       !! mixed aerosol

  integer, parameter              :: I_ELEM_PRSUL     = 1       !! sulfate aerosol;  nameing needs to only have 2 charaters  before the element name to work with
                                                                !! partsof the code reading different elements
  integer, parameter              :: I_ELEM_MXAER     = 2       !! aerosol
  integer, parameter              :: I_ELEM_MXOC      = 3       !! organics aerosol
  integer, parameter              :: I_ELEM_MXBC      = 4       !! black carbon
  integer, parameter              :: I_ELEM_MXDUST    = 5       !! dust aerosol
  integer, parameter              :: I_ELEM_MXSALT    = 6       !! sea salt aerosol

  integer, parameter              :: I_GAS_H2O        = 1              !! water vapor
  integer, parameter              :: I_GAS_H2SO4      = 2              !! sulphuric acid

  real(kind=f),parameter         :: Kappa_OC = 0.5_f      !! hygroscopicity of OC
  real(kind=f),parameter         :: Kappa_BC = 0.1_f
  real(kind=f),parameter         :: Kappa_DUST = 0.2_f
  real(kind=f),parameter         :: Kappa_SALT = 1.0_f
  real(kind=f),parameter         :: Kappa_SULF = 0.5_f

  real(kind=f), parameter        :: WTMOL_H2SO4 = 98.078479_f    !! molecular weight of sulphuric acid

  character(len=32)    :: carma_seasalt_emis  = 'CMS'
  real(kind=f)         :: ocnfrac = 1.
  real(kind=f)         :: u10in(100),u10out(100)
  real(kind=f)         :: SaltFlux(100,NBIN)

  do i = 1,100
    u10in(i) = 0.+i/80.
  end do

  ! Open the output text file
  open(unit=lun,file="carma_kappawetrtest.txt",status="unknown")

  ! Allocate the arrays that we need for the model
  allocate(xc(NZ), dx(NZ), yc(NZ), dy(NZ), &
           zc(NZ), zl(NZP1), p(NZ), pl(NZP1), &
           t(NZ), rho(NZ), rlheat(NZ))
  allocate(mmr(NZ,NELEM,NBIN))
  allocate(mmr_gas(NZ,NGAS))
  allocate(satliq(NZ,NGAS))
  allocate(satice(NZ,NGAS))
  allocate(r(NBIN))
  allocate(rmass(NBIN))
  allocate(dr(NBIN))
  allocate(kappa(NBIN))
  allocate(r_wet(NBIN))
  allocate(rhop_wet(NBIN))
  allocate(relhum(NZ))

  ! Define the particle-grid extent of the CARMA test
  call CARMA_Create(carma, NBIN, NELEM, NGROUP, NSOLUTE, NGAS, NWAVE, rc, LUNOPRT=LUNOPRT)
  if (rc /=0) stop "    *** CARMA_Create FAILED ***"
  carma_ptr => carma

  relhum(:) =0.3

  ! Define the groups
  rmin_PRSUL     = 3.43e-8_f
  rmrat_PRSUL    = 3.67_f
  rmin_MXAER     = 5e-6_f
  rmrat_MXAER    = 2.2588_f

  RHO_obc  = 1.35_f                  !! dry density of smoke aerosol
  RHO_DUST = 2.65_f                  !! dry density of dust particles (g/cm^3) -Lin Su
  RHO_SALT = 2.65_f                  !! dry density of sea salt particles (g/cm)
  RHO_SULFATE  = 1.923_f     !! dry density of sulfate particles (g/cm3)

  call CARMAGROUP_Create(carma, I_GRP_PRSUL, "sulfate", rmin_PRSUL, rmrat_PRSUL, I_SPHERE, 1._f, .false., &
        rc, irhswell=I_WTPCT_H2SO4, do_wetdep=.true., do_drydep=.true., solfac=0.3_f, &
        scavcoef=0.1_f, is_sulfate=.true., shortname="PRSUL")
  if (rc /=0) stop "    *** CARMAGROUP_Create FAILED ***"

  call CARMAGROUP_Create(carma, I_GRP_MXAER, "mixed aerosol", rmin_MXAER, rmrat_MXAER, I_SPHERE, 1._f, .false., &
        rc, do_wetdep=.true., do_drydep=.true., solfac=0.2_f, &
        scavcoef=0.1_f, shortname="MXAER", irhswell=I_PETTERS,neutral_volfrc=-1._f)
  if (rc /=0) stop "    *** CARMAGROUP_Create FAILED ***"

  ! Define the elements
  call CARMAELEMENT_Create(carma, I_ELEM_PRSUL, I_GRP_PRSUL, "Sulfate", &
                           RHO_SULFATE, I_VOLATILE, I_H2SO4, rc, shortname="PRSUL")
  if (rc /=0) stop "    *** CARMAELEMENT_Create FAILED ***"

  call CARMAELEMENT_Create(carma, I_ELEM_MXAER,  I_GRP_MXAER, "Sulfate in mixed sulfate", &
                           RHO_SULFATE, I_VOLATILE, I_H2SO4, rc,  kappa=Kappa_SULF, shortname="MXAER")
  if (rc /=0) stop "    *** CARMAELEMENT_Create FAILED ***"

  call CARMAELEMENT_Create(carma, I_ELEM_MXOC,   I_GRP_MXAER, "organic carbon", &
                           RHO_obc, I_COREMASS, I_OC, rc, kappa=Kappa_OC, shortname="MXOC")
  if (rc /=0) stop "    *** CARMAELEMENT_Create FAILED ***"

  call CARMAELEMENT_Create(carma, I_ELEM_MXBC,   I_GRP_MXAER, "black carbon", &
                           RHO_obc, I_COREMASS, I_BC, rc, kappa=Kappa_BC, shortname="MXBC")
  if (rc /=0) stop "    *** CARMAELEMENT_Create FAILED ***"

  call CARMAELEMENT_Create(carma, I_ELEM_MXDUST, I_GRP_MXAER, "dust", &
                           RHO_DUST, I_COREMASS, I_DUST, rc,  kappa=Kappa_DUST, shortname="MXDUST")
  if (rc /=0) stop "    *** CARMAELEMENT_Create FAILED ***"

  call CARMAELEMENT_Create(carma, I_ELEM_MXSALT, I_GRP_MXAER, "SALT in mixed sulfate", &
                           RHO_SALT, I_COREMASS, I_SALT, rc, kappa=Kappa_SALT, shortname="MXSALT")
  if (rc /=0) stop "    *** CARMAELEMENT_Create FAILED ***"

  ! Define the gases
  call CARMAGAS_Create(carma, I_GAS_H2O, "Water Vapor", WTMOL_H2O, I_VAPRTN_H2O_MURPHY2005, I_GCOMP_H2O, &
        rc, shortname = "Q", ds_threshold=-0.2_f)
  if (rc /=0) stop "    *** CARMAGAS_Create FAILED ***"

  call CARMAGAS_Create(carma, I_GAS_H2SO4, "Sulfuric Acid", WTMOL_H2SO4, I_VAPRTN_H2SO4_AYERS1980, &
        I_GCOMP_H2SO4, rc, shortname = "H2SO4", ds_threshold=-0.2_f)
  if (rc /=0) stop "    *** CARMAGAS_Create FAILED ***"

  ! Setup the CARMA processes to exercise
  !call CARMA_AddGrowth(carma, I_ELEM_PRSUL, I_GAS_H2SO4, rc)
  !if (rc /=0) stop "    *** CARMA_AddGrowth FAILED ***"

  !call CARMA_AddGrowth(carma, I_ELEM_MXAER, I_GAS_H2SO4, rc)
  !if (rc /=0) stop "    *** CARMA_AddGrowth FAILED ***"

  !call CARMA_AddNucleation(carma, I_ELEM_PRSUL, I_ELEM_PRSUL, I_HOMNUC, 0._f, rc, igas=I_GAS_H2SO4)
  !if (rc /=0) stop "    *** CARMA_AddGrowth FAILED ***"

  !call CARMA_AddCoagulation(carma, I_GRP_PRSUL, I_GRP_PRSUL, I_GRP_PRSUL, I_COLLEC_FUCHS, rc)
  !if (rc /=0) stop "    *** CARMA_AddGrowth FAILED ***"

  !call CARMA_AddCoagulation(carma, I_GRP_PRSUL, I_GRP_MXAER, I_GRP_MXAER, I_COLLEC_DATA, rc)
  !if (rc /=0) stop "    *** CARMA_AddGrowth FAILED ***"

  !call CARMA_AddCoagulation(carma, I_GRP_MXAER, I_GRP_MXAER, I_GRP_MXAER, I_COLLEC_DATA, rc)
  !if (rc /=0) stop "    *** CARMA_AddGrowth FAILED ***"

  call CARMA_Initialize(carma, rc, do_grow=.true., do_thermo=.false.)
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

  ! Start with some initial water drops in the smallest bin, which can then grow
  ! to larger sizes in the presence of the water vapor.


  ! Write output for the test
  !write(lun,*) NGROUP, NELEM, NBIN, NGAS

  do igroup = 1, NGROUP
    call CARMAGROUP_Get(carma, igroup, rc, r=r, rmass=rmass)
    if (rc /=0) stop "    *** CARMAGROUP_Get FAILED ***"

    do ibin = 1, NBIN
      !write(lun,'(2i4,2e10.3)') igroup, ibin, r(ibin) * 1e4_f, rmass(ibin)
    end do
  end do


  ! Try TTL Conditions ...
  !
  ! p, T, z, mmrgas, rmin, particle concentration
  ! 90 hPa, 190 K, 17 km, 3.5e-6 g/g, 5 um, 0.1/cm^3.
  p(1)         = 120._f * 100._f
  zc(1)        = 15000._f
  t(1)         = 200._f
  zl(1)        = zc(1) - deltaz
  zl(2)        = zc(1) + deltaz
  rho(1)       = (p(1) * 10._f) / (R_AIR * t(1)) * (1e-3_f * 1e6_f)
  pl(1)        = p(1) - (zl(1) - zc(1)) * rho(1) * (GRAV / 100._f)
  pl(2)        = p(1) - (zl(2) - zc(1)) * rho(1) * (GRAV / 100._f)
  mmr_gas(:,1) = 1.e-6_f  ! water mixing ratio
  mmr_gas(:,2) = 1.e-9_f ! sulfate mixing ratio

  mmr(:,:,:)   = 1.e-10_f

  t_orig = t(1)

  do ielem = 1, NELEM
    do ibin = 1, NBIN
      !write(lun,'(2i4,e10.3)') ielem, ibin, real(mmr(1,ielem,ibin))
    end do
  end do

  do igas = 1, NGAS
    !write(lun,'(i4,3e10.3)') igas, real(mmr_gas(1,igas)), 0., 0.
  end do


  ! Iterate the model over a few time steps.
  do istep = 1, nstep

    ! Calculate the model time.
    time = (istep - 1) * dtime

    t(1) = t(1)-0.5
    !write(*,*) "t(1)",t(1)

    !bin 1 no variation

    !bin 2 assume all aerosol close to 1e-10, sulfate veries from 0, 1e-8
    mmr(:,I_ELEM_MXAER,2) = 1.e-12_f*1.5**(istep)

    !bin 3 assume all aerosol close to 1e-10, OC veries from 0, 1e-8
    mmr(:,I_ELEM_MXOC,3) = 1.e-12_f*1.5**(istep)

    !bin 4 assume all aerosol close to 1e-10, BC veries from 0, 1e-8
    mmr(:,I_ELEM_MXBC,4) = 1.e-12_f*1.5**(istep)

    !bin 5 assume all aerosol close to 1e-10, salt veries from 0, 1e-8
    mmr(:,I_ELEM_MXSALT,5) = 1.e-12_f*1.5**(istep)

    !bin 6 assume all aerosol close to 1e-10, dust veries from 0, 1e-8
    mmr(:,I_ELEM_MXDUST,6) = 1.e-12_f*1.5**(istep)

      ! Create a CARMASTATE for this column.
      call CARMASTATE_Create(cstate, carma_ptr, time, dtime, NZ, &
                          I_CART, I_CART, lat, lon, &
                          xc(:), dx(:), &
                          yc(:), dy(:), &
                          zc(:), zl(:), &
                          p(:),  pl(:), &
                          t(:), rc, &
                          qh2o = mmr_gas(:,1))
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
      call CARMASTATE_SetGas(cstate, igas, mmr_gas(:,igas), rc)
      if (rc /=0) stop "    *** CARMASTATE_SetGas FAILED ***"
    end do

    ! Execute the step
    call CARMASTATE_Step(cstate, rc)
    if (rc /=0) stop "    *** CARMASTATE_Step FAILED ***"

    ! Get the updated bin mmr.
    do ielem = 1, NELEM
      do ibin = 1, NBIN
        call CARMASTATE_GetBin(cstate, ielem, ibin, mmr(:,ielem,ibin), rc)
        if (rc /=0) stop "    *** CARMASTATE_GetBin FAILED ***"
      end do
    end do

    do ibin = 1,NBIN
      call CARMASTATE_GetBin(cstate, 2, ibin, mmr(:,2,ibin), rc, r_wet=r_wet(ibin),rhop_wet=rhop_wet(ibin),kappa=kappa(ibin))
      if (rc /=0) stop "    *** CARMASTATE_GetBin FAILED ***"
    end do

    ! Get the updated gas mmr.
    do igas = 1, NGAS
      call CARMASTATE_GetGas(cstate, igas, &
                             mmr_gas(:,igas), rc, &
                             satliq=satliq(:,igas), &
                             satice=satice(:,igas))
      if (rc /=0) stop "    *** CARMASTATE_GetGas FAILED ***"
    end do

    ! Get the updated temperature.
    call CARMASTATE_GetState(cstate, rc, t=t(:), rlheat=rlheat(:))
    if (rc /=0) stop "    *** CARMASTATE_Get FAILED ***"

    ! Write output for the falltest
    write(lun,'(f12.0,f12.3)') istep*dtime,t(1)


    do ielem = 2, NELEM
      do ibin = 1, 6
        write(lun,'(2i4,e12.3)') ielem, ibin, real(mmr(1,ielem,ibin))
      end do
    end do

    do ibin = 1, 6
      write(lun,'(i4,4e12.4)') ibin,kappa(ibin),r(ibin),r_wet(ibin),rhop_wet(ibin)
    end do

    !write(*,*) "satliq(:,1)",satliq(:,1)
    write(lun,'(e12.3)') satliq(:,1)

  end do   ! time loop


  ! Close the output file
  close(unit=lun)

  ! Cleanup the carma state objects
  call CARMASTATE_Destroy(cstate, rc)
  if (rc /=0) stop "    *** CARMASTATE_Destroy FAILED ***"

  call CARMA_Destroy(carma, rc)
  if (rc /=0) stop "    *** CARMA_Destroy FAILED ***"
end subroutine
