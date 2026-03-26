!! This code is to test the impact of particle swelling from
!! relative humidity on sedimentation. Upon execution, a text
!! file (carma_swelltest.txt) is generated.  The text file can
!! be read with the IDL procedure read_swelltest.pro.
!!

program carma_trop_strat_NH4NO3_2nc_test
  implicit none

  write(*,*) "trop_strat_NH4NO3 Test"

  call test_trop_strat_NH4NO3_simple()
end program carma_trop_strat_NH4NO3_2nc_test

!! Just have one grid box. In that grid box, but an initial concentration
!! of drops at the smallest size, then allow that to grow using a gas. The
!! total mas of drops + gas should be conserved.
subroutine test_trop_strat_NH4NO3_simple()
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
  use nc_types_mod
  use test2nc_mod
  use carmadiags_mod
  implicit none

  ! Model dimension
  integer, parameter        :: NZ           = 1
  integer, parameter        :: NY           = 1 ! lat
  integer, parameter        :: NX           = 1 ! lon
  integer, parameter        :: NZP1         = NZ+1
  integer, parameter        :: NELEM        = 13
  integer, parameter        :: NGROUP       = 2
  integer, parameter        :: NBIN         = 20
  integer, parameter        :: NSOLUTE      = 0
  integer, parameter        :: NGAS         = 4
  integer, parameter        :: NWAVE        = 0
  integer, parameter        :: NGROWCORE    = 11   !! total Number of cores in the group that has growcore
 
  integer, parameter        :: LUNOPRT      = 6

  ! Model Grid
  real(kind=f)          :: lat(NY)
  real(kind=f)          :: lon(NX)
  real(kind=f)          :: time
  integer               :: bins(NBIN)

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
  real(kind=f), parameter   :: deltaz = 10000._f
  real(kind=f), parameter   :: zmin   = 145000._f
  integer, parameter        :: nstep  = 72000 / int(dtime)
  real(kind=f), parameter   :: secsperday = 24._f * 3600._f

  integer                    :: maxsubsteps = 32
  integer                    :: maxretries  = 25

  integer               :: nsubsteps
  integer               :: lastsub = 0
  real(kind=f)          :: nretries
  real(kind=f)          :: lastret = 0._f

  ! Indexes for iterations
  integer               :: i
  integer               :: iy = 1
  integer               :: ix = 1
  integer               :: istep
  integer               :: igas
  integer               :: igroup
  integer               :: ielem
  integer               :: ibin

  ! CARMA object
  type(carma_type), target       :: carma
  type(carma_type), pointer      :: carma_ptr
  type(carmastate_type)          :: cstate
  integer                        :: rc = 0

  ! Input for Group/element
  integer, parameter      :: I_GRP_PRSUL     = 1        !! sulfate aerosol
  integer, parameter      :: I_GRP_MXAER     = 2        !! mixed aerosol

  integer, parameter      :: I_ELEM_PRSUL     = 1       !! sulfate aerosol;  nameing needs to only have 2 charaters  before the element name to work with
                                                                !! partsof the code reading different elements
  integer, parameter      :: I_ELEM_MXAER     = 2       !! aerosol
  integer, parameter      :: I_ELEM_MXOC      = 3       !! organics aerosol
  integer, parameter      :: I_ELEM_MXSOA1    = 4       !! secondary organic aerosol
  integer, parameter      :: I_ELEM_MXSOA2    = 5       !! secondary organic aerosol
  integer, parameter      :: I_ELEM_MXSOA3    = 6       !! secondary organic aerosol
  integer, parameter      :: I_ELEM_MXSOA4    = 7       !! secondary organic aerosol
  integer, parameter      :: I_ELEM_MXSOA5    = 8       !! secondary organic aerosol
  integer, parameter      :: I_ELEM_MXBC      = 9       !! black carbon
  integer, parameter      :: I_ELEM_MXDUST    = 10      !! dust aerosol
  integer, parameter      :: I_ELEM_MXSALT    = 11      !! sea salt aerosol
  integer, parameter      :: I_ELEM_MXHNO3    = 12      !! nitrate
  integer, parameter      :: I_ELEM_MXNH3     = 13      !! ammonium

  integer, parameter      :: I_GAS_H2O        = 1       !! water vapor
  integer, parameter      :: I_GAS_H2SO4      = 2       !! sulphuric acid
  integer, parameter      :: I_GAS_HNO3       = 3       !! nitric acid
  integer, parameter      :: I_GAS_NH3        = 4       !! NH3

  real(kind=f), parameter         :: Kappa_OC = 0.5_f      !! hygroscopicity of OC
  real(kind=f), parameter         :: Kappa_SOA = 0.5_f     !! hygroscopicity of SOA
  real(kind=f), parameter         :: Kappa_BC = 0.1_f
  real(kind=f), parameter         :: Kappa_DUST = 0.2_f
  real(kind=f), parameter         :: Kappa_SALT = 1.0_f
  real(kind=f), parameter         :: Kappa_SULF = 0.5_f

  real(kind=f), parameter         :: RHO_obc  = 1.35_f          !! dry density of smoke aerosol
  real(kind=f), parameter         :: RHO_DUST = 2.65_f          !! dry density of dust particles (g/cm^3) -Lin Su
  real(kind=f), parameter         :: RHO_SALT = 2.65_f          !! dry density of sea salt particles (g/cm)
  real(kind=f), parameter         :: RHO_SULFATE  = 1.923_f     !! dry density of sulfate particles (g/cm3)
  real(kind=f), parameter         :: RHO_NITRATE = 1.51_f       !! dry density of nitrate particle (g/cm3)
  real(kind=f), parameter         :: RHO_AMMONIUM = 0.9_f       !! dry density of ammonium solution (g/cm3)

  real(kind=f), parameter     :: rmin_PRSUL     = 3.43e-8_f  ! minimum radius (cm)
  real(kind=f), parameter     :: vmrat_PRSUL    = 3.67_f     ! volume ratio
  real(kind=f), parameter     :: rmin_MXAER     = 5e-7_f     ! minimum radius (cm)
  real(kind=f), parameter     :: vmrat_MXAER    = 2.5588_f    !2.4610_f        ! volume ratio

  real(kind=f), parameter     :: WTMOL_H2SO4    = 98.078479_f    !! molecular weight of sulphuric acid
  real(kind=f), parameter     :: WTMOL_HNO3     = 63.0_f         !! molecular weight of Nitrate acid
  real(kind=f), parameter     :: WTMOL_NH3      = 17.0_f         !! molecular weight of Ammonia
  
  ! Output for Atmospheric variables
  real(kind=f) :: zc(NZ), zl(NZP1), p(NZ), pl(NZP1), t(NZ), rhoa(NZ), rlheat(NZ), deltaT(NZ)
  real(kind=f) :: relhum(NZ)

  ! Output for Gas variables
  real(kind=f) :: mmr_gas(NZ,NGAS), satliq(NZ,NGAS), satice(NZ,NGAS)
  real(kind=f) :: ei(NZ,NGAS), el(NZ,NGAS), wt(NZ,NGAS)

  ! Output for Group variables
  real(kind=f) :: r(NBIN, NGROUP), rlow(NBIN,NGROUP), rup(NBIN,NGROUP), dr(NBIN,NGROUP)
  real(kind=f) :: arat(NBIN, NGROUP), rrat(NBIN, NGROUP), rprat(NBIN, NGROUP)
  real(kind=f) :: rmass(NBIN,NGROUP), dm(NBIN,NGROUP), vol(NBIN,NGROUP)
  real(kind=f) :: df(NBIN,NGROUP), nmon(NBIN,NGROUP)
  real(kind=f) :: kappa(NBIN)

  real(kind=f) :: nd(NZ, NGROUP), ad(NZ, NGROUP), md(NZ, NGROUP), re(NZ, NGROUP), rew(NZ, NGROUP)
  real(kind=f) :: rm(NZ, NGROUP), jn(NZ, NGROUP), mr(NZ, NGROUP), pa(NZ, NGROUP)
  real(kind=f) :: ar(NZ, NGROUP), vm(NZ, NGROUP), ex(NZ, NGROUP), od(NZ, NGROUP)
  real(kind=f) :: wr_bin(NZ, NGROUP, NBIN), nd_bin(NZ, NGROUP, NBIN), ro_bin(NZ, NGROUP, NBIN)
  real(kind=f) :: mr_bin(NZ, NGROUP, NBIN), vd_bin(NZ, NGROUP, NBIN)

  real(kind=f)          :: wave(NWAVE)
  complex(kind=f)       :: refidx(NWAVE)
  real(kind=f)          :: r_refidx(NWAVE)
  real(kind=f)          :: i_refidx(NWAVE)

  ! Output for Element variables
  real(kind=f) :: mr_elm_bin(NZ, NELEM, NBIN)

  ! Netcdf variable
  type(nc_type)     :: ncflgs
  character(len=80) :: filename_out = "trop_strat_test.nc"
  integer           :: outid
  integer           :: date = 00010101
  real(kind=f)      :: fill_value =-9.255963134931783e+061_f

  ! Working variables
  character(len=80) :: sname,sname_gas(NGAS),sname_elm(NELEM), sname_grp(NGROUP)
  real(kind=f)      :: t_orig(NZ)
  real(kind=f)      :: new_gas(NZ, NGAS) ! If you want to increment mmr_gas
  real(kind=f)      :: mmr(NZ, NELEM, NBIN)
  real(kind=f)      :: qext(NWAVE, NBIN, NGROUP)

! Set the wavelengths
  data wave /  0.532_f /  
    
  data r_refidx / 1.448_f/    
      
  data i_refidx / 0._f /
  
  ! Convert wavelength um -> cm
  wave = wave * 1e-4_f
  
  ! Construct a complex refractive index.
  refidx = cmplx(r_refidx, i_refidx)

  ! Define the particle-grid extent of the CARMA test
  call CARMA_Create(carma, NBIN, NELEM, NGROUP, NSOLUTE, NGAS, NWAVE, rc, NGROWCORE=NGROWCORE, &
      LUNOPRT=LUNOPRT)
  if (rc /=0) stop "    *** CARMA_Create FAILED ***"

  carma_ptr => carma

  ! Create Groups
  call CARMAGROUP_Create(carma, I_GRP_PRSUL, "sulfate", rmin_PRSUL, vmrat_PRSUL, I_SPHERE, &
          1._f, .false., rc, irhswell=I_WTPCT_H2SO4, do_wetdep=.false., do_drydep=.true., &
          solfac=0.3_f, scavcoef=0.1_f, is_sulfate=.true., shortname="PURSUL", do_mie=.false., &
          imiertn=I_MIERTN_TOON1981, iopticstype = I_OPTICS_SULFATE)
  if (rc /=0) stop "    *** CARMAGROUP_Create PRSUL FAILED ***"

    call CARMAGROUP_Create(carma, I_GRP_MXAER, "mixed aerosol", rmin_MXAER, vmrat_MXAER, I_SPHERE, &
          1._f, .false., rc, do_wetdep=.false., do_drydep=.true., solfac=0.2_f, &
          scavcoef=0.1_f, shortname="MIXAER", irhswell=I_PETTERS, do_mie=.true., imiertn=I_MIERTN_TOON1981, &
          iopticstype = I_OPTICS_MIXED_YU_H2O, irhswcomp=I_SWG_URBAN, neutral_volfrc=-1._f)
    if (rc /=0) stop "   *** CARMAGROUP_Create MIXAER FAILED ***"

  ! Define the elements
    call CARMAELEMENT_Create(carma, I_ELEM_PRSUL, I_GRP_PRSUL, "Sulfate", &
          RHO_SULFATE, I_VOLATILE, I_H2SO4, rc, shortname="PRSUL", refidx=refidx)
    if (rc /=0) stop " CARMAELEMENT_Create failed."

    call CARMAELEMENT_Create(carma, I_ELEM_MXAER,  I_GRP_MXAER, "Sulfate in mixed sulfate", &
                             RHO_SULFATE, I_VOLATILE, I_H2SO4, rc,  kappa=Kappa_SULF, shortname="MXAER", refidx=refidx)
    if (rc /=0) stop " CARMAELEMENT_Create failed."

    call CARMAELEMENT_Create(carma, I_ELEM_MXOC,   I_GRP_MXAER, "organic carbon", &
                             RHO_obc, I_COREMASS, I_OC, rc, kappa=Kappa_OC, shortname="MXOC")
    if (rc /=0) stop " CARMAELEMENT_Create failed."

    call CARMAELEMENT_Create(carma, I_ELEM_MXSOA1,   I_GRP_MXAER, "secondary organic aerosol1", &
                             RHO_obc, I_COREMASS, I_SOA1, rc, kappa=Kappa_SOA, shortname="MXSOA1")
    if (rc /=0) stop " CARMAELEMENT_Create failed."

    call CARMAELEMENT_Create(carma, I_ELEM_MXSOA2,   I_GRP_MXAER, "secondary organic aerosol2", &
                             RHO_obc, I_COREMASS, I_SOA2, rc, kappa=Kappa_SOA, shortname="MXSOA2")
    if (rc /=0) stop " CARMAELEMENT_Create failed."

    call CARMAELEMENT_Create(carma, I_ELEM_MXSOA3,   I_GRP_MXAER, "secondary organic aerosol3", &
                             RHO_obc, I_COREMASS, I_SOA3, rc, kappa=Kappa_SOA, shortname="MXSOA3")
    if (rc /=0) stop " CARMAELEMENT_Create failed."

    call CARMAELEMENT_Create(carma, I_ELEM_MXSOA4,   I_GRP_MXAER, "secondary organic aerosol4", &
                             RHO_obc, I_COREMASS, I_SOA4, rc, kappa=Kappa_SOA, shortname="MXSOA4")
    if (rc /=0) stop " CARMAELEMENT_Create failed."

    call CARMAELEMENT_Create(carma, I_ELEM_MXSOA5,   I_GRP_MXAER, "secondary organic aerosol5", &
                             RHO_obc, I_COREMASS, I_SOA5, rc, kappa=Kappa_SOA, shortname="MXSOA5")
    if (rc /=0) stop " CARMAELEMENT_Create failed."

    call CARMAELEMENT_Create(carma, I_ELEM_MXBC,   I_GRP_MXAER, "black carbon", &
                             RHO_obc, I_COREMASS, I_BC, rc, kappa=Kappa_BC, shortname="MXBC", refidx=refidx)
    if (rc /=0) stop " CARMAELEMENT_Create failed."

    call CARMAELEMENT_Create(carma, I_ELEM_MXDUST, I_GRP_MXAER, "dust", &
                             RHO_DUST, I_COREMASS, I_DUST, rc,  kappa=Kappa_DUST, shortname="MXDUST", refidx=refidx)
    if (rc /=0) stop " CARMAELEMENT_Create failed."

    call CARMAELEMENT_Create(carma, I_ELEM_MXSALT, I_GRP_MXAER, "SALT in mixed sulfate", &
                             RHO_SALT, I_COREMASS, I_SALT, rc, kappa=Kappa_SALT, shortname="MXSALT")
    if (rc /=0) stop " CARMAELEMENT_Create failed."

    call CARMAELEMENT_Create(carma, I_ELEM_MXHNO3, I_GRP_MXAER, "Nitrate in mixed sulfate", &
                             RHO_NITRATE, I_GROWCORE, I_HNO3, rc, kappa=Kappa_SULF, shortname="MXHNO3")
    if (rc /=0) stop " CARMAELEMENT_Create failed."

    call CARMAELEMENT_Create(carma, I_ELEM_MXNH3, I_GRP_MXAER, "Ammonium in mixed sulfate", &
                             RHO_AMMONIUM, I_GROWCORE, I_NH3, rc, kappa=Kappa_SULF,shortname="MXNH3")
    if (rc /=0) stop "    *** CARMAELEMENT_Create FAILED ***"

  ! Define the gases
  call CARMAGAS_Create(carma, I_GAS_H2O, "Water Vapor", WTMOL_H2O, I_VAPRTN_H2O_MURPHY2005, I_GCOMP_H2O, &
        rc, shortname = "Q")
  if (rc /=0) stop "    *** CARMAGAS_Create FAILED ***"

  call CARMAGAS_Create(carma, I_GAS_H2SO4, "Sulfuric Acid", WTMOL_H2SO4, I_VAPRTN_H2SO4_AYERS1980, &
        I_GCOMP_H2SO4, rc, shortname = "H2SO4")
  if (rc /=0) stop "    *** CARMAGAS_Create FAILED ***"

  call CARMAGAS_Create(carma, I_GAS_HNO3, "Nitrate Acid", WTMOL_HNO3, I_VAPRTN_HNO3_LUO1995_HANSON1988, & 
        I_GCOMP_HNO3, rc, shortname = "HNO3")!, ds_threshold=0.2_f)
  if (rc /=0) stop "CARMA_DefineModel::CARMAGAS_Create failed. "

  call CARMAGAS_Create(carma, I_GAS_NH3, "Ammonia", WTMOL_NH3, I_VAPRTN_H2SO4_AYERS1980, &
                   I_GCOMP_NH3, rc, shortname = "NH3")
  if (rc /=0) stop "CARMA_DefineModel::CARMAGAS_Create failed."
  
  ! Setup the CARMA processes
  !call CARMA_AddGrowth(carma, I_ELEM_PRSUL, I_GAS_H2SO4, rc)
  !if (rc /=0) stop "    *** CARMA_AddGrowth FAILED ***"

  call CARMA_AddGrowth(carma, I_ELEM_MXAER, I_GAS_H2SO4, rc)
  if (rc /=0) stop "    *** CARMA_AddGrowth FAILED ***"

  call CARMA_AddGrowth(carma, I_ELEM_MXHNO3, I_GAS_HNO3, rc)
  if (rc /=0) stop "    *** CARMA_AddGrowth FAILED ***"

  call CARMA_AddGrowth(carma, I_ELEM_MXNH3, I_GAS_NH3, rc)
  if (rc /=0) stop "    *** CARMA_AddGrowth FAILED ***"

  !call CARMA_AddNucleation(carma, I_ELEM_PRSUL, I_ELEM_PRSUL, I_HOMNUC, 0._f, rc, igas=I_GAS_H2SO4)
  !if (rc /=0) stop "    *** CARMA_AddGrowth FAILED ***"

  call CARMA_AddCoagulation(carma, I_GRP_PRSUL, I_GRP_PRSUL, I_GRP_PRSUL, I_COLLEC_FUCHS, rc)
  if (rc /=0) stop "    *** CARMA_AddGrowth FAILED ***"

  call CARMA_AddCoagulation(carma, I_GRP_PRSUL, I_GRP_MXAER, I_GRP_MXAER, I_COLLEC_DATA, rc)
  if (rc /=0) stop "    *** CARMA_AddGrowth FAILED ***"

  call CARMA_AddCoagulation(carma, I_GRP_MXAER, I_GRP_MXAER, I_GRP_MXAER, I_COLLEC_DATA, rc)
  if (rc /=0) stop "    *** CARMA_AddGrowth FAILED ***"

  call CARMA_Initialize(carma, rc, do_grow=.true., do_thermo=.false.,do_substep=.true., minsubsteps=16, maxsubsteps=maxsubsteps, maxretries=maxretries, dt_threshold=1._f, sulfnucl_method='ZhaoTurco')
  if (rc /=0) stop "    *** CARMA_Initialize FAILED ***"

  ! For simplicity of setup, do a case with Cartesian coordinates,
  ! which are specified in this interface in meters.
  !
  ! NOTE: For Cartesian coordinates, the first level is the bottom
  ! of the model (e.g. z = 0), while for sigma and hybrid coordinates
  ! the first level is the top of the model.
  lat(iy) = -40.0_f
  lon(ix) = -105.0_f

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


  do ibin = 1,NBIN
        bins(ibin) = ibin
  end do


   ! p, T, z, mmrgas
   ! 90 hPa, 190 K, 17 km, H2O mmr 3.5e-6 g/g
   p(1)         = 100._f * 100._f ! hPa
   zc(1)        = 17000._f       ! m
   t(1)         = 200._f         ! K
   t_orig = t(1)
   zl(1)        = zc(1) - deltaz
   zl(2)        = zc(1) + deltaz
   rhoa(1)       = (p(1) * 10._f) / (R_AIR * t(1)) * (1e-3_f * 1e6_f)
   pl(1)        = p(1) - (zl(1) - zc(1)) * rhoa(1) * (GRAV / 100._f)
   pl(2)        = p(1) - (zl(2) - zc(1)) * rhoa(1) * (GRAV / 100._f)

   ! Initial H2O and H2SO4 concentrations
   mmr_gas(:,1)  = 3.e-6_f     ! H2O g/g
   mmr_gas(:,2)  = 0.1e-9_f * (98._f / 29._f)    ! H2SO4
   mmr_gas(:,3)  = 10.e-9_f * (63._f / 29._f) ! hno3 mixing ratio
   mmr_gas(:,4)  = 5.e-10_f ! ammonia mixing ratio

   mmr(:,:,:)   = 0._f
   mmr(:,2,:)   = 1.e-15_f ! Initial element concentration
   mmr(:,1,:)   = 1.e-15_f ! Initial element concentration

   ! ---------------- Create a NETCDF file  --------------------
   ! Get flags for definition/creation of variables in the netcdf file
   ! Flags are all set to FALSE. Set to TRUE the flag of the variable you want to save
   call nctype_create(ncflgs, NGAS, NELEM,NGROUP)

   ! Atm flags
   ncflgs%f_nc_atm%f_nc_p    = .TRUE.
   ncflgs%f_nc_atm%f_nc_t    = .TRUE.
   ncflgs%f_nc_atm%f_nc_dt   = .TRUE.
   ncflgs%f_nc_atm%f_nc_z    = .TRUE.
   ncflgs%f_nc_atm%f_nc_rhoa = .TRUE.
   ncflgs%f_nc_atm%f_nc_rlh  = .TRUE.

   ! GAS flags
   do igas = 1, NGAS
           ncflgs%f_nc_gas(igas)%f_nc_mmr = .TRUE.
           ncflgs%f_nc_gas(igas)%f_nc_sl  = .TRUE.
           ncflgs%f_nc_gas(igas)%f_nc_si  = .TRUE.
           ncflgs%f_nc_gas(igas)%f_nc_ei  = .TRUE.
           ncflgs%f_nc_gas(igas)%f_nc_el  = .TRUE.
           ncflgs%f_nc_gas(igas)%f_nc_wt  = .TRUE.
  end do

  ! ELEM flags
  do ielem = 1, NELEM
           ncflgs%f_nc_elem(ielem)%f_nc_mr_elm_bin = .TRUE.
  end do

  ! GROUP flags
  do igroup = 1, NGROUP
           ncflgs%f_nc_group(igroup)%f_nc_r       = .TRUE.
           ncflgs%f_nc_group(igroup)%f_nc_rmass   = .TRUE.
           ncflgs%f_nc_group(igroup)%f_nc_vol     = .TRUE.
           ncflgs%f_nc_group(igroup)%f_nc_dr      = .TRUE.
           ncflgs%f_nc_group(igroup)%f_nc_dm      = .TRUE.
           ncflgs%f_nc_group(igroup)%f_nc_rup     = .TRUE.
           ncflgs%f_nc_group(igroup)%f_nc_rlw     = .TRUE.
           ncflgs%f_nc_group(igroup)%f_nc_rrat    = .TRUE.
           ncflgs%f_nc_group(igroup)%f_nc_rprat   = .TRUE.
           ncflgs%f_nc_group(igroup)%f_nc_arat    = .TRUE.
           ncflgs%f_nc_group(igroup)%f_nc_df      = .TRUE.
           ncflgs%f_nc_group(igroup)%f_nc_nmon    = .TRUE.

           ncflgs%f_nc_group(igroup)%f_nc_nd    = .TRUE.
           ncflgs%f_nc_group(igroup)%f_nc_ad    = .TRUE.
           ncflgs%f_nc_group(igroup)%f_nc_md    = .TRUE.
           ncflgs%f_nc_group(igroup)%f_nc_re    = .TRUE.
           ncflgs%f_nc_group(igroup)%f_nc_rew   = .TRUE.
           ncflgs%f_nc_group(igroup)%f_nc_rm    = .TRUE.
           ncflgs%f_nc_group(igroup)%f_nc_jn    = .TRUE.
           ncflgs%f_nc_group(igroup)%f_nc_mr    = .TRUE.
           ncflgs%f_nc_group(igroup)%f_nc_pa    = .TRUE.
           ncflgs%f_nc_group(igroup)%f_nc_ar    = .TRUE.
           ncflgs%f_nc_group(igroup)%f_nc_vm    = .TRUE.
           ncflgs%f_nc_group(igroup)%f_nc_exvis = .TRUE.
           ncflgs%f_nc_group(igroup)%f_nc_odvis = .TRUE.

           ncflgs%f_nc_group(igroup)%f_nc_wr_bin = .TRUE.
           ncflgs%f_nc_group(igroup)%f_nc_nd_bin = .TRUE.
           ncflgs%f_nc_group(igroup)%f_nc_ro_bin = .TRUE.
           ncflgs%f_nc_group(igroup)%f_nc_mr_bin = .TRUE.
           ncflgs%f_nc_group(igroup)%f_nc_vd_bin = .TRUE.
  end do



  ! Get variables for nc definition
  do igas = 1, NGAS
          call CARMAGAS_Get(carma, igas, rc, shortname=sname_gas(igas))
  end do

  do ielem = 1, NELEM
          call CARMAELEMENT_Get(carma, ielem, rc, shortname=sname_elm(ielem))
  end do

  do igroup = 1, NGROUP
        call CARMAGROUP_Get(carma, igroup, rc, shortname=sname_grp(igroup), &
                r=r(:,igroup), rlow=rlow(:,igroup), rup=rup(:,igroup), dr=dr(:,igroup), &
                arat=arat(:,igroup), rrat=rrat(:,igroup), rprat=rprat(:,igroup), &
                rmass=rmass(:,igroup), dm=dm(:,igroup), vol=vol(:,igroup), &
                df=df(:,igroup), nmon=nmon(:,igroup))
  end do


  ! Creation of ncfile and variable definition
  call ncdef(ncflgs,filename_out, fill_value, NZ, NY, NX, NBIN, NGAS, NELEM, NGROUP, &
         sname_gas(:), sname_elm(:), sname_grp(:), &
         outid)

  ! Fill nc variables that are not time dependent (out of the loop)
  call ncwrt_olp(ncflgs, outid, NZ, NY, NX, NGAS, NELEM, NGROUP, NBIN, &
          p(:), lat, lon, bins, sname_grp, &
          r=r(:,:), rmass=rmass(:,:), vol=vol(:,:), &
          dr=dr(:,:), dm=dm(:,:), rup=rup(:,:), rlow=rlow(:,:), &
          rrat=rrat(:,:), arat=arat(:,:),rprat=rprat(:,:), &
          df=df(:,:), nmon=nmon(:,:))
  ! ---------------------------------------------------------------------------------

  ! Iterate the model over a few time steps.
  do istep = 2, nstep + 1

    ! Calculate the model time.
    time = (istep - 1) * dtime
    t(:) = t(:)-0.5
    
    ! Create a CARMASTATE for this column.
    call CARMASTATE_Create(cstate, carma_ptr, time, dtime, NZ, &
                        I_CART, lat(iy), lon(ix), &
                        zc(:), zl(:), &
                        p(:),  pl(:), &
                        t(:), rc, &
                        told=t(:), &
                        qh2o=mmr_gas(1,:))
    if (rc /=0) stop "    *** CARMASTATE_Create FAILED ***"

    deltaT(:) = t(:) - t_orig(:)


    ! Send the bin mmrs to CARMA
    do ielem = 1, NELEM
      do ibin = 1, NBIN
        !write(*,*) "mmr(:,ielem,ibin)",ielem, mmr(:,ielem,ibin)
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

    ! if (istep == 1) then
    ! new_gas(:,2) = new_gas(:,2) + 0.10 * new_gas(:,2)  ! H2SO4 - add a source of H2SO4, X% of the initial value
    ! end if:wq


!    mmr_gas(:,2) = 100.e-9_f                       ! H2SO4 - reset to the initial condition (i.e. is constant)


    do igas = 1, NGAS
      call CARMASTATE_SetGas(cstate, igas, new_gas(:,igas), rc, &
              mmr_old=mmr_gas(:,igas),&
              satice_old=satice(:,igas), &
              satliq_old=satliq(:,igas))
      if (rc /=0) stop "    *** CARMASTATE_SetGas FAILED ***"
    end do

    ! Write initial conditions in NetCDF file
    if (istep == 2) then
        call carmadiags(cstate, carma, NZ, NBIN, NELEM, NGROUP,NWAVE,  &
                deltaz, &
                nd=nd(:, :), ad=ad(:, :), md=md(:, :), re=re(:, :), rew=rew(:, :), &
                rm=rm(:, :), jn=jn(:, :), mr=mr(:, :), pa=pa(:, :), &
                ar=ar(:, :), vm=vm(:, :), ex=ex(:,:), od=od(:,:), &
                wr_bin=wr_bin(:,:,:), nd_bin=nd_bin(:,:,:), &
                ro_bin=ro_bin(:,:,:), mr_bin=mr_bin(:,:,:), &
                vd_bin=vd_bin(:,:,:), mr_elm_bin=mr_elm_bin(:,:,:) )

        call ncwrt_ilp(ncflgs, outid, NZ, NY, NX, NBIN, NGAS, NELEM, NGROUP, 1,  &
                zc, iy, ix, bins, 0._f, secsperday, date, &
                sname_gas, sname_elm, sname_grp, &
                p, t, zc, deltaT, rhoa, rlheat, &
                mmr=mmr_gas, si=satice, sl=satliq, ei=ei, el=el, wt=wt, &
                nd=nd, ad=ad, md=md, re=re, rew=rew, rm=rm, jn=jn, mr=mr, pa=pa, &
                ar=ar, vm=vm, ex_vis=ex, od_vis=od, &
                vd_bin=vd_bin, nd_bin=nd_bin, wr_bin=wr_bin, ro_bin=ro_bin, mr_bin=mr_bin, &
                mr_elm_bin=mr_elm_bin)
    end if

    ! Execute the step
    call CARMASTATE_Step(cstate, rc)
    if (rc /=0) stop "    *** CARMASTATE_Step FAILED ***"


    ! Get the retry stats and the updated temperature.
    call CARMASTATE_Get(cstate, rc, nsubstep=nsubsteps, nretry=nretries)
    if (rc /=0) stop "    *** CARMASTATE_Get FAILED ***"

    call CARMASTATE_GetState(cstate, rc, t=t(:), rhoa_wet=rhoa(:), rlheat=rlheat(:))
    if (rc /=0) stop "    *** CARMASTATE_GetState FAILED ***"


    ! Get the updated bin mmr.
    do ielem = 1, NELEM
      do ibin = 1, NBIN
        call CARMASTATE_GetBin(cstate, ielem, ibin, mmr(:,ielem,ibin), rc)
        if (rc /=0) stop "    *** CARMASTATE_GetBin FAILED ***"
      end do
    end do

    !write(*,*) "mmr H2SO4", sum(mmr(:,1,:)+mmr(:,2,:)),"mmr HNO3",sum(mmr(:,12,:))

    ! Get the updated gas mmr.
    do igas = 1, NGAS
      call CARMASTATE_GetGas(cstate, igas, NELEM, &
                             mmr_gas(:,igas), rc, &
                             satliq=satliq(:,igas), satice=satice(:,igas), &
                             eqice=ei(:,igas), eqliq=el(:,igas), wtpct=wt(:,igas))
      if (rc /=0) stop "    *** CARMASTATE_GetGas FAILED ***"
      write(*,*) "igas",igas,"wtpct",wt(:,igas)
    end do
    write(*,*) "t(:)",t(:)

    ! Writing output in NETCDF

    call carmadiags(cstate, carma, NZ, NBIN, NELEM, NGROUP, NWAVE, &
                deltaz, &
                nd=nd(:, :), ad=ad(:, :), md=md(:, :), re=re(:, :), rew=rew(:, :), &
                rm=rm(:, :), jn=jn(:, :), mr=mr(:, :), pa=pa(:, :), &
                ar=ar(:, :), vm=vm(:, :), ex=ex(:,:), od=od(:,:), &
                wr_bin=wr_bin(:,:,:), nd_bin=nd_bin(:,:,:), &
                ro_bin=ro_bin(:,:,:), mr_bin=mr_bin(:,:,:), &
                vd_bin=vd_bin(:,:,:), mr_elm_bin=mr_elm_bin(:,:,:))

    call ncwrt_ilp(ncflgs, outid, NZ, NY, NX, NBIN, NGAS, NELEM, NGROUP, istep,  &
                zc, iy, ix, bins, time, secsperday, date, &
                sname_gas, sname_elm, sname_grp, &
                p, t, zc, deltaT, rhoa, rlheat, &
                mmr=mmr_gas, si=satice, sl=satliq, ei=ei, el=el, wt=wt, &
                nd=nd, ad=ad, md=md, re=re, rew=rew, rm=rm, jn=jn, mr=mr, pa=pa, &
                ar=ar, vm=vm, ex_vis=ex, od_vis=od,&
                vd_bin=vd_bin, nd_bin=nd_bin, wr_bin=wr_bin, ro_bin=ro_bin, mr_bin=mr_bin, &
                mr_elm_bin=mr_elm_bin)

  end do   ! time loop


  ! Close the output NETCDF file
  call ncclose(outid)
  if (rc /=0) write(*, *) "    *** FAILED ***, rc=", rc

  ! Cleanup the carma state objects
  call CARMASTATE_Destroy(cstate, rc)
  if (rc /=0) stop "    *** CARMASTATE_Destroy FAILED ***"

  call CARMA_Destroy(carma, rc)
  if (rc /=0) stop "    *** CARMA_Destroy FAILED ***"
end subroutine
