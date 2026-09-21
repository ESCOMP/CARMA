!! This code is to test condensational growth.
!!
!! Upon execution, a text file (carma_growtest.txt) is generated.
!! The text file can be read with the IDL procedure read_growtest.pro.
!!
!! @author  Chuck Bardeen
!! @version May-2009

program carma_alumsulf_2nc_test
  implicit none

  write(*,*) "alumsulftest Test"

  call test_alumsulftest_simple()
  
end program

!! Just have one grid box. In that grid box, put an initial concentration
!! of drops at the smallest size, then allow that to grow using a gas. The
!! total mass of drops + gas should be conserved.
subroutine test_alumsulftest_simple()
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
  use cangle_mod

  implicit none

  integer, parameter        :: NZ           = 1
  integer, parameter        :: NY           = 1 ! lat
  integer, parameter        :: NX           = 1 ! lon
  integer, parameter        :: NZP1         = NZ+1
  integer, parameter        :: NELEM        = 3
  integer, parameter        :: NBIN         = 20
  integer, parameter        :: NGROUP       = 2
  integer, parameter        :: NSOLUTE      = 0
  integer, parameter        :: NGAS         = 2
  integer, parameter        :: NWAVE        = 0
  integer, parameter        :: LUNOPRT      = 6

    ! Model Grid
  real(kind=f)          :: lat(NY)
  real(kind=f)          :: lon(NX)
  real(kind=f)          :: time
  integer               :: bins(NBIN)

  
  ! Different sizes for time steps provide different results
  ! because of the satbility issues.
!  real(kind=f), parameter   :: dtime  = .1_f
!  real(kind=f), parameter   :: dtime  = 1._f
!  real(kind=f), parameter   :: dtime  = 5._f
!  real(kind=f), parameter   :: dtime  = 10._f
!  real(kind=f), parameter   :: dtime  = 100._f
!  real(kind=f), parameter   :: dtime  = 1000._f
  !real(kind=f), parameter   :: zmin   = 21000._f
  !integer, parameter        :: nstep  = 432000/ int(dtime)!180000 / int(dtime)
  !integer, parameter        :: nstep        = 54000 / dtime
  !real(kind=f), parameter   :: secsperday = 24._f * 3600._f
  !integer, parameter        :: nstep        = 5400 / dtime
  ! real(kind=f), parameter   :: dtime  = 7200._f
  !real(kind=f), parameter   :: dtime  = 36000._f
  !integer, parameter        :: nstep        = 10 
  !real(kind=f), parameter   :: deltaz = 500._f


  !real(kind=f), parameter   :: dtime  = 1800._f
  real(kind=f), parameter   :: dtime  = 3000._f
  real(kind=f), parameter   :: deltaz = 1000._f
  real(kind=f), parameter   :: zmin   = 21000._f
  real(kind=f), parameter   :: secsperday = 24._f * 3600._f
  !integer, parameter        :: nstep  = 180000 / int(dtime)
  integer, parameter        :: nstep  = 360000 / int(dtime)



  !real(kind=f), parameter         :: Kappa_DUST = 0.0_f
  real(kind=f), parameter         :: Kappa_zero = 0.0_f
  real(kind=f), parameter         :: Kappa_DUST = 0.2_f
  real(kind=f), parameter         :: Kappa_SULF = 0.5_f


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
  

  ! Indexes for contact angle calcuation
   real(kind=f)          :: vliq
   real(kind=f)          :: vol_alu_mon1 
   real(kind=f)          :: vfull 
   real(kind=f)          :: area_h2so4_out
   real(kind=f)          :: area_solid_out
   real(kind=f)          :: area_h2so4
   real(kind=f)          :: area_solid
   real(kind=f)          :: hliq
   real(kind=f)          :: rliq

  ! CARMA object
  type(carma_type), target            :: carma
  type(carma_type), pointer           :: carma_ptr
  type(carmastate_type)               :: cstate
  integer                             :: rc = 0


  ! Input for Group/element
   real(kind=f)              :: rmin_PRSUL = 3.43e-8_f   ! radius cm
   real(kind=f)              :: rmrat_PRSUL = 3.67_f      ! bin2 radius = bin1 radius * rmrat^(1./3.), so bin2 ~ 5.29e-8 cm
   !real(kind=f)              :: rmrat_ALUM = 2.0_f       !2.4610_f        ! volume ratio
   real(kind=f)              :: rmrat_ALUM = 1.5_f       !2.4610_f        ! volume ratio
   real(kind=f)              :: rmon_ALUM = 24.e-6_f    ! 240 monomer radius
   real(kind=f)              :: rmin_ALUM = 24.e-6_f    ! 240 monomer radius


  real(kind=f)              :: df(NBIN,NGROUP) = 2.6_f
  !real(kind=f)              :: df(NBIN,NGROUP) = 1.6_f
  !real(kind=f)              :: rmon = 3.43e-8_f  !! cm
  real(kind=f)              :: rmon = 24.0e-6_f  !! cm
  real(kind=f)              :: falpha = 1._f
  real(kind=f)              :: alpha = 30._f    !! contact angle (degrees)

  real(kind=f)              :: RHO_SULFATE  = 1.923_f     !! dry density of sulfate particles (g/cm3)
  real(kind=f)              :: RHO_ALUMINUM = 3.95_f
  real(kind=f)              :: rhogroup(NBIN) = 3.95_f  !alumina 

  integer, parameter      :: I_H2SO4          = 1       !! H2SO4 coposition
  integer, parameter      :: I_ALUMINUM       = 2       !! AL2O3 composition

  integer, parameter      :: I_GRP_PRSUL     = 1       !! sulfate aerosol
  integer, parameter      :: I_GRP_ALSUL     = 2       !! mixed aerosol

  integer, parameter      :: I_ELEM_PRSUL     = 1       !! sulfate aerosol;  nameing needs to only have 2 charaters  before the element name to work with
                                                                !! partsof the code reading different elements
  integer, parameter     :: I_ELEM_ALSUL     = 2       !! sulfate coating around alsulf aerosol
  integer, parameter     :: I_ELEM_ALUM      = 3       !! aluminum aerosol

  integer, parameter      :: I_GAS_H2O        = 1              !! water vapor
  integer, parameter      :: I_GAS_H2SO4      = 2              !! sulphuric acid

  real(kind=f), parameter         :: WTMOL_H2SO4 = 98.078479_f    !! molecular weight of sulphuric acid
  
   real(kind=f)          :: rdry

 
   ! Output for Atmospheric variables
  real(kind=f) :: zc(NZ), zl(NZP1), p(NZ), pl(NZP1), t(NZ), rhoa(NZ), rlheat(NZ), deltaT(NZ)

  ! Output for Gas variables
  real(kind=f) :: mmr_gas(NZ,NGAS), satliq(NZ,NGAS), satice(NZ,NGAS)
  real(kind=f) :: ei(NZ,NGAS), el(NZ,NGAS), wt(NZ,NGAS)

  ! Output for Group variables
  real(kind=f) :: r(NBIN, NGROUP), rlow(NBIN,NGROUP), rup(NBIN,NGROUP), dr(NBIN,NGROUP)
  real(kind=f) :: arat(NBIN, NGROUP), rrat(NBIN, NGROUP), rprat(NBIN, NGROUP)
  real(kind=f) :: rmass(NBIN,NGROUP), dm(NBIN,NGROUP), vol(NBIN,NGROUP), vol_sulf(NBIN,NGROUP), vol_alum(NBIN,NGROUP)
  real(kind=f) :: nmon(NBIN,NGROUP)

  real(kind=f) :: nd(NZ, NGROUP), ad(NZ, NGROUP), md(NZ, NGROUP), re(NZ, NGROUP), rew(NZ, NGROUP)
  real(kind=f) :: rm(NZ, NGROUP), jn(NZ, NGROUP), mr(NZ, NGROUP), pa(NZ, NGROUP)
  real(kind=f) :: ar(NZ, NGROUP), vm(NZ, NGROUP), ex(NZ, NGROUP), od(NZ, NGROUP)
  real(kind=f) :: wr_bin(NZ, NGROUP, NBIN), nd_bin(NZ, NGROUP, NBIN), ro_bin(NZ, NGROUP, NBIN)
  real(kind=f) :: mr_bin(NZ, NGROUP, NBIN), vd_bin(NZ, NGROUP, NBIN)
  real(kind=f) :: adh2so4_bin(NZ, NGROUP, NBIN), adsolid_bin(NZ, NGROUP, NBIN)
  real(kind=f) :: corefrac_bin(NZ, NGROUP, NBIN), shellthick_bin(NZ, NGROUP, NBIN)

  ! Output for Element variables
  real(kind=f) :: mr_elm_bin(NZ, NELEM, NBIN)


  ! Netcdf variable
  type(nc_type)     :: ncflgs
  character(len=80) :: filename_out = "alumsulftest_20bin_e10_1.5_sed.nc"
  !character(len=80) :: filename_out = "alumsulftest_sulf.nc"
  integer           :: outid
  integer           :: date = 00010101
  real(kind=f)      :: fill_value =-9.255963134931783e+061_f

  ! Working variables
  character(len=80) :: sname,sname_gas(NGAS),sname_elm(NELEM), sname_grp(NGROUP)
  real(kind=f)      :: t_orig(NZ)
  real(kind=f)      :: new_gas(NZ, NGAS) ! If you want to increment mmr_gas
  real(kind=f)      :: mmr(NZ, NELEM, NBIN)
  real(kind=f)      :: rhop_dry(NBIN, NGROUP)
  real(kind=f)      :: rhop_wet(NBIN, NGROUP)
  real(kind=f)      :: r_wet(NBIN, NGROUP)
  real(kind=f)      :: numberDensity(NBIN, NGROUP)
  real(kind=f)      :: qext(NWAVE, NBIN, NGROUP)

 ! Initial H2O and H2SO4 concentrations
!   mmr_gas(:,1)  = 100.e-6_f     ! H2O g/g
!   mmr_gas(:,2)  = 0.1e-9_f * (98._f / 29._f)    ! H2SO4
   mmr(:,:,:)   = 0._f ! Initial element concentration

   mmr_gas(:,1) = 3.e-6_f  ! background water gas mixing ratio (kg/kg?)
   !mmr_gas(:,2) = 2.e-14_f  ! background sulfuric acid gas mixing ratio
   mmr_gas(:,2) = 2.e-10_f  ! background sulfuric acid gas mixing ratio

  !specify the initial size distribution
   !mmr(1,1,:)   = (/1.e-10_f,1.e-10_f,1.e-10_f,1.e-10_f,1.e-10_f,1.e-10_f,1.e-10_f,1.e-10_f,1.e-10_f,1.e-10_f/)
   mmr(1,1,:)   = (/2.e-13_f,2.e-13_f,2.e-13_f,1.e-13_f,1.e-13_f,1.e-14_f,1.e-15_f,1.e-14_f,1.e-14_f,1.e-14_f,1.e-12_f,1.e-11_f,1.e-11_f,1.e-10_f,1.e-10_f,1.e-10_f,1.e-10_f,1.e-10_f,1.e-10_f,1.e-10_f/)
   !mmr(1,2,:)   = 1.2  * mmr(1,3,:)  ! mixed alum/sulfate element
   !mmr(1,3,:)   = (/1.e-10_f,1.e-10_f,1.e-10_f,1.e-10_f,1.e-10_f,1.e-10_f,1.e-10_f,0._f,0._f,0._f,0._f,0._f,0._f,0._f,0._f,0._f,0._f,0._f,0._f,0._f/) ! alumina 
   !mmr(1,3,:)   = (/1.e-10_f,1.e-10_f,1.e-10_f,1.e-10_f,1.e-10_f,1.e-10_f,1.e-10_f,1.e-10_f,1.e-10_f,1.e-10_f,1.e-10_f,1.e-10_f,0._f,0._f,0._f,0._f,0._f,0._f,0._f,0._f/) ! al element
   !mmr(1,3,:)   = (/4.e-10_f,1.e-13_f,1.e-15_f,1.e-16_f,1.e-17_f,1.e-18_f,1.e-19_f,0._f,0._f,0._f,0._f,0._f,0._f,0._f,0._f,0._f,0._f,0._f,0._f,0._f/)
   !mmr(1,3,:)   = (/1.e-14_f,0._f,0._f,0._f,0._f,0._f,0._f,0._f,0._f,0._f,0._f,0._f,0._f,0._f,0._f,0._f,0._f,0._f,0._f,0._f/)
   !mmr(1,1,1)   = 5.e-10_f
   mmr(1,3,1:2)   = (/3.e-09_f,1.e-11_f/)
   !mmr(1,3,:)   = (/1.e-12_f,1.e-12_f,0._f,0._f,0._f,0._f,0._f,0._f,0._f,0._f,0._f,0._f,0._f,0._f,0._f,0._f,0._f,0._f,0._f,0._f/) ! alumina 
   !mmr(1,2,:)   =  0.0*mmr(1,3,:)  ! sulfate coating element
   !mmr(1,1,:)   = 0.0  * mmr(1,3,:)  ! sulfate coating element


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

  write(*,*) 'zl ',zl

  do ibin = 1,NBIN
        bins(ibin) = ibin
  end do


   ! p, T, z, mmrgas
   ! 90 hPa, 190 K, 17 km, H2O mmr 3.5e-6 g/g
   p(1)         = 90._f * 100._f ! hPa
   zc(1)        = 20000._f       ! m
   t(1)         = 250._f         ! K
   t_orig = t(1)
   zl(1)        = zc(1) - deltaz
   zl(2)        = zc(1) + deltaz
   rhoa(1)       = (p(1) * 10._f) / (R_AIR * t(1)) * (1e-3_f * 1e6_f)
   pl(1)        = p(1) - (zl(1) - zc(1)) * rhoa(1) * (GRAV / 100._f)
   pl(2)        = p(1) - (zl(2) - zc(1)) * rhoa(1) * (GRAV / 100._f)

  
  !check if still needed
 ! real(kind=f), allocatable   :: mmr(:,:,:)
 ! real(kind=f), allocatable   :: mmr_gas(:,:)
 ! real(kind=f), allocatable   :: new_gas(:,:)
 ! real(kind=f), allocatable   :: satliq(:,:)
 ! real(kind=f), allocatable   :: satice(:,:)
 ! real(kind=f), allocatable   :: r_wet(:,:)
 ! real(kind=f), allocatable   :: nmr(:,:)
 ! real(kind=f), allocatable   :: numberDensity(:,:)
 ! real(kind=f), allocatable   :: rhop_wet(:,:)
 ! real(kind=f), allocatable   :: rhop_dry(:,:)
  
 ! real(kind=f), allocatable   :: r(:,:)
 ! real(kind=f), allocatable   :: rmass(:,:)
 ! real(kind=f), allocatable   :: nmon(:)
 ! real(kind=f), allocatable   :: dr(:)
!
!  real(kind=f)          :: lat
!  real(kind=f)          :: lon
!  real(kind=f)          :: mmr_sulf 
!  real(kind=f)          :: r_alu
!  real(kind=f)          :: r_wet_mon
!  real(kind=f)          :: r_alu_mon
!   real(kind=f)          :: rmon_b
!  real(kind=f)          :: ad_alu
!  real(kind=f)          :: ad_alusulf
!  real(kind=f)          :: vol_alu
!  real(kind=f)          :: vol_alu_sulf
!  real(kind=f)          :: vol_alu_mon
!  real(kind=f)          :: vol_alu_mon1
!   real(kind=f)          :: mass_mon 
!  real(kind=f)          :: mass_sulf_mon
!  real(kind=f)          :: vol_sulf_mon
!  real(kind=f)          :: rho_test
!  real(kind=f)          :: vol_dry
!  real(kind=f)          :: vol_tot
!  real(kind=f)          :: nmon1 
!  real(kind=f)          :: area_h2so4 
!  real(kind=f)          :: area_h2so4_out
!  real(kind=f)          :: area_solid
!  real(kind=f)          :: area_solid_out
!  real(kind=f)          :: hliq, rliq 
!  real(kind=f)          :: xt_ca 
  
  !integer, parameter    :: lun = 42

  ! Define the particle-grid extent of the CARMA test
  call CARMA_Create(carma, NBIN, NELEM, NGROUP, NSOLUTE, NGAS, NWAVE, rc, LUNOPRT=LUNOPRT)
  if (rc /=0) stop "    *** CARMA_Create FAILED ***"
  carma_ptr => carma

   call CARMAGROUP_Create(carma, I_GRP_PRSUL, "sulfate", rmin_PRSUL, rmrat_PRSUL, I_SPHERE, 1._f, .false., &
                        rc, irhswell=I_WTPCT_H2SO4, do_drydep=.false., &
                        shortname="PRSUL", is_sulfate=.true.)
  if (rc /=0) stop "    *** CARMAGROUP_Create FAILED ***"

  call CARMAGROUP_Create(carma, I_GRP_ALSUL, "aluminum", rmin_ALUM, rmrat_ALUM, &
                 I_SPHERE, 1._f, .false., rc, irhswell=I_PETTERS, do_drydep=.false.,  &
                 is_sulfate=.true., shortname="ALSUL",  &
                 is_fractal=.TRUE., rmon=rmon_ALUM, df=df, falpha=falpha, &
                 rho=rhogroup)
   if (rc /=0) stop "    *** CARMAGROUP_Create FAILED ***"

  ! Define the elements
   call CARMAELEMENT_Create(carma, I_ELEM_PRSUL, I_GRP_PRSUL, "Sulfate", &
                            RHO_SULFATE, I_VOLATILE, I_H2SO4, rc, shortname="PRSULF")
   if (rc /=0) stop "    *** CARMAELEMENT_Create FAILED ***"

   call CARMAELEMENT_Create(carma, I_ELEM_ALSUL, I_GRP_ALSUL, "Sulfate in mixed ALSUL", RHO_SULFATE,  &
                       I_VOLATILE, I_H2SO4, rc, kappa=Kappa_SULF, shortname="ALSULF")
   if (rc /=0) stop "    *** CARMAELEMENT_Create FAILED ***"

   call CARMAELEMENT_Create(carma, I_ELEM_ALUM, I_GRP_ALSUL, "Aluminum", RHO_ALUMINUM, I_COREMASS, I_ALUMINUM, &
                       rc, kappa=Kappa_zero, shortname="ALUM")
   if (rc /=0) stop "    *** CARMAELEMENT_Create FAILED ***"

  ! Define the gases
  call CARMAGAS_Create(carma, I_GAS_H2O, "Water Vapor", WTMOL_H2O, I_VAPRTN_H2O_MURPHY2005, I_GCOMP_H2O, &
        rc, shortname = "Q")
        !rc, shortname = "Q", ds_threshold=-0.2_f)
  if (rc /=0) stop "    *** CARMAGAS_Create FAILED ***"
 
   call CARMAGAS_Create(carma, I_GAS_H2SO4, "Sulfuric Acid", WTMOL_H2SO4, I_VAPRTN_H2SO4_AYERS1980, &
         I_GCOMP_H2SO4, rc, shortname = "H2SO4", ds_threshold=-0.2_f)
   if (rc /=0) stop "    *** CARMAGAS_Create FAILED ***"
  
  ! Setup the CARMA processes to exercise
  ! Ming can comment out one or several of these processes to do sensitivity tests
     call CARMA_AddGrowth(carma, I_ELEM_PRSUL, I_GAS_H2SO4, rc)
     if (rc /=0) stop "    *** CARMA_AddGrowth FAILED ***"

  !  call CARMA_AddGrowth(carma, I_ELEM_ALSUL, I_GAS_H2SO4, rc)
  !  if (rc /=0) stop "    *** CARMA_AddGrowth FAILED ***"
 
    write(*,*) "nuc ", I_HOMNUC, I_ELEM_PRSUL, I_GAS_H2SO4
    call CARMA_AddNucleation(carma, I_ELEM_PRSUL, I_ELEM_PRSUL, I_HOMNUC, 0._f, rc, igas=I_GAS_H2SO4)

    if (rc /=0) stop "    *** CARMA_AddNuc FAILED ***"

    call CARMA_AddCoagulation(carma, I_GRP_PRSUL, I_GRP_PRSUL, I_GRP_PRSUL, I_COLLEC_FUCHS, rc)
    if (rc /=0) stop "    *** CARMA_Addcoag FAILED ***"
 
   call CARMA_AddCoagulation(carma, I_GRP_ALSUL, I_GRP_ALSUL, I_GRP_ALSUL, I_COLLEC_FUCHS, rc)
   if (rc /=0) stop "    *** CARMA_Addcoag FAILED ***"
! 
    call CARMA_AddCoagulation(carma, I_GRP_PRSUL, I_GRP_ALSUL, I_GRP_ALSUL, I_COLLEC_FUCHS, rc)
    if (rc /=0) stop "    *** CARMA_Addcoag FAILED ***"
 
  ! don't touch
  call CARMA_Initialize(carma, rc, do_grow=.true., do_coag=.true., do_substep=.true., &
          do_thermo=.false., do_vtran=.false., maxretries=16, maxsubsteps=32, sulfnucl_method='ZhaoTurco')
  

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
           ncflgs%f_nc_group(igroup)%f_nc_adh2so4_bin = .TRUE.
           ncflgs%f_nc_group(igroup)%f_nc_adsolid_bin = .TRUE.
           ncflgs%f_nc_group(igroup)%f_nc_corefrac_bin = .TRUE.
           ncflgs%f_nc_group(igroup)%f_nc_shellthick_bin = .TRUE.
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

  adh2so4_bin(:,:,:) = 0._f
  adsolid_bin(:,:,:) = 0._f
  corefrac_bin(:,:,:) = 0._f
  shellthick_bin(:,:,:) = 0._f

  ! Iterate the model over a few time steps.
  do istep = 2, nstep + 1
  
    ! Calculate the model time.
    time = (istep - 1) * dtime

    ! Create a CARMASTATE for this column.
    call CARMASTATE_Create(cstate, carma_ptr, time, dtime, NZ, &
                          I_CART, lat(iy), lon(ix), &
                          zc(:), zl(:), &
                          p(:),  pl(:), &
                          t(:), rc, &
                          told=t(:), &
                          qh2o = mmr_gas(:,1))
    if (rc /=0) stop "    *** CARMASTATE_Create FAILED ***"

    deltaT(:) = t(:) - t_orig(:)
    
    ! Send the bin mmrs to CARMA
    do ielem = 1, NELEM
        do ibin = 1, NBIN
         call CARMASTATE_SetBin(cstate, ielem, ibin, mmr(:,ielem,ibin), rc)
         if (rc /=0) stop "    *** CARMASTATE_SetBin FAILED ***"
        end do
    end do
      
    new_gas(:,:) = mmr_gas(:,:)
    mmr_gas(:,2) = new_gas(:,2) + 0.10 * new_gas(:,2)  ! H2SO4 - add a source of H2SO4, 5% of the initial value

    ! Send the gas mmrs to CARMA
    do igas = 1, NGAS
      call CARMASTATE_SetGas(cstate, igas, mmr_gas(:,igas), rc, &   ! mmr_gas need to set to new_gas if the gas change is included.
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
                mr_elm_bin=mr_elm_bin, &
                adh2so4_bin=adh2so4_bin, adsolid_bin=adsolid_bin, &
                corefrac_bin=corefrac_bin, shellthick_bin=shellthick_bin)
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

!    do ibin = 1,NBIN
!     call CARMASTATE_GetBin(cstate, 1, ibin, mmr(:,1,ibin), rc, r_wet=r_wet(ibin,1),rhop_wet=rhop_wet(ibin,1))
!     call CARMASTATE_GetBin(cstate, 2, ibin, mmr(:,2,ibin), rc, numberDensity=numberDensity(ibin,2), nmr=nmr(ibin,2), \
!          r_wet=r_wet(ibin,2),rhop_wet=rhop_wet(ibin,2),rhop_dry=rhop_dry(ibin,2))

!     if (rc /=0) stop "    *** CARMASTATE_GetBin FAILED ***"
!   end do

    ! Get the updated gas mmr.
    do igas = 1, NGAS
      call CARMASTATE_GetGas(cstate, igas, &
                             mmr_gas(:,igas), rc, &
                             satliq=satliq(:,igas), &
                             satice=satice(:,igas))
      if (rc /=0) stop "    *** CARMASTATE_GetGas FAILED ***"
    end do

    ! Get the updated temperature.
    call CARMASTATE_GetState(cstate, rc, t=t(:))
    if (rc /=0) stop "    *** CARMASTATE_Get FAILED ***"


    ! Calculate surface area density 
    ! required input for the routines:
    !alpha   !degrees
    !Vliq    !cm3   => the liquid H2SO4-H2O volume on the particle/momoner
    !Rp      !cm      => the radius of the solid particle monomer on which the H2SO4-H2O is sitting on (r)

    !mmr_sulf = mmr_alusulf - mmr_alu  [kg/kg]
    ! vol = mmr [g/g] / rho_sulf [g/cm3] = [cm3/g]
    
    vliq = 0._f
    vfull = 0._f
    vol_sulf(:,:) = 0._f
    vol_alum(:,:) = 0._f
    vol_alu_mon1= 0._f

    rdry= 0._f
    adh2so4_bin(1,I_GRP_ALSUL,:) = 0._f
    adsolid_bin(1,I_GRP_ALSUL,:) = 0._f
    !rmon_b= 0._f
    !r_wet_mon = 0._f
    !r_alu_mon = 0._f
    !ad_alu = 0._f
    !ad_alusulf = 0._f
    !vol_dry = 0._f
    !vol_tot = 0._f
    !area_h2so4 = 0._f
    !area_h2so4_out = 0._f
    !area_solid = 0._f
    !area_solid_out = 0._f
    !xt_ca = 0._f
    !mass_mon = 0._f
    !mass_sulf_mon= 0._f
    !vol_sulf_mon= 0._f

    ! Calculate volume of the sulfate/H2O per monomer
    ! Volume of the total aerosol per air
    ! Alsulf Group: vol total = 4/3. pi * r_wet^3 * nmr_total -> cm3/cm3 total air 
    ! mmr(1,3,ibin) ! volume of aluminum: 
    ! Alum. vol [cm3/cm3] : mmr_al [g/g] * air density [g/cm3] / rho_alum [g/cm3]
    ! Notes:  mass mixing ratio of any core: core element times #of particles/cm3 / air density

    do ibin = 1, NBIN
       !mass of sulfate coating
       !vol = mass/rho
       call CARMASTATE_GetBin(cstate, 2, ibin, mmr(:,2,ibin), rc, rhop_dry=rhop_dry(ibin,2),r_wet=r_wet(ibin,2),numberDensity=numberDensity(ibin,2))
       call CARMASTATE_GetBin(cstate, 3, ibin, mmr(:,3,ibin), rc)
       !vol_sulf(ibin,2) = mmr(1,2,ibin) * rhop_dry / RHO_SULFATE
       !vol_alum(ibin,2) = mmr(1,3,ibin) * rhop_dry / RHO_ALUMINUM
       rdry = ((3._f * rmass(ibin,2)) / (4._f * PI * rhop_dry(ibin,2)))**(1._f/3._f)

       ! Estimate the actual simulated alumina core volume fraction and H2SO4/H2O
       ! shell thickness for this bin, from the live sulfate:alumina mass ratio
       ! (same ratio already used for vliq below) and the wet/dry radii.
       if (mmr(1,3,ibin) > 0._f) then
         corefrac_bin(1,I_GRP_ALSUL,ibin) = (mmr(1,3,ibin) / RHO_ALUMINUM) / &
              ((mmr(1,3,ibin) / RHO_ALUMINUM) + (mmr(1,2,ibin) / RHO_SULFATE))
         shellthick_bin(1,I_GRP_ALSUL,ibin) = r_wet(ibin,2) - rdry * corefrac_bin(1,I_GRP_ALSUL,ibin)**(1._f/3._f)
       else
         corefrac_bin(1,I_GRP_ALSUL,ibin) = 0._f
         shellthick_bin(1,I_GRP_ALSUL,ibin) = 0._f
       end if


       vol_alu_mon1 = 4._f/3._f *  PI * rmon **3._f   ! monomer mass [kg] per particle
       !vliq/vol_alu_mon1 = (mmr(1,2,ibin) / mmr(1,3,ibin) *  (RHO_ALUMINUM / RHO_SULFATE)
       ! Guard against 0./0. (NaN) when a bin has not yet received any aluminum core mass.
       if (mmr(1,3,ibin) > 0._f) then
         vliq = vol_alu_mon1 * (mmr(1,2,ibin) / mmr(1,3,ibin)) *  (RHO_ALUMINUM / RHO_SULFATE)    !total volume of a single particle
       else
         vliq = 0._f
       end if

       vfull = volfull(rmon, alpha)

       if (ibin.le.6) then
        write(*,*) 'ibin, vol_sulf',ibin, mmr(1,2,ibin)  / RHO_SULFATE 
        write(*,*) 'ibin, vol_alum',ibin, mmr(1,3,ibin)  / RHO_ALUMINUM
        write(*,*) 'rmon',rmon
        write(*,*) 'rdry',rdry
        write(*,*) 'vol_alu_mon1',vol_alu_mon1
        write(*,*) 'vliq',vliq
        write(*,*) 'vfull',vfull
       end if


 


       ! We need the mass of H2SO4/H2O around alumina for the contact angle calculation.
       ! For this we use the ratio between dry and wet radius.

      ! vol_tot = (rdry**3/rmon_b**3) * vol_alu_mon1  !total volume of a single particle assuming rmon as radius , scaled to alumina volume

       !vol_tot = (r_wet(ibin,2)**3/rdry**3) * vol_alu_mon1  !total volume of a single particle assuming rmon as radius , scaled to alumina volume



       ! calculate the frac_wet alum plus sulfate  to the dry alumina radius to get the total volume per monomer
      ! vol_tot = (r_wet(ibin,2)**3/rmon**3) * vol_alu_mon1  !total volume of a single particle assuming rmon as dry radius

      ! vliq = vol_tot - vol_alu_mon1  ! volume of the particle without alumina (therefore H2SO4/H2O)

      ! ad_alusulf  = 4._f * PI * numberDensity(ibin,2) * (r_wet(ibin,2))**2._f *100000000._f  ! #/cm3 * r_alu surface area density 

       !Calculate volume threshold when particles is immerged by sulfates/H2O

      ! vfull = volfull(rmon, alpha)


       if (vliq.gt.0._f) then
 
         if (vliq.ge.vfull) then
            ! Particle is fully immersed: r_wet already describes the whole
            ! wetted aggregate, so no per-monomer -> per-particle scaling is needed.
             area_h2so4  = 4._f * PI * numberDensity(ibin,2) * (r_wet(ibin,2))**2._f *100000000._f  ! #/cm3 * r_alu surface area density
             area_solid=0._f
             adh2so4_bin(1,I_GRP_ALSUL,ibin) = area_h2so4
             adsolid_bin(1,I_GRP_ALSUL,ibin) = area_solid
             write(*,*) 'immersion: bin, area_h2so4, area_solid ',ibin, area_h2so4, area_solid
         else
           write(*,*) 'before cangle'

          call cangle(alpha,vliq,rmon,area_h2so4,area_solid,hliq,rliq)
           write(*,*) 'start cangle'

          ! area_h2so4/area_solid from cangle() are per-monomer areas, so scale by
          ! both the particle number density and the number of monomers per
          ! aggregate (nmon) to get the bin's total surface area density.
          area_h2so4_out = area_h2so4*100000000._f * numberDensity(ibin,2) * nmon(ibin,I_GRP_ALSUL)
          area_solid_out = area_solid*100000000._f * numberDensity(ibin,2) * nmon(ibin,I_GRP_ALSUL)
          adh2so4_bin(1,I_GRP_ALSUL,ibin) = area_h2so4_out
          adsolid_bin(1,I_GRP_ALSUL,ibin) = area_solid_out
          write(*,*) 'bin, area_h2so4, area_solid ',ibin, area_h2so4_out, area_solid_out

         end if

        end if


    end do

    ! Writing output in NETCDF
    write(*,*) 'start CARMA diag'

    call carmadiags(cstate, carma, NZ, NBIN, NELEM, NGROUP, NWAVE, &
                deltaz, &
                nd=nd(:, :), ad=ad(:, :), md=md(:, :), re=re(:, :), rew=rew(:, :), &
                rm=rm(:, :), jn=jn(:, :), mr=mr(:, :), pa=pa(:, :), &
                ar=ar(:, :), vm=vm(:, :), ex=ex(:,:), od=od(:,:), &
                wr_bin=wr_bin(:,:,:), nd_bin=nd_bin(:,:,:), &
                ro_bin=ro_bin(:,:,:), mr_bin=mr_bin(:,:,:), &
                vd_bin=vd_bin(:,:,:), mr_elm_bin=mr_elm_bin(:,:,:))

    write(*,*) 'start nc write'
    call ncwrt_ilp(ncflgs, outid, NZ, NY, NX, NBIN, NGAS, NELEM, NGROUP, istep,  &
                zc, iy, ix, bins, time, secsperday, date, &
                sname_gas, sname_elm, sname_grp, &
                p, t, zc, deltaT, rhoa, rlheat, &
                mmr=mmr_gas, si=satice, sl=satliq, ei=ei, el=el, wt=wt, &
                nd=nd, ad=ad, md=md, re=re, rew=rew, rm=rm, jn=jn, mr=mr, pa=pa, &
                ar=ar, vm=vm, ex_vis=ex, od_vis=od,&
                vd_bin=vd_bin, nd_bin=nd_bin, wr_bin=wr_bin, ro_bin=ro_bin, mr_bin=mr_bin, &
                mr_elm_bin=mr_elm_bin, &
                adh2so4_bin=adh2so4_bin, adsolid_bin=adsolid_bin, &
                corefrac_bin=corefrac_bin, shellthick_bin=shellthick_bin)

    write(*,*) 'finished nc write'
  end do   ! time loop
 
  ! Close the output NETCDF file
    write(*,*) 'before ncclose'
  call ncclose(outid)       
  if (rc /=0) write(*, *) "    *** FAILED ***, rc=", rc

    write(*,*) 'after ncclose'
  ! Cleanup the carma state objects
  call CARMASTATE_Destroy(cstate, rc)
  if (rc /=0) stop "    *** CARMASTATE_Destroy FAILED ***"
    write(*,*) 'after carmastate_destry'
  
  call CARMA_Destroy(carma, rc)
  if (rc /=0) stop "    *** CARMA_Destroy FAILED ***"
    write(*,*) 'after carma_destry'
end subroutine

