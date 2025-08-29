!! This code is to test the impact of particle swelling from
!! relative humidity on sedimentation. Upon execution, a text
!! file (carma_swelltest.txt) is generated.  The text file can
!! be read with the IDL procedure read_swelltest.pro.
!!

program carma_sulfate_2nc_test
  implicit none

  write(*,*) "Sulfate Test"

  call test_sulfate_simple()
end program carma_sulfate_2nc_test

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
  use nc_types_mod
  use test2nc_mod
  use carmadiags_mod
  implicit none

  ! Model dimension
  integer, parameter        :: NZ           = 1
  integer, parameter        :: NY           = 1 ! lat
  integer, parameter        :: NX           = 1 ! lon
  integer, parameter        :: NZP1         = NZ+1
  integer, parameter        :: NELEM        = 1
  integer, parameter        :: NGROUP       = 1
  integer, parameter        :: NBIN         = 22
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
  integer, parameter        :: nstep  = 180000 / int(dtime)
  real(kind=f), parameter   :: secsperday = 24._f * 3600._f

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
  integer, parameter        :: I_H2SO4  = 1
  real(kind=f)              :: rmin = 2.e-8_f
  real(kind=f)              :: rmrat = 2._f
  real(kind=f)              :: RHO_SULFATE = 1.923_f  ! dry density of sulfate particles (g/cm3)

  ! Output for Atmospheric variables
  real(kind=f) :: zc(NZ), zl(NZP1), p(NZ), pl(NZP1), t(NZ), rhoa(NZ), rlheat(NZ), deltaT(NZ)

  ! Output for Gas variables
  real(kind=f) :: mmr_gas(NZ,NGAS), satliq(NZ,NGAS), satice(NZ,NGAS)
  real(kind=f) :: ei(NZ,NGAS), el(NZ,NGAS), wt(NZ,NGAS)

  ! Output for Group variables
  real(kind=f) :: r(NBIN, NGROUP), rlow(NBIN,NGROUP), rup(NBIN,NGROUP), dr(NBIN,NGROUP)
  real(kind=f) :: arat(NBIN, NGROUP), rrat(NBIN, NGROUP), rprat(NBIN, NGROUP)
  real(kind=f) :: rmass(NBIN,NGROUP), dm(NBIN,NGROUP), vol(NBIN,NGROUP)
  real(kind=f) :: df(NBIN,NGROUP), nmon(NBIN,NGROUP)

  real(kind=f) :: nd(NZ, NGROUP), ad(NZ, NGROUP), md(NZ, NGROUP), re(NZ, NGROUP), rew(NZ, NGROUP)
  real(kind=f) :: rm(NZ, NGROUP), jn(NZ, NGROUP), mr(NZ, NGROUP), pa(NZ, NGROUP)
  real(kind=f) :: ar(NZ, NGROUP), vm(NZ, NGROUP), ex(NZ, NGROUP), od(NZ, NGROUP)
  real(kind=f) :: wr_bin(NZ, NGROUP, NBIN), nd_bin(NZ, NGROUP, NBIN), ro_bin(NZ, NGROUP, NBIN)
  real(kind=f) :: mr_bin(NZ, NGROUP, NBIN), vd_bin(NZ, NGROUP, NBIN)

  ! Output for Element variables
  real(kind=f) :: mr_elm_bin(NZ, NELEM, NBIN)


  ! Netcdf variable
  type(nc_type)     :: ncflgs
  character(len=80) :: filename_out = "sulfatetest_test.nc"
  integer           :: outid
  integer           :: date = 00010101
  real(kind=f)      :: fill_value =-9.255963134931783e+061_f

  ! Working variables
  character(len=80) :: sname,sname_gas(NGAS),sname_elm(NELEM), sname_grp(NGROUP)
  real(kind=f)      :: t_orig(NZ)
  real(kind=f)      :: new_gas(NZ, NGAS) ! If you want to increment mmr_gas
  real(kind=f)      :: mmr(NZ, NELEM, NBIN)
  real(kind=f)      :: qext(NWAVE, NBIN, NGROUP)


  ! Define the particle-grid extent of the CARMA test
  call CARMA_Create(carma, NBIN, NELEM, NGROUP, NSOLUTE, NGAS, NWAVE, rc, &
      LUNOPRT=LUNOPRT)
  if (rc /=0) stop "    *** CARMA_Create FAILED ***"

  carma_ptr => carma

  ! Create Groups
  call CARMAGROUP_Create(carma, 1, "sulfate", rmin, rmrat, I_SPHERE, 1._f, .false., &
          rc, irhswell=I_WTPCT_H2SO4, do_drydep=.true., do_mie=.false., &
                        shortname="PRSUL", is_sulfate=.true.)
  if (rc /=0) stop "    *** CARMAGROUP_Create FAILED ***"

  ! Define the elements
  call CARMAELEMENT_Create(carma, 1, 1, "Sulfate", RHO_SULFATE, I_VOLATILE, I_H2SO4, rc, shortname="SULF")
  if (rc /=0) stop "    *** CARMAELEMENT_Create FAILED ***"

  ! Define the gases
  call CARMAGAS_Create(carma, 1, "Water Vapor", WTMOL_H2O, I_VAPRTN_H2O_MURPHY2005, &
    I_GCOMP_H2O, rc, shortname = "Q", dgc_threshold=0.1_f, ds_threshold=0.1_f)
  if (rc /=0) stop "    *** CARMAGAS_Create FAILED ***"

  call CARMAGAS_Create(carma, 2, "Sulphuric Acid", 98.078479_f, I_VAPRTN_H2SO4_AYERS1980, &
    I_GCOMP_H2SO4, rc, shortname = "H2SO4", dgc_threshold=0.1_f, ds_threshold=0.1_f)
  if (rc /=0) stop "    *** CARMAGAS_Create FAILED ***"


  ! Setup the CARMA processes
  call CARMA_AddGrowth(carma, 1, 2, rc)   ! set H2SO4 to be the condensing gas
  if (rc /=0) stop "    *** CARMA_AddGrowth FAILED ***"

   call CARMA_AddNucleation(carma, 1, 1, I_HOMNUC, 0._f, rc, igas=2)
  if (rc /=0) stop "    *** CARMA_AddNucleation FAILED ***"

   call CARMA_AddCoagulation(carma, 1, 1, 1, I_COLLEC_FUCHS, rc)
  if (rc /=0) stop "    *** CARMA_AddCoagulation FAILED ***"


  call CARMA_Initialize(carma, rc, do_grow=.true., do_coag=.true., do_substep=.true., &
          do_thermo=.true., maxretries=16, maxsubsteps=32, dt_threshold=1._f, sulfnucl_method='ZhaoTurco')
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
   p(1)         = 90._f * 100._f ! hPa
   zc(1)        = 17000._f       ! m
   t(1)         = 250._f         ! K
   t_orig = t(1)
   zl(1)        = zc(1) - deltaz
   zl(2)        = zc(1) + deltaz
   rhoa(1)       = (p(1) * 10._f) / (R_AIR * t(1)) * (1e-3_f * 1e6_f)
   pl(1)        = p(1) - (zl(1) - zc(1)) * rhoa(1) * (GRAV / 100._f)
   pl(2)        = p(1) - (zl(2) - zc(1)) * rhoa(1) * (GRAV / 100._f)

   ! Initial H2O and H2SO4 concentrations
   mmr_gas(:,1)  = 100.e-6_f     ! H2O g/g
   mmr_gas(:,2)  = 0.1e-9_f * (98._f / 29._f)    ! H2SO4
   mmr(:,:,:)   = 0._f ! Initial element concentration

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

    ! Get the updated gas mmr.
    do igas = 1, NGAS
      call CARMASTATE_GetGas(cstate, igas, &
                             mmr_gas(:,igas), rc, &
                             satliq=satliq(:,igas), satice=satice(:,igas), &
                             eqice=ei(:,igas), eqliq=el(:,igas), wtpct=wt(:,igas))
      if (rc /=0) stop "    *** CARMASTATE_GetGas FAILED ***"
    end do

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
