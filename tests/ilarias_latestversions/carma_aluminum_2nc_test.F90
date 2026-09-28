!! This code is to test the impact of particle swelling from
!! relative humidity on sedimentation. Upon execution, a text
!! file (carma_swelltest.txt) is generated.  The text file can
!! be read with the IDL procedure read_swelltest.pro.
!!

program carma_aluminum_2nc_test
  implicit none

  write(*,*) "Aluminum Test"

  call test_aluminum_simple()
end program


subroutine check(status)
  use netcdf

  implicit none
  integer, intent ( in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
end subroutine check

!! Just have one grid box. In that grid box, but an initial concentration
!! of drops at the smallest size, then allow that to grow using a gas. The
!! total mas of drops + gas should be conserved.
subroutine test_aluminum_simple()
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
  use netcdf
  implicit none

  ! Model dimension 
  integer, parameter        :: NZ           = 3
  integer, parameter        :: NY           = 1 ! lat
  integer, parameter        :: NX           = 1 ! lon
  integer, parameter        :: NZP1         = NZ+1
  integer, parameter        :: NELEM        = 1
  integer, parameter        :: NGROUP       = 1
  integer, parameter        :: NBIN         = 8
  integer, parameter        :: NSOLUTE      = 0
  integer, parameter        :: NGAS         = 0
  integer, parameter        :: NWAVE        = 30 !22

  ! Model Grid
  real(kind=f), parameter       :: lat(NY) = (/0.0_f/) 
  real(kind=f), parameter       :: lon(NX) = (/30.0_f/)
  real(kind=f)                  :: time 
  integer                       :: bins(NBIN) 

  ! Different sizes for time steps provide different results
  ! because of the satbility issues.
  real(kind=f), parameter   :: dtime  = 1800._f
  real(kind=f), parameter   :: dArea  = 56000000_f! For regional 58561162288._f  ! gridbox size [m2]
  real(kind=f), parameter   :: deltaz = 0.56_f! For regional 1000._f         ! gridbox height [m]
  real(kind=f), parameter   :: zmin   = 19500._f
  
  integer, parameter        :: nstep  = 672 !240 
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
  integer               :: iwave

  ! CARMA object
  type(carma_type), target       :: carma
  type(carma_type), pointer      :: carma_ptr
  type(carmastate_type)          :: cstate
  integer                        :: rc = 0

  ! Input for Group/element
  integer, parameter        :: I_ALUMINUM  = 1
  integer, parameter        :: I_GRP_ALUM  = 1
  integer, parameter        :: I_ELEM_ALUM = 1

  real(kind=f)              :: rmrat = 2._f
  real(kind=f)              :: rmin = 21.5e-6_f
  real(kind=f)              :: rmon = 21.5e-6_f 
  real(kind=f)              :: df(NBIN,NGROUP) = 1.6_f   ! satellite aerosol fractal dimension, aluminum oxide, Weisenstein 2015, Karasev et al., 2001, 2004
  real(kind=f)              :: falpha = 1._f             ! satellite aerosol fractal packing coefficient, UPDATED VALUE NEEDED!!!
  real(kind=f)              :: RHO_ALUMINUM = 3.95_f     ! dry density of sulfate particles (g/cm3)
  real(kind=f)              :: vf_const = 0.01919_f
  real(kind=f)              :: WTMOL_SOLID = 101.961_f   ! molecular weight (g/mol)
 
  ! Input for optics
  character(len=100)     :: file_optics = "/glade/u/home/iquaglia/carma_box/run/carma/alumina_215nm_avg_fractaloptics_test.nc"
!  character(len=100)     :: file_optics = "/glade/u/home/iquaglia/carma_box/run/carma/alumina_215nm_avg_SOCOL_fractaloptics_test.nc"  
  !character(len=100)     :: file_optics = "/glade/u/home/iquaglia/carma_box/run/carma/Alumina_SOCOL.nc"
  integer               :: idx_vis = 4 !3
  integer               :: inid, wave_varid, bin_varid, qext_varid
  real(kind=f)          :: wave(NWAVE)
  real(kind=f)          :: qext(NWAVE, NBIN, NGROUP)
  
  
  ! Input for SAI
  real(kind=f) :: inj_rate = 2809.28_f * 1e10_f  ! Injection rate [molecules/cm3/s] (previous value: 1831.16_f)
  real(kind=f) :: f_inj(nstep+1, NZ, NGROUP, NBIN) = 0.0_f  ! Time-and-space-dependent scaling factor applied to inj_rate. If zero, it means no injection

  ! Netcdf variable
  type(nc_type)      :: ncflgs
  character(len=80)  :: filename_out = "alumina-215nm_30min_vf_5Tg_plume_test.nc"
  integer            :: outid 
  integer            :: date = 00010101
  real(kind=f)       :: fill_value = -9.255963134931783e+061_f
 
  ! Output for Atmospheric variables
  real(kind=f) :: zc(NZ), zl(NZP1), p(NZ), pl(NZP1), t(NZ), rhoa(NZ), rlheat(NZ), deltaT(NZ) 

  ! Output for Gas variables
  real(kind=f) :: mmr_gas(NZ,NGAS), satliq(NZ,NGAS), satice(NZ,NGAS)
  real(kind=f) :: ei(NZ,NGAS), el(NZ,NGAS), wt(NZ,NGAS)
  
  ! Output for Group variables
  real(kind=f) :: r(NBIN, NGROUP), rlow(NBIN,NGROUP), rup(NBIN,NGROUP), dr(NBIN,NGROUP)
  real(kind=f) :: arat(NBIN, NGROUP), rrat(NBIN, NGROUP), rprat(NBIN, NGROUP) 
  real(kind=f) :: rmass(NBIN,NGROUP), dm(NBIN,NGROUP), vol(NBIN,NGROUP)
  real(kind=f) :: nmon(NBIN,NGROUP) 
    
  real(kind=f) :: nd(NZ, NGROUP), ad(NZ, NGROUP), md(NZ, NGROUP), re(NZ, NGROUP), rew(NZ, NGROUP)
  real(kind=f) :: rm(NZ, NGROUP), jn(NZ, NGROUP), mr(NZ, NGROUP), pa(NZ, NGROUP)
  real(kind=f) :: ar(NZ, NGROUP), vm(NZ, NGROUP), ex(NZ, NGROUP), od(NZ, NGROUP)
  real(kind=f) :: wr_bin(NZ, NGROUP, NBIN), nd_bin(NZ, NGROUP, NBIN), ro_bin(NZ, NGROUP, NBIN)
  real(kind=f) :: tm_bin(NZ, NGROUP, NBIN), mr_bin(NZ, NGROUP, NBIN), vd_bin(NZ, NGROUP, NBIN)

  ! Working variables 
  character(len=80) :: sname,sname_gas(NGAS),sname_elm(NELEM), sname_grp(NGROUP)  
  real(kind=f)      :: t_orig(NZ)
  real(kind=f)      :: new_gas(NZ, NGAS) ! If you want to increment mmr_gas
  real(kind=f)      :: mmr(NZ, NELEM, NBIN)  
  real(kind=f)      :: mmr_inj(NZ, NELEM, NBIN) = 0._f   ! injection rate for SAI studies in (kg/kg)  
  integer           :: val
  logical           :: inject
  
  ! Open netCDF file  
  CALL check(nf90_open(file_optics, nf90_nowrite, inid))
 
  ! Get the varids of the wavelength, bin and qext index
  call check( nf90_inq_varid(inid, "wave", wave_varid) ) 
  call check( nf90_inq_varid(inid, "bin", bin_varid) ) 
  call check( nf90_inq_varid(inid, "fractal_qext", qext_varid) )

  ! Get the values
  call check(nf90_get_var(inid,wave_varid,wave))
  call check(nf90_get_var(inid,qext_varid,qext(:,:,1)))

  ! Close the file
  call check( nf90_close(inid) ) 
 
  ! Define the particle-grid extent of the CARMA test
  call CARMA_Create(carma, NBIN, NELEM, NGROUP, NSOLUTE, NGAS, NWAVE, rc)
  if (rc /=0) stop "    *** CARMA_Create FAILED ***"

  carma_ptr => carma

  ! Define Groups
  call CARMAGROUP_Create(carma, I_GRP_ALUM, "pure_aluminum", rmin, rmrat, &
                        I_SPHERE, 1._f, .false., rc,&
                        is_fractal=.TRUE., rmon=rmon, df=df, falpha=falpha, &
                        irhswell=I_NO_SWELLING, do_drydep=.true., & 
                        shortname="PRALUM", is_sulfate=.false.) ! solfac=0.1_f, scavcoef=0.1_f, 
  if (rc /=0) stop "    *** CARMAGROUP_Create FAILED ***"

  ! Define the elements
  call CARMAELEMENT_Create(carma, I_ELEM_ALUM, I_GRP_ALUM, "Aluminum", RHO_ALUMINUM, I_INVOLATILE, I_ALUMINUM, rc, shortname="ALUM")
  
  if (rc /=0) stop "    *** CARMAELEMENT_Create FAILED ***"


  ! Setup the CARMA processes
  call CARMA_AddCoagulation(carma, I_GRP_ALUM, I_GRP_ALUM, I_GRP_ALUM, I_COLLEC_FUCHS, rc)
  if (rc /=0) stop "    *** CARMA_AddCoagulation FAILED ***"


  call CARMA_Initialize(carma, rc, do_grow=.false., do_coag=.true., do_substep=.false., do_vtran=.TRUE.) !, vf_const=vf_const) 

  if (rc /=0) stop "    *** CARMA_Initialize FAILED ***"

  ! For simplicity of setup, do a case with Cartesian coordinates,
  ! which are specified in this interface in meters.
  !
  ! NOTE: For Cartesian coordinates, the first level is the bottom
  ! of the model (e.g. z = 0), while for sigma and hybrid coordinates
  ! the first level is the top of the model.
 
  ! Vertical center
  do i = 1, NZ
        zc(i) = zmin + deltaz * (i - 0.5_f)
  end do
  call GetStandardAtmosphere(zc, p=p, t=t)
 
  !! For global model comparison, we want to redefine the pressure and temperature at a specific lat and lon
  p(2)    = 52._f * 100._f  ! Redifine pressure of selected level    [Pa]
  t(2)    = 206 ! 202.5 !209._f          ! Redifine temperature of selected level [K]
  rhoa(2) = p(2) / (R_AIR * 1e-4_f * t(2)) ! kg/m3   
 
  do i = 1, NZ
       p(i) = p(2) * exp(-GRAV / 100._f * (zc(i) - zc(2))/(R_AIR * 1e-4_f * t(2)))
       t(i) = t(2) - (-0.001_f) * (zc(i) - zc(2))
  end do

  do i = 1, NZ
        rhoa(i) = p(i) / (R_AIR * 1e-4_f * t(i)) ! kg/m3  
  end do

 
  ! Vertical edge
  do i = 1, NZP1
        zl(i) = zmin + (i - 1) * deltaz
  end do
  call GetStandardAtmosphere(zl, p=pl)
 
  !! For global model comparison, we want to redefine the pressure and temperature at a specific lat and lon
  ! rewrite it in reference to level 2
  do i = 1, NZP1 -1 
          pl(i) = p(i)  - (zl(i) - zc(i)) * rhoa(i) * (GRAV / 100._f)
  end do 
  i = NZP1
  pl(i) = p(i-1) - (zl(i) - zc(i-1)) * rhoa(i-1) * (GRAV / 100._f) 


  ! T and rho_air
  t_orig = t
 
  
  ! Bin structure
    do ibin = 1,NBIN
        bins(ibin) = ibin
  end do


  ! Initial concentrations   
  mmr(:,:,:) = 0._f ! Initial element concentration
  
  ! Injection of more mass (SAI) 
  ! mmr_inj(:,:,:) = 0._f ! Initialization of mmr_inj
  mmr_inj(2,1,1) = inj_rate * 1e6_f * dtime * WTMOL_SOLID * 1e-3_f / AVG /rhoa(2) 
  
  do istep = 2, nstep
        val = modulo(istep - 49, 48)
        if (istep <= 5 .or. (val >= 0 .and. val < 4)) then
                f_inj(istep, 2, 1, 1) = 1._f
        end if
  end do

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
           ncflgs%f_nc_group(igroup)%f_nc_mr    = .TRUE.
           ncflgs%f_nc_group(igroup)%f_nc_pa    = .TRUE.
           ncflgs%f_nc_group(igroup)%f_nc_ar    = .TRUE.
           ncflgs%f_nc_group(igroup)%f_nc_vm    = .TRUE. 
           ncflgs%f_nc_group(igroup)%f_nc_exvis = .TRUE.
           ncflgs%f_nc_group(igroup)%f_nc_odvis = .TRUE.

           ncflgs%f_nc_group(igroup)%f_nc_wr_bin = .TRUE.
           ncflgs%f_nc_group(igroup)%f_nc_nd_bin = .TRUE.
           ncflgs%f_nc_group(igroup)%f_nc_ro_bin = .TRUE.
           ncflgs%f_nc_group(igroup)%f_nc_tm_bin = .TRUE.
           ncflgs%f_nc_group(igroup)%f_nc_mr_bin = .TRUE.
           ncflgs%f_nc_group(igroup)%f_nc_vd_bin = .TRUE.
  end do

  
  ! Get variables for nc definition
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
  call ncdef(ncflgs, filename_out, fill_value, NZ, NY, NX, NBIN, NGAS, NELEM, NGROUP, &
         sname_grp=sname_grp(:), sname_elm=sname_elm(:),&
         outid=outid)

  ! Fill nc variables that are not time dependent (out of the loop) 
  call ncwrt_olp(ncflgs, outid, NZ, NY, NX, NGAS, NELEM, NGROUP, NBIN, &  
            p(:), lat, lon, bins, sname_grp, &
          r=r(:,:), rmass=rmass(:,:), vol=vol(:,:), &
          dr=dr(:,:), dm=dm(:,:), rup=rup(:,:), rlow=rlow(:,:), &
          rrat=rrat(:,:), arat=arat(:,:),rprat=rprat(:,:), &
          df=df(:,:), nmon=nmon(:,:))
  ! ---------------------------------------------------------------------------------

  !  ------------------ Iterate the model over time steps   -------------------------
  do istep = 2, nstep + 1

    ! Calculate the model time.
    time = (istep - 1) * dtime

    ! Create a CARMASTATE for this column. 
    call CARMASTATE_Create(cstate, carma_ptr, time, dtime, NZ, &
                        I_CART, lat(iy), lon(ix), &
                        zc(:), zl(:), &
                        p(:),  pl(:), &
                        t(:), rc, &
                        told=t(:))

    if (rc /=0) stop "    *** CARMASTATE_Create FAILED ***"

    deltaT(:) = t(:) - t_orig(:)

    
    if (istep == 2) then
        
        call carmadiags(cstate, carma, NZ, NBIN, NELEM, NGROUP, NWAVE, &
                dArea, deltaz, qext=qext(:,:,:), idx_wave=idx_vis, &
                nd=nd(:, :), ad=ad(:, :), md=md(:, :), re=re(:, :), rew=rew(:, :), &
                rm=rm(:, :), mr=mr(:, :), pa=pa(:, :), &
                ar=ar(:, :), vm=vm(:, :), ex=ex(:,:), od=od(:,:), &
                wr_bin=wr_bin(:,:,:), nd_bin=nd_bin(:,:,:), &
                ro_bin=ro_bin(:,:,:), tm_bin=tm_bin(:,:,:), mr_bin=mr_bin(:,:,:), &
                vd_bin=vd_bin(:,:,:) )

        call ncwrt_ilp(ncflgs, outid, NZ, NY, NX, NBIN, NGAS, NELEM, NGROUP, 1,  &
                zc, iy, ix, bins, 0._f, secsperday, date, &
                sname_grp=sname_grp, sname_elm=sname_elm(:), &
                p=p, t=t, z=zc, dT=deltaT, rhoa=rhoa, rlheat=rlheat, & 
                nd=nd, ad=ad, md=md, re=re, rew=rew, rm=rm, mr=mr, pa=pa, &
                ar=ar, vm=vm, ex_vis=ex, od_vis=od, &
                vd_bin=vd_bin, nd_bin=nd_bin, wr_bin=wr_bin,& 
                ro_bin=ro_bin, tm_bin=tm_bin, mr_bin=mr_bin)
    end if


    ! Add more mass
    mmr(:,:,:) = mmr(:,:,:) + mmr_inj(:,:,:) * f_inj(istep,:,:,:)  
 

    ! Send the bin mmrs to CARMA
    do ielem = 1, NELEM
      do ibin = 1, NBIN
        call CARMASTATE_SetBin(cstate, ielem, ibin, mmr(:,ielem,ibin), rc)
        if (rc /=0) stop "    *** CARMASTATE_SetBin FAILED ***"
      end do
    end do


    ! Execute the step
    call CARMASTATE_Step(cstate, rc)
    if (rc /=0) stop "    *** CARMASTATE_Step FAILED ***"


    ! Get the retry stats and the updated temperature.
    call CARMASTATE_Get(cstate, rc, nsubstep=nsubsteps, nretry=nretries)
    if (rc /=0) stop "    *** CARMASTATE_Get FAILED ***"

    call CARMASTATE_GetState(cstate, rc, t=t(:), rhoa_wet=rhoa(:)) 
    if (rc /=0) stop "    *** CARMASTATE_GetState FAILED ***" 
    
    ! Get the updated bin mmr. 
    do ielem = 1, NELEM
      do ibin = 1, NBIN
        call CARMASTATE_GetBin(cstate, ielem, ibin, mmr(:,ielem,ibin), rc)
        if (rc /=0) stop "    *** CARMASTATE_GetBin FAILED ***"
      end do
    end do

    call carmadiags(cstate, carma, NZ, NBIN, NELEM, NGROUP, NWAVE, &
                dArea, deltaz, qext=qext(:,:,:), idx_wave=idx_vis, &
                nd=nd(:, :), ad=ad(:, :), md=md(:, :), re=re(:, :), rew=rew(:, :), &
                rm=rm(:, :), mr=mr(:, :), pa=pa(:, :), &
                ar=ar(:, :), vm=vm(:, :), ex=ex(:,:), od=od(:,:), &
                wr_bin=wr_bin(:,:,:), nd_bin=nd_bin(:,:,:), &
                ro_bin=ro_bin(:,:,:), tm_bin=tm_bin(:,:,:), mr_bin=mr_bin(:,:,:), &
                vd_bin=vd_bin(:,:,:) )
    
    
    call ncwrt_ilp(ncflgs, outid, NZ, NY, NX, NBIN, NGAS, NELEM, NGROUP, istep,  &
                zc, iy, ix, bins, time, secsperday, date, &
                sname_grp=sname_grp, sname_elm=sname_elm(:), &
                p=p, t=t, z=zc, dT=deltaT, rhoa=rhoa, rlheat=rlheat, &
                nd=nd, ad=ad, md=md, re=re, rew=rew, rm=rm, mr=mr, pa=pa, &
                ar=ar, vm=vm, ex_vis=ex, od_vis=od, &
                vd_bin=vd_bin, nd_bin=nd_bin, wr_bin=wr_bin, &
                ro_bin=ro_bin, tm_bin=tm_bin, mr_bin=mr_bin)

  end do   ! time loop
  ! --------------------------------------------------------------------------------

  ! Close the output NETCDF file
  call ncclose(outid)
  if (rc /=0) write(*, *) "    *** FAILED ***, rc=", rc

  ! Cleanup the carma state objects
  call CARMASTATE_Destroy(cstate, rc)
  if (rc /=0) stop "    *** CARMASTATE_Destroy FAILED ***"

  call CARMA_Destroy(carma, rc)
  if (rc /=0) stop "    *** CARMA_Destroy FAILED ***"
end subroutine

