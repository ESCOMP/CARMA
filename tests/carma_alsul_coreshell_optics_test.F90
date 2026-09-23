!! This test compares different ways of representing the optical properties
!! of the ALSUL group (a mixed alumina-core / sulfate-and-water-shell aerosol)
!! as a function of composition, all using compact-sphere Mie theory
!! (I_MIERTN_TOON1981), so that geometry (sphere vs. fractal) is NOT a
!! confounding variable here - that comparison is done separately in
!! carma_alsul_fractal_optics_test.F90.
!!
!! Composition-mixing methods compared (see carma_alsul_optics_common_mod.F90
!! for the shared refractive index data and the Maxwell-Garnett helper):
!!   ALSUL_fixed_*             - today's production behavior (I_OPTICS_FIXED):
!!                               a single, alumina-only refractive index is
!!                               used for the whole particle, regardless of how
!!                               much sulfate/water coating is actually
!!                               present. Radius still swells with RH.
!!   ALSUL_coreshell_cfNN_*    - literal two-layer core/shell Mie (miess): an
!!                               alumina core plus a sulfate+water shell with
!!                               its own volume-mixed refractive index.
!!   ALSUL_maxwell_cfNN_*      - Maxwell-Garnett effective-medium mixing:
!!                               alumina treated as inclusions within the
!!                               sulfate+water shell, producing a single
!!                               effective refractive index for an equivalent
!!                               homogeneous sphere.
!!
!! "cfNN" is the assumed dry-particle alumina (core) volume fraction (e.g.
!! cf90 = 90% alumina by volume, a thin sulfate coating; cf30 = 30% alumina,
!! a thick coating), swept over core_frac_vals below so that the comparison
!! shows how each method's optics respond to the AMOUNT of sulfate coating,
!! not just to whether a coating is present at all. "_dry"/"_wet" is whether
!! the sulfate has additionally taken up water at rh_wet.
!!
!! A plain PRSUL (pure sulfate) group is included as a stable reference.
!! Each "run" is written to carma_alsul_coreshell_optics_test.nc via
!! test2nc_mod's ncdef_optics/ncwrt_optics, following the same NRUN pattern
!! as carma_fractaloptics_2nc_test.F90.
!!
!! Upon execution, a netcdf file (carma_alsul_coreshell_optics_test.nc) is
!! generated, to be plotted similarly to ./tests/python/plots_optics.ipynb.

program carma_alsul_coreshell_optics_test
  implicit none

  write(*,*) "ALSUL Core/Shell Optics Comparison Test"

  call test_alsul_coreshell_optics()

  write(*,*) "Done"
end program carma_alsul_coreshell_optics_test


subroutine test_alsul_coreshell_optics()
  use carma_precision_mod
  use carma_constants_mod
  use carma_enums_mod
  use carma_types_mod
  use carmaelement_mod
  use carmagroup_mod
  use carma_mod
  use test2nc_mod
  use carma_alsul_optics_common_mod
  use wetr, only: getwetr

  implicit none

  ! Model dimensions.
  integer, parameter :: NELEM   = 3
  integer, parameter :: NBIN    = 20
  integer, parameter :: NGROUP  = 2
  integer, parameter :: NSOLUTE = 0
  integer, parameter :: NGAS    = 0
  integer, parameter :: NREFIDX = 1

  ! Use one assumed dry-particle alumina (core) volume fraction initially.
  integer, parameter      :: NCOREFRAC = 4
  real(kind=f), parameter :: core_frac_vals(NCOREFRAC) = (/ 0.9_f, 0.7_f, 0.5_f, 0.3_f /)

  ! 1 PRSUL reference + 2 ALSUL_fixed (dry/wet) + coreshell (dry/wet)
  integer, parameter :: NRUN = 1 + 2 + 2 * NCOREFRAC

  ! CARMA object.
  type(carma_type), target  :: carma
  type(carma_type), pointer :: carma_ptr
  integer                   :: rc = 0

  ! Group/element indices.
  integer, parameter :: I_GRP_PRSUL  = 1
  integer, parameter :: I_GRP_ALSUL  = 2
  integer, parameter :: I_ELEM_PRSUL = 1
  integer, parameter :: I_ELEM_ALSUL = 2
  integer, parameter :: I_ELEM_ALUM  = 3
  integer, parameter :: I_H2SO4      = 1
  integer, parameter :: I_ALUMINUM   = 2

  integer :: bins(NBIN)
  integer :: ibin, iwave, irun, imethod, icf

  ! Group geometry.
  real(kind=f) :: df(NBIN, NGROUP)
  real(kind=f) :: rhogroup(NBIN)
  real(kind=f) :: r(NBIN,NGROUP), rmass(NBIN,NGROUP), vol(NBIN,NGROUP)
  real(kind=f) :: rlow(NBIN,NGROUP), rup(NBIN,NGROUP), dr(NBIN,NGROUP), dm(NBIN,NGROUP)
  real(kind=f) :: nmon(NBIN,NGROUP)

  ! Refractive indices.
  complex(kind=f) :: refidx_alum(NWAVE_RRTMG)
  complex(kind=f) :: refidx_sulf(NWAVE_RRTMG)
  complex(kind=f) :: refidx_water(NWAVE_RRTMG)
  complex(kind=f) :: m_prsul

  real(kind=f) :: core_frac, rho_dry_mix

  real(kind=f), parameter :: rh_wet     = 0.9_f
  real(kind=f), parameter :: temp_swell = 250._f

  ! Output arrays.
  character(len=80) :: name_run(NRUN)
  character(len=2)  :: cf_label
  real(kind=f)       :: out_r(NBIN,NRUN), out_rmass(NBIN,NRUN), out_vol(NBIN,NRUN)
  real(kind=f)       :: qext(NWAVE_RRTMG,NBIN,NRUN), ssa(NWAVE_RRTMG,NBIN,NRUN), asy(NWAVE_RRTMG,NBIN,NRUN)
  real(kind=f)       :: out_refidx_r(NWAVE_RRTMG,NRUN), out_refidx_i(NWAVE_RRTMG,NRUN)

  ! Working variables.
  real(kind=f) :: rwet_bin, rhopwet
  real(kind=f) :: radius, Qext_l, Qsca_l, asym_l
  logical      :: is_wet, is_coreshell

  integer            :: outid
  character(len=80)  :: filename_out = "carma_SDPSul_CS_optics_test_cf30,50,70,90.nc"

  rc = 0

  ! ---------------- Build the CARMA object (both groups, once) --------------------

  call CARMA_Create(carma, NBIN, NELEM, NGROUP, NSOLUTE, NGAS, NWAVE_RRTMG, rc, &
       wave=wave_rrtmg, NREFIDX=NREFIDX)
  if (rc < 0) stop "    *** CARMA_Create FAILED ***"
  carma_ptr => carma

  call CARMAGROUP_Create(carma, I_GRP_PRSUL, "sulfate", rmin_PRSUL, rmrat_PRSUL, I_SPHERE, 1._f, .false., &
       rc, irhswell=I_WTPCT_H2SO4, do_drydep=.false., is_sulfate=.true., shortname="PRSUL")
  if (rc < 0) stop "    *** CARMAGROUP_Create (PRSUL) FAILED ***"

  df(:,I_GRP_ALSUL)       = df_ALSUL_value
  rhogroup(:)             = RHO_ALUMINUM
  call CARMAGROUP_Create(carma, I_GRP_ALSUL, "alumina", rmin_ALUM, rmrat_ALUM, I_SPHERE, 1._f, .false., &
       rc, irhswell=I_PETTERS, do_drydep=.false., is_sulfate=.true., shortname="ALSUL", &
       is_fractal=.TRUE., rmon=rmon_ALUM, df=df(:,I_GRP_ALSUL), falpha=falpha_ALUM, rho=rhogroup)
  if (rc < 0) stop "    *** CARMAGROUP_Create (ALSUL) FAILED ***"

  call CARMAELEMENT_Create(carma, I_ELEM_PRSUL, I_GRP_PRSUL, "Sulfate", &
       RHO_SULFATE, I_VOLATILE, I_H2SO4, rc, shortname="PRSULF")
  if (rc < 0) stop "    *** CARMAELEMENT_Create (PRSUL) FAILED ***"

  call CARMAELEMENT_Create(carma, I_ELEM_ALSUL, I_GRP_ALSUL, "Sulfate in mixed ALSUL", RHO_SULFATE, &
       I_VOLATILE, I_H2SO4, rc, kappa=Kappa_SULF, shortname="ALSULF")
  if (rc < 0) stop "    *** CARMAELEMENT_Create (ALSULF) FAILED ***"

  call CARMAELEMENT_Create(carma, I_ELEM_ALUM, I_GRP_ALSUL, "Aluminum", RHO_ALUMINUM, I_COREMASS, I_ALUMINUM, &
       rc, kappa=Kappa_zero, shortname="ALUM")
  if (rc < 0) stop "    *** CARMAELEMENT_Create (ALUM) FAILED ***"

  call CARMA_Initialize(carma, rc, do_grow=.FALSE., do_coag=.FALSE., do_vtran=.FALSE., do_substep=.FALSE.)
  if (rc < 0) stop "    *** CARMA_Initialize FAILED ***"

  call CARMAGROUP_Get(carma, I_GRP_PRSUL, rc, r=r(:,I_GRP_PRSUL), rlow=rlow(:,I_GRP_PRSUL), rup=rup(:,I_GRP_PRSUL), &
       dr=dr(:,I_GRP_PRSUL), rmass=rmass(:,I_GRP_PRSUL), dm=dm(:,I_GRP_PRSUL), vol=vol(:,I_GRP_PRSUL))
  if (rc < 0) stop "    *** CARMAGROUP_Get (PRSUL) FAILED ***"

  call CARMAGROUP_Get(carma, I_GRP_ALSUL, rc, r=r(:,I_GRP_ALSUL), rlow=rlow(:,I_GRP_ALSUL), rup=rup(:,I_GRP_ALSUL), &
       dr=dr(:,I_GRP_ALSUL), rmass=rmass(:,I_GRP_ALSUL), dm=dm(:,I_GRP_ALSUL), vol=vol(:,I_GRP_ALSUL), &
       nmon=nmon(:,I_GRP_ALSUL))
  if (rc < 0) stop "    *** CARMAGROUP_Get (ALSUL) FAILED ***"

  do ibin = 1, NBIN
    bins(ibin) = ibin
  end do

  ! ---------------- Refractive indices ---------------------------------------------

  call get_alumina_refidx("/home/n02/n02/hqmp0145/carma_input_data/rindex/Alumina.nc", refidx_alum, rc)
  if (rc < 0) stop "    *** get_alumina_refidx FAILED ***"

  refidx_sulf(:)  = cmplx(shellreal(:), shellimag(:), kind=f)
  refidx_water(:) = cmplx(waterreal(:), waterimag(:), kind=f)

  ! ---------------- Run 1: PRSUL reference (dry, homogeneous sulfate sphere) -------

  irun = 1
  name_run(irun) = "PRSUL"
  out_r(:,irun)     = r(:,I_GRP_PRSUL)
  out_rmass(:,irun) = rmass(:,I_GRP_PRSUL)
  out_vol(:,irun)   = vol(:,I_GRP_PRSUL)
  out_refidx_r(:,irun) = real(refidx_sulf(:))
  out_refidx_i(:,irun) = aimag(refidx_sulf(:))


  do ibin = 1, NBIN
    do iwave = 1, NWAVE_RRTMG
      m_prsul = refidx_sulf(iwave)
      call mie(carma, I_MIERTN_TOON1981, r(ibin,I_GRP_PRSUL), wave_rrtmg(iwave), &
               0._f, 3.0_f, 0.0_f, 1.0_f, m_prsul, 0.0_f, m_prsul, Qext_l, Qsca_l, asym_l, rc)
      if (rc < 0) stop "    *** mie (PRSUL) FAILED ***"
      qext(iwave,ibin,irun) = Qext_l
      ssa (iwave,ibin,irun) = Qsca_l / Qext_l
      asy (iwave,ibin,irun) = asym_l
    end do
  end do

  ! ---------------- Runs 2-3: ALSUL_fixed baseline (dry, wet) ----------------------
  !
  ! Today's production behavior: alumina-only refractive index, homogeneous
  ! sphere, but the radius still responds to RH (matches I_OPTICS_FIXED).

  do irun = 2, 3
    is_wet = (irun == 3)
    name_run(irun) = "ALSUL_fixed_" // merge("wet", "dry", is_wet)

    if (is_wet) then
      call getwetr(carma, I_GRP_ALSUL, 1, rh_wet, r(1,I_GRP_ALSUL), rwet_bin, RHO_ALUMINUM, rhopwet, rc, &
           kappa=Kappa_SULF, temp=temp_swell)
      if (rc < 0) stop "    *** getwetr FAILED ***"
    else
      rwet_bin = r(1,I_GRP_ALSUL)
      rhopwet  = RHO_ALUMINUM
    end if

    out_refidx_r(:,irun) = real(refidx_alum(:))
    out_refidx_i(:,irun) = aimag(refidx_alum(:))

    do ibin = 1, NBIN
      radius = r(ibin,I_GRP_ALSUL) * (rwet_bin / r(1,I_GRP_ALSUL))

      do iwave = 1, NWAVE_RRTMG
        call mie(carma, I_MIERTN_TOON1981, radius, wave_rrtmg(iwave), &
                 0._f, 3.0_f, 0.0_f, 1.0_f, refidx_alum(iwave), 0.0_f, refidx_alum(iwave), &
                 Qext_l, Qsca_l, asym_l, rc)
        if (rc < 0) stop "    *** mie (ALSUL_fixed) FAILED ***"
        qext(iwave,ibin,irun) = Qext_l
        ssa (iwave,ibin,irun) = Qsca_l / Qext_l
        asy (iwave,ibin,irun) = asym_l
      end do

      out_r(ibin,irun)     = radius
      out_vol(ibin,irun)   = (4._f/3._f) * PI * radius**3
      out_rmass(ibin,irun) = out_vol(ibin,irun) * rhopwet
    end do
  end do

  ! ---------------- Runs 4+: coreshell and maxwell, swept over coating amount ------
  !
  ! NOTE: for the I_PETTERS swelling law, rwet/rdry is a pure ratio (independent
  ! of the dry radius itself), so wet/dry volume ratios - and therefore m_shell
  ! and f_incl - are the same for every bin in a given run. Only the absolute
  ! radius (rwet_bin, rcore_bin) genuinely varies bin-to-bin, so the shell
  ! mixing and the Maxwell-Garnett inclusion fraction are computed once per
  ! run, using bin 1's dry/wet volume ratio, and reused for every bin.

  irun = 3

  do imethod = 1, 1   ! 1 = coreshell (literal core/shell Mie), 2 = maxwell (effective medium)
    is_coreshell = (imethod == 1)

    do icf = 1, NCOREFRAC
      core_frac   = core_frac_vals(icf)
      rho_dry_mix = core_frac * RHO_ALUMINUM + (1._f - core_frac) * RHO_SULFATE
      write(cf_label, '(I2.2)') nint(core_frac * 100._f)

      ! --- dry run ---
      irun = irun + 1
      if (is_coreshell) then
        name_run(irun) = "ALSUL_coreshell_cf" // cf_label // "_dry"
      else
        name_run(irun) = "ALSUL_maxwell_cf" // cf_label // "_dry"
      end if
      call fill_alsul_run(irun, is_coreshell, .FALSE., core_frac, rho_dry_mix, 0._f)

      ! --- wet run ---
      irun = irun + 1
      if (is_coreshell) then
        name_run(irun) = "ALSUL_coreshell_cf" // cf_label // "_wet"
      else
        name_run(irun) = "ALSUL_maxwell_cf" // cf_label // "_wet"
      end if
      call fill_alsul_run(irun, is_coreshell, .TRUE., core_frac, rho_dry_mix, rh_wet)

    end do
  end do

  ! ---------------- Write the NetCDF output -----------------------------------------

  call ncdef_optics(filename_out, NWAVE_RRTMG, NBIN, NRUN, name_run, outid)

  call ncwrt_optics(outid, NWAVE_RRTMG, NBIN, NRUN, wave_rrtmg, bins, name_run, &
       r_refidx=real(refidx_alum(:)), i_refidx=aimag(refidx_alum(:)), &
       r=out_r, rmass=out_rmass, vol=out_vol, &
       qext=qext, asy=asy, ssa=ssa, &
       r_refidx_run=out_refidx_r, i_refidx_run=out_refidx_i)

  call ncclose(outid)

  call CARMA_Destroy(carma, rc)
  if (rc /=0) stop "    *** CARMA_Destroy FAILED ***"

contains

  !! Computes one ALSUL run (either the coreshell or the maxwell treatment) at
  !! the given core_frac and RH, filling qext/ssa/asy/out_r/out_rmass/out_vol
  !! and out_refidx_r/i for column "irun_l".
  subroutine fill_alsul_run(irun_l, is_coreshell_l, is_wet_l, core_frac_l, rho_dry_mix_l, rh_l)
    integer, intent(in)      :: irun_l
    logical, intent(in)      :: is_coreshell_l, is_wet_l
    real(kind=f), intent(in) :: core_frac_l, rho_dry_mix_l, rh_l

    real(kind=f) :: vol_dry_tot_l, vol_core_l, vol_sulf_dry_l
    real(kind=f) :: rwet_bin_l, rhopwet_l, vol_wet_tot_l, vol_shell_wet_l, vol_water_added_l, f_incl_l
    real(kind=f) :: radius_l, rcore_bin_l
    complex(kind=f) :: m_shell_l, m_eff_l
    integer :: ibin_l, iwave_l, rc_l

    rc_l = 0

    vol_dry_tot_l  = (4._f/3._f) * PI * r(1,I_GRP_ALSUL)**3
    vol_core_l     = core_frac_l * vol_dry_tot_l
    vol_sulf_dry_l = vol_dry_tot_l - vol_core_l

    if (is_wet_l) then
      call getwetr(carma, I_GRP_ALSUL, 1, rh_l, r(1,I_GRP_ALSUL), rwet_bin_l, rho_dry_mix_l, rhopwet_l, rc_l, &
           kappa=Kappa_SULF, temp=temp_swell)
      if (rc_l < 0) stop "    *** getwetr FAILED ***"
    else
      rwet_bin_l = r(1,I_GRP_ALSUL)
      rhopwet_l  = rho_dry_mix_l
    end if

    vol_wet_tot_l     = (4._f/3._f) * PI * rwet_bin_l**3
    vol_shell_wet_l   = vol_wet_tot_l - vol_core_l
    vol_water_added_l = vol_wet_tot_l - vol_dry_tot_l
    f_incl_l          = vol_core_l / vol_wet_tot_l

    do iwave_l = 1, NWAVE_RRTMG
      if (is_wet_l .and. (vol_shell_wet_l > 0._f)) then
        m_shell_l = (vol_water_added_l / vol_shell_wet_l) * refidx_water(iwave_l) + &
                    (vol_sulf_dry_l    / vol_shell_wet_l) * refidx_sulf(iwave_l)
      else
        m_shell_l = refidx_sulf(iwave_l)
      end if

      if (is_coreshell_l) then
        out_refidx_r(iwave_l,irun_l) = real(m_shell_l)
        out_refidx_i(iwave_l,irun_l) = aimag(m_shell_l)
      else
        m_eff_l = maxwell_garnett_sym(refidx_alum(iwave_l), m_shell_l, f_incl_l)
        out_refidx_r(iwave_l,irun_l) = real(m_eff_l)
        out_refidx_i(iwave_l,irun_l) = aimag(m_eff_l)
      end if
    end do

    do ibin_l = 1, NBIN
      radius_l    = r(ibin_l,I_GRP_ALSUL) * (rwet_bin_l / r(1,I_GRP_ALSUL))
      rcore_bin_l = r(ibin_l,I_GRP_ALSUL) * core_frac_l**(1._f/3._f)

      do iwave_l = 1, NWAVE_RRTMG
        if (is_wet_l .and. (vol_shell_wet_l > 0._f)) then
          m_shell_l = (vol_water_added_l / vol_shell_wet_l) * refidx_water(iwave_l) + &
                      (vol_sulf_dry_l    / vol_shell_wet_l) * refidx_sulf(iwave_l)
        else
          m_shell_l = refidx_sulf(iwave_l)
        end if

        if (is_coreshell_l) then
          call mie(carma, I_MIERTN_TOON1981, radius_l, wave_rrtmg(iwave_l), &
                   0._f, 3.0_f, 0.0_f, 1.0_f, m_shell_l, rcore_bin_l, refidx_alum(iwave_l), &
                   Qext_l, Qsca_l, asym_l, rc_l)
        else
          m_eff_l = maxwell_garnett_sym(refidx_alum(iwave_l), m_shell_l, f_incl_l)
          call mie(carma, I_MIERTN_TOON1981, radius_l, wave_rrtmg(iwave_l), &
                   0._f, 3.0_f, 0.0_f, 1.0_f, m_eff_l, 0.0_f, m_eff_l, &
                   Qext_l, Qsca_l, asym_l, rc_l)
        end if
        if (rc_l < 0) stop "    *** mie (ALSUL) FAILED ***"

        qext(iwave_l,ibin_l,irun_l) = Qext_l
        ssa (iwave_l,ibin_l,irun_l) = Qsca_l / Qext_l
        asy (iwave_l,ibin_l,irun_l) = asym_l
      end do

      out_r(ibin_l,irun_l)   = radius_l
      out_vol(ibin_l,irun_l) = (4._f/3._f) * PI * radius_l**3
      if (is_wet_l) then
        out_rmass(ibin_l,irun_l) = out_vol(ibin_l,irun_l) * rhopwet_l
      else
        out_rmass(ibin_l,irun_l) = out_vol(ibin_l,irun_l) * rho_dry_mix_l
      end if
    end do

  end subroutine fill_alsul_run

end subroutine test_alsul_coreshell_optics
