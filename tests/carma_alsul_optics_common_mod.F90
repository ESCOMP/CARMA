!! Shared constants and helper routines for the ALSUL (mixed alumina/sulfate)
!! core/shell optics comparison tests:
!!   carma_alsul_coreshell_optics_test.F90  - compares composition-mixing methods
!!   carma_alsul_fractal_optics_test.F90    - compares sphere vs. fractal geometry
!!
!! The wavelength grid and the sulfate/water refractive index tables below are
!! copied from the CAM SourceMods carma_model_mod.F90 (fixed_sst_paper case),
!! so that results here are directly comparable to what CAM/RRTMG would see.
!! The wave_rrtmg grid is derived from RRTMG's actual hardcoded band-boundary
!! wavenumber tables (source/physics/rrtmg/radconstants.F90: wavenumber1/2_longwave,
!! wavenum_low/high), matching exactly what CAM's get_lw/sw_spectral_boundaries
!! would return in 'cm' units (band low edge is what CAM stores into its "wave"
!! array, despite the "band center" naming in the SourceMods code).
module carma_alsul_optics_common_mod

  use carma_precision_mod
  use carma_constants_mod, only: PI, CARMA_NAME_LEN, CARMA_SHORT_NAME_LEN

  implicit none

  private
  public :: NWAVE_RRTMG, NLWBANDS_RRTMG, NSWBANDS_RRTMG
  public :: wave_rrtmg
  public :: shellreal, shellimag, waterreal, waterimag
  public :: get_alumina_refidx
  public :: maxwell_garnett, maxwell_garnett_sym

  ! Geometry / group parameters shared by both tests (mirrors carma_alumsulf_2nc_test.F90).
  public :: rmin_PRSUL, rmrat_PRSUL, rmin_ALUM, rmrat_ALUM
  public :: RHO_SULFATE, RHO_ALUMINUM
  public :: Kappa_SULF, Kappa_zero
  public :: rmon_ALUM, falpha_ALUM, df_ALSUL_value

  integer, parameter :: NLWBANDS_RRTMG = 16
  integer, parameter :: NSWBANDS_RRTMG = 14
  integer, parameter :: NWAVE_RRTMG    = NLWBANDS_RRTMG + NSWBANDS_RRTMG   ! = 30

  ! Band low-edge wavelength (cm), in the same order/units CAM's carma_model_mod.F90
  ! stores into its "wave" array: LW bands 1-16 first, then SW bands 17-30.
  real(kind=f), parameter :: wave_rrtmg(NWAVE_RRTMG) = (/ &
       2.8571428571e-03_f, 2.0000000000e-03_f, 1.5873015873e-03_f, 1.4285714286e-03_f, &
       1.2195121951e-03_f, 1.0204081633e-03_f, 9.2592592593e-04_f, 8.4745762712e-04_f, &
       7.1942446043e-04_f, 6.7567567568e-04_f, 5.5555555556e-04_f, 4.8076923077e-04_f, &
       4.4444444444e-04_f, 4.1841004184e-04_f, 3.8461538462e-04_f, 3.0769230769e-04_f, &
       3.0769230769e-04_f, 2.5000000000e-04_f, 2.1505376344e-04_f, 1.9417475728e-04_f, &
       1.6260162602e-04_f, 1.2987012987e-04_f, 1.2422360248e-04_f, 7.7821011673e-05_f, &
       6.2500000000e-05_f, 4.4150110375e-05_f, 3.4482758621e-05_f, 2.6315789474e-05_f, &
       2.0000000000e-05_f, 3.8461538462e-04_f /)

  ! Sulfate (H2SO4/H2O shell) refractive index, from carma_model_mod.F90 (MAM-derived values).
  real(kind=f), parameter :: shellreal(NWAVE_RRTMG) = (/ 1.89_f, 1.912857_f, 1.932063_f, 1.586032_f, &
               1.677979_f, 1.757825_f, 1.855336_f, 1.596767_f, 1.146559_f, 1.261314_f, 1.424219_f, &
               1.351645_f, 1.378697_f, 1.385_f, 1.385_f, 1.366909_f, 1.366909_f, 1.314577_f, &
               1.357978_f, 1.380309_f, 1.392645_f, 1.404506_f, 1.412181_f, 1.421632_f, &
               1.427968_f, 1.430335_f, 1.441641_f, 1.467642_f, 1.484_f, 1.164128_f /)

  real(kind=f), parameter :: shellimag(NWAVE_RRTMG) = (/ 0.22_f, 0.15185711_f, 0.08457167_f, 0.22250789_f, 0.19499999_f, &
              0.44068847_f, 0.69594361_f, 0.69466153_f, 0.45876573_f, 0.16060575_f, &
              0.1715766_f , 0.14352135_f, 0.12025213_f, 0.12222873_f, 0.12581848_f, 0.15793008_f, &
              1.57930076e-01_f, 5.66869128e-02_f, 2.88634387e-03_f, 1.49071286e-03_f, &
              5.30385233e-04_f, 1.02977119e-04_f, 1.61967358e-05_f, 1.75122678e-06_f, &
              2.21435655e-08_f, 9.99999994e-09_f, 9.99999994e-09_f, 9.99999994e-09_f, &
              9.99999994e-09_f, 5.51133746e-01_f /)

  ! Water refractive index, from carma_model_mod.F90.
  real(kind=f), parameter :: waterreal(NWAVE_RRTMG) = &
   (/ 1.532_f, 1.523857_f, 1.420063_f, 1.274308_f, &
      1.161387_f, 1.142222_f, 1.232189_f, 1.266436_f, 1.295687_f, 1.320659_f, 1.341516_f, &
      1.315192_f, 1.330235_f, 1.339058_f, 1.350425_f, 1.408042_f, 1.408042_f, 1.324462_f, &
      1.276726_f, 1.301847_f, 1.312051_f, 1.321301_f, 1.322836_f, 1.326836_f, 1.330968_f, &
      1.33367_f, 1.339547_f, 1.348521_f, 1.362_f, 1.290783_f /)

  real(kind=f), parameter :: waterimag(NWAVE_RRTMG) = &
   (/ 0.336_f, 0.36000001_f, 0.42623809_f, 0.40341724_f, &
      0.32062717_f, 0.11484398_f, 0.04710282_f, 0.03901278_f, 0.03373134_f, 0.03437707_f, &
      0.09216518_f, 0.0121094_f, 0.01314786_f, 0.01013119_f, 0.00486624_f, 0.0142042_f, &
      1.42042044e-02_f, 1.57659209e-01_f, 1.51634401e-03_f, 1.15906247e-03_f, &
      2.35527521e-04_f, 1.71196912e-04_f, 2.43626002e-05_f, 3.12758360e-06_f, &
      3.74323598e-08_f, 1.63841034e-09_f, 2.49434956e-09_f, 1.52413800e-08_f, &
      3.35000010e-08_f, 3.43825518e-02_f /)

  ! Group/element geometry, matching carma_alumsulf_2nc_test.F90.
  real(kind=f), parameter :: rmin_PRSUL  = 3.43e-8_f    ! sulfate bin 1 radius [cm]
  real(kind=f), parameter :: rmrat_PRSUL = 3.67_f        ! sulfate volume ratio between bins
  real(kind=f), parameter :: rmin_ALUM   = 24.e-6_f      ! alumina bin 1 radius = monomer radius [cm]
  real(kind=f), parameter :: rmrat_ALUM  = 1.5_f         ! alumina volume ratio between bins
  real(kind=f), parameter :: rmon_ALUM   = 24.e-6_f      ! alumina monomer radius [cm]
  real(kind=f), parameter :: falpha_ALUM = 1._f          ! fractal packing coefficient
  real(kind=f), parameter :: df_ALSUL_value = 2.6_f      ! fractal dimension

  real(kind=f), parameter :: RHO_SULFATE  = 1.923_f      ! dry density of sulfate [g/cm3]
  real(kind=f), parameter :: RHO_ALUMINUM = 3.95_f       ! dry density of alumina [g/cm3]

  real(kind=f), parameter :: Kappa_SULF = 0.5_f
  real(kind=f), parameter :: Kappa_zero = 0.0_f

contains

  !! Reads alumina's lab-measured refractive index (real, imag vs. wavelength in
  !! cm) from a netcdf file with variables "wavelength", "real", "imag" - the
  !! John Dykema lab-measurement dataset referenced by the CAM SourceMods'
  !! al2o3lab_file - and nearest-neighbor matches it onto the
  !! wave_rrtmg grid, exactly as carma_model_mod.F90 does for refidx_labaer.
  subroutine get_alumina_refidx(filename, refidx_alum, rc)
    use netcdf

    character(len=*), intent(in)     :: filename
    complex(kind=f), intent(out)     :: refidx_alum(NWAVE_RRTMG)
    integer, intent(out)             :: rc

    integer :: inid, wave_varid, real_varid, imag_varid, wave_dimid
    integer :: lab_nwave, i, lab_idx
    real(kind=f), allocatable :: lab_wave(:), lab_real(:), lab_imag(:)

    rc = 0

    call ncheck(nf90_open(filename, nf90_nowrite, inid), rc)
    if (rc < 0) return

    call ncheck(nf90_inq_dimid(inid, "wavelength", wave_dimid), rc)
    if (rc < 0) return
    call ncheck(nf90_inquire_dimension(inid, wave_dimid, len=lab_nwave), rc)
    if (rc < 0) return

    allocate(lab_wave(lab_nwave), lab_real(lab_nwave), lab_imag(lab_nwave))

    call ncheck(nf90_inq_varid(inid, "wavelength", wave_varid), rc)
    if (rc == 0) call ncheck(nf90_get_var(inid, wave_varid, lab_wave), rc)

    call ncheck(nf90_inq_varid(inid, "real", real_varid), rc)
    if (rc == 0) call ncheck(nf90_get_var(inid, real_varid, lab_real), rc)

    call ncheck(nf90_inq_varid(inid, "imag", imag_varid), rc)
    if (rc == 0) call ncheck(nf90_get_var(inid, imag_varid, lab_imag), rc)

    call ncheck(nf90_close(inid), rc)

    if (rc == 0) then
      ! lab_wave and wave_rrtmg are both in cm - no unit conversion needed.
      do i = 1, NWAVE_RRTMG
        lab_idx = minloc(abs(lab_wave(:) - wave_rrtmg(i)), dim=1)
        ! Floor a near-/exactly-zero imaginary refractive index at a tiny
        ! epsilon - a perfectly lossless material at a non-trivial size
        ! parameter can trigger extremely sharp, narrow Mie resonances that
        ! send adaptive-quadrature Mie/fractal solvers into very slow (or
        ! effectively non-terminating) convergence. The production
        ! shellimag/waterimag tables in carma_alsul_optics_common_mod already
        ! carry this same ~1e-8 floor for identical reasons; the raw
        ! Alumina_dykema.nc lab file does not, so it's applied here instead.
        refidx_alum(i) = cmplx(lab_real(lab_idx), max(lab_imag(lab_idx), 1.0e-8_f), kind=f)
      end do
    end if

    deallocate(lab_wave, lab_real, lab_imag)

  end subroutine get_alumina_refidx


  !! Maxwell-Garnett effective-medium mixing rule: returns the effective
  !! refractive index of a host medium (m_host) containing a volume fraction
  !! f_incl of spherical inclusions with refractive index m_incl.
  complex(kind=f) function maxwell_garnett(m_incl, m_host, f_incl) result(m_eff)
    complex(kind=f), intent(in) :: m_incl
    complex(kind=f), intent(in) :: m_host
    real(kind=f), intent(in)    :: f_incl

    complex(kind=f) :: eps_incl, eps_host, eps_eff, b

    eps_incl = m_incl * m_incl
    eps_host = m_host * m_host

    b = f_incl * (eps_incl - eps_host) / (eps_incl + 2._f * eps_host)
    eps_eff = eps_host * (1._f + 2._f * b) / (1._f - b)

    m_eff = sqrt(eps_eff)

  end function maxwell_garnett


  !! Symmetrized Maxwell-Garnett mixing: maxwell_garnett() above is only valid
  !! for a DILUTE inclusion phase (f_incl well below 0.5) - for materials like
  !! alumina with strong far-IR absorption resonances, calling it with a
  !! majority-volume "inclusion" (e.g. f_a > 0.5) can produce an unphysical
  !! effective refractive index (very large real/imaginary parts) that then
  !! sends downstream Mie/fractal solvers into non-convergent behavior. This
  !! wrapper always treats whichever phase is the minority by volume as the
  !! inclusion, which keeps the mixing formula within its valid regime.
  complex(kind=f) function maxwell_garnett_sym(m_a, m_b, f_a) result(m_eff)
    complex(kind=f), intent(in) :: m_a   !! refractive index of phase A
    complex(kind=f), intent(in) :: m_b   !! refractive index of phase B
    real(kind=f), intent(in)    :: f_a   !! volume fraction of phase A (0-1)

    if (f_a <= 0.5_f) then
      m_eff = maxwell_garnett(m_a, m_b, f_a)
    else
      m_eff = maxwell_garnett(m_b, m_a, 1._f - f_a)
    end if

  end function maxwell_garnett_sym


  !! Local netcdf error checker (mirrors the "check" subroutine in
  !! carma_fractaloptics_2nc_test.F90, but returns rc instead of stopping,
  !! so callers can report the failing group/variant).
  subroutine ncheck(status, rc)
    use netcdf
    integer, intent(in)    :: status
    integer, intent(inout) :: rc

    if (status /= nf90_noerr) then
      print *, 'carma_alsul_optics_common_mod: ', trim(nf90_strerror(status))
      rc = -1
    end if
  end subroutine ncheck

end module carma_alsul_optics_common_mod
