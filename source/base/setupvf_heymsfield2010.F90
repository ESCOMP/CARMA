! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine evaluates particle fall velocities, vf(k) [cm s^-1]
!! and reynolds' numbers based on fall velocities, re(j,i,k) [dimensionless].
!! indices correspond to vertical level <k>, bin index <i>, and aerosol
!! group <j>.
!!
!! Method: Use the routined from Heymsfield and Westbrook [2010], which is
!! designed only for ice particles. Thus this routine uses the dry mass and
!! radius, not the wet mass and radius. The area ration (Ar) is determined
!! based upon the formulation of Schmitt and Heymsfield [JAS, 2009].
!!
!! @author  Chuck Bardeen
!! @version Mar-2010 
subroutine setupvf_heymsfield2010(carma, cstate, j, rc)

  ! types
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod

  implicit none

  type(carma_type), intent(in)         :: carma   !! the carma object
  type(carmastate_type), intent(inout) :: cstate  !! the carma state object
  integer, intent(in)                  :: j       !! group index
  integer, intent(inout)               :: rc      !! return code, negative indicates failure

  ! Local declarations
  integer                  :: i, k
  real(kind=f)             :: rhoa_cgs, vg, rmfp, rkn, expon, x
  real(kind=f), parameter  :: c0      = 0.35_f
  real(kind=f), parameter  :: delta0  = 8.0_f
  
  real(kind=f)             :: ar                  ! area ratio
  real(kind=f)             :: dmax                ! maximum diameter
  
  
  ! The area ratio is the ratio of the area of the shape to the area of the
  ! circumscribing circle.
  !
  ! Default to a sphere.
  if (ishape(j) .eq. I_HEXAGON) then
    ar = 0.8270_f
  else if (ishape(j) .eq. I_CYLINDER) then
    ar = 1.0_f
  else
    ar = 1.0_f
  end if
  
  ! Loop over all atltitudes.
  do k = 1, NZ

    ! This is <rhoa> in cartesian coordinates (good old cgs units)
    rhoa_cgs = rhoa(k) / (xmet(k)*ymet(k)*zmet(k))

    ! <vg> is mean thermal velocity of air molecules [cm/s]
    vg = sqrt(8._f / PI * R_AIR * t(k))

    ! <rmfp> is mean free path of air molecules [cm]
    rmfp = 2._f * rmu(k) / (rhoa_cgs * vg)

    ! Loop over particle size bins.
    do i = 1,NBIN
    
      ! <rkn> is knudsen number
      rkn = rmfp / r_wet(k,i,j)

      ! <bpm> is the slip correction factor, the correction term for
      ! non-continuum effects.  Also used to calculate coagulation kernels
      ! and diffusion coefficients.
      expon = -.87_f / rkn
      expon = max(-POWMAX, expon)
      bpm(k,i,j) = 1._f + (1.246_f*rkn + 0.42_f*rkn*exp(expon))

      dmax = 2._f * r(i,j)
      
      if (ishape(j) .eq. I_HEXAGON) then
      
        ! Get the minimum diameter of the hexagon.
        dmax = dmax * 0.8456_f * eshape(j)**(-ONE/3._f)
        
        ! And now get the maximum diameter.
        dmax = 2._f / sqrt(3._f) * dmax
      else if (ishape(j) .eq. I_CYLINDER) then
        dmax = dmax * 0.8736_f * eshape(j)**(-ONE/3._f)
      end if

      ! Determine the area ratio based on the formulation given in Schmitt and Heymsfield
      ! [2009].
!      if (dmax <= 200.e-4_f) then
!        ar = exp(-38._f * dmax)
!      else
!        ar = 0.16_f * (dmax ** (-0.27_f))
!      end if

  
      x = (rhoa_cgs / (rmu(k)**2)) * ((8._f * rmass(i,j) * GRAV) / (PI * (ar**0.5_f)))
      
      ! Apply the slip correction factor. This is not included in the formulation
      ! from Heymsfield and Westbrook [2010].
      !
      ! NOTE: This is applied according to eq 8.46 and surrounding discussion in
      ! Seinfeld and Pandis [1998].
      x = x * bpm(k,i,j)
      
      re(k,i,j) = ((delta0**2) / 4._f) * (sqrt(1._f + (4._f * sqrt(x) / (delta0**2 * sqrt(c0)))) - 1._f)**2
      
      
      vf(k,i,j) = rmu(k) * re(k,i,j) / (rhoa_cgs * dmax)
    enddo      ! <i=1,NBIN>
  enddo      ! <k=1,NZ>
  
  ! Return to caller with particle fall velocities evaluated.
  return
end
