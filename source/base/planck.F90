! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

module planck

contains
  !! This routine calculates the planck intensity.
  !!
  !! This algorithm is based upon eqn 1.2.4 from Liou[2002].
  !!
  !! @author Chuck Bardeen
  !! @version Jan-2010
  function planckIntensity(wvl, temp)
  
    ! types
    use carma_precision_mod
    use carma_enums_mod
    use carma_constants_mod
    use carma_types_mod
    use carmastate_mod
    use carma_mod
  
    implicit none
  
    real(kind=f), intent(in)             :: wvl     !! wavelength (cm)
    real(kind=f), intent(in)             :: temp    !! temperature (K)
    real(kind=f)                         :: planckIntensity  !! Planck intensity (erg/s/cm2/sr/cm)
  
    ! Local declarations
    
    real(kind=f), parameter              :: C = 2.9979e10_f     ! Speed of light [cm/s]       
    real(kind=f), parameter              :: H = 6.62608e-27_f   ! Planck constant [erg s]
    
    ! Calculate the planck intensity.
    planckIntensity = 2._f * H * C**2 / ((wvl**5) * (exp(H * C / (BK * wvl * temp)) - 1._f))
  
    ! Return the planck intensity to the caller.
    return
  end function
end
