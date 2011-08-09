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
  

  !! This routine calculates the total planck intensity from the specified
  !! wavelength to a wavelength of 0.
  !!
  !! This algorithm is based upon Widger and Woodall, BAMS, 1976 as
  !! indicated at http://www.spectralcalc.com/blackbody/appendixA.html.
  !!
  !! @author Chuck Bardeen
  !! @version Aug-2011
  function planckIntensityWidger1976(wvl, temp, miniter)
  
    ! types
    use carma_precision_mod
    use carma_enums_mod
    use carma_constants_mod
    use carma_types_mod
    use carmastate_mod
    use carma_mod
  
    implicit none
 
    real(kind=f), intent(in)             :: wvl     !! band center wavelength (cm)
    real(kind=f), intent(in)             :: temp    !! temperature (K)
    integer, intent(in)                  :: miniter !! minimum iterations
    real(kind=f)                         :: planckIntensityWidger1976  !! Planck intensity (erg/s/cm2/sr/cm)

    ! Local Variables
    real(kind=f), parameter              :: C  = 299792458.0_f     ! Speed of light [m/s]       
    real(kind=f), parameter              :: H  = 6.6260693e-34_f   ! Planck constant [J s]
    real(kind=f), parameter              :: BZ = 1.380658e-23_f    ! Boltzman constant

    real(kind=f)                         :: c1, x, x2, x3, sumJ, dn, sigma
    integer                              :: iter, n

    sigma = 1._f / wvl
    
    c1 = H * C / BZ
    x  = c1 * 100._f * sigma / temp
    x2 = x * x
    x3 = x2 * x

    ! Use fewer iterations, since speed is more important than accuracy for
    ! the particle heating code, and even with fewer iterations the results
    ! with CAM bands still show good accuracy.
!    iter = min(512, int(2._f + 20._f / x))
    iter = min(miniter, int(2._f + 20._f / x))

    sumJ = 0._f

    do n = 1, iter
      dn  = 1._f / n
      sumJ = sumJ + exp(-n*x) * (x3 + (3.0_f * x2 + 6.0_f * (x + dn) * dn) * dn) * dn
    end do

    ! Convert results from W/m2/sr to erg/cm2/s/sr 
    planckIntensityWidger1976 = 2.0_f * H * (C**2) * ((temp / c1) ** 4) * sumJ * 1e7_f / 1e4_f

    return
  end function

 
  !! This routine calculates the average planck intensity in the wavelength
  !! band defined by wvl and dwvl.
  !!
  !! This algorithm is based upon Widger and Woodall, BAMS, 1976 as
  !! indicated at http://www.spectralcalc.com/blackbody/appendixA.html.
  !!
  !! @author Chuck Bardeen
  !! @version Aug-2011
  function planckBandIntensityWidger1976(wvl, dwvl, temp, miniter)
  
    ! types
    use carma_precision_mod
    use carma_enums_mod
    use carma_constants_mod
    use carma_types_mod
    use carmastate_mod
    use carma_mod
  
    implicit none
  
    real(kind=f), intent(in)             :: wvl     !! band center wavelength (cm)
    real(kind=f), intent(in)             :: dwvl    !! band width (cm)
    real(kind=f), intent(in)             :: temp    !! temperature (K)
    integer, intent(in)                  :: miniter !! minimum iterations
    real(kind=f)                         :: planckBandIntensityWidger1976  !! Planck intensity (erg/s/cm2/sr/cm)
    
    ! Calculate the integral from the edges to 0 and subtract.
    planckBandIntensityWidger1976 = (planckIntensityWidger1976(wvl + (dwvl / 2._f), temp, miniter) - planckIntensityWidger1976(wvl - (dwvl / 2._f), temp, miniter)) / dwvl

    return 
  end function


  !! This routine calculates the average planck intensity in the wavelength
  !! band defined by wvl and dwvl.
  !!
  !! This algorithm does a brute force integral by dividing the band into
  !! small sub-bands. This routine can be slow.
  !!
  !! @author Chuck Bardeen
  !! @version Aug-2011
  function planckBandIntensity(wvl, dwvl, temp, iter)
  
    ! types
    use carma_precision_mod
    use carma_enums_mod
    use carma_constants_mod
    use carma_types_mod
    use carmastate_mod
    use carma_mod

    implicit none
  
    real(kind=f), intent(in)             :: wvl     !! band center wavelength (cm)
    real(kind=f), intent(in)             :: dwvl    !! band width (cm)
    real(kind=f), intent(in)             :: temp    !! temperature (K)
    integer, intent(in)                  :: iter    !! number of iterations
    real(kind=f)                         :: planckBandIntensity  !! Planck intensity (erg/s/cm2/sr/cm)
    
    ! Local Variables
    real(kind=f)                         :: wstart    ! Starting wavelength (cm)
    real(kind=f)                         :: ddwave    ! sub-band width (cm)
    integer                              :: i

    wstart = wvl - (dwvl / 2._f)
    ddwave = dwvl / iter

    planckBandIntensity = 0._f

    do i = 1, iter
      planckBandIntensity = planckBandIntensity + planckIntensity(wstart + (i - 0.5) * ddwave, temp) * ddwave
    end do
 
    planckBandIntensity = planckBandIntensity / dwvl

    return
  end function
  
  
  !! This routine calculates the average planck intensity in the wavelength
  !! band defined by wvl and dwvl.
  !!
  !! error computed on full spectrum compared to planck function.  Band-levels may be different
  !! 8.9%   error with   5 quadrature points in [100 micrometer, 1 millimeter]
  !! 1.7%   error with  10 quadrature points in [100 micrometer, 1 millimeter]
  !! 0.001% error with 100 quadrature points in [100 micrometer, 1 millimeter]
  !!
  !! For most RRTMG bands, 3 quadrature points are probably sufficient, but testing is
  !! left to the reader.
  !!
  !! @author Andrew Conley, Chuck Bardeen
  !! @version Aug-2011
  function planckBandIntensityConley2011(wvl, dwvl, temp, iter)

    ! types
    use carma_precision_mod
    use carma_enums_mod
    use carma_constants_mod
    use carma_types_mod
    use carmastate_mod
    use carma_mod

    implicit none
  
    real(kind=f), intent(in)             :: wvl     !! band center wavelength (cm)
    real(kind=f), intent(in)             :: dwvl    !! band width (cm)
    real(kind=f), intent(in)             :: temp    !! temperature (K)
    integer, intent(in)                  :: iter    !! number of iterations
    real(kind=f)                         :: planckBandIntensityConley2011  !! Planck intensity (erg/s/cm2/sr/cm)
    
    real(kind=f) :: k = 1.3806488e-23_f     ! boltzmann J/K
    real(kind=f) :: c = 2.99792458e8_f      ! light m/s
    real(kind=f) :: h = 6.62606957e-34_f    ! planck J s
    real(kind=f) :: lambda1                 ! low wavelength m
    real(kind=f) :: lambda2                 ! high wavelength m
    real(kind=f) :: sigma = 5.670373e-8_f   ! stef-bolt W/m/m/k/k/k/k

    ! quadrature iteration
    integer :: i,inumber 
    real(kind=f) :: fnumber
    
    
    real(kind=f) :: fr1, fr2         ! frequency bounds of partition
    real(kind=f) :: kt               ! k_boltzmann * temperature
    real(kind=f) :: l1,l2, ll1, ll2  ! lower and upper bounds of (wavelength) log(wavelength)
    real(kind=f) :: t1,t3            ! 1st and 3rd order terms
    real(kind=f) :: total, total2    ! 1st and 3rd order cumulative partial integral
    real(kind=f) :: dfr,m,e,d,a,o,tt ! terms appearing in integral
    real(kind=f) :: coeff            ! front coefficient of integral
    real(kind=f) :: planck           ! planck function
    
    total   = 0._f
    total2  = 0._f
    inumber = iter
    fnumber = 1._f * inumber
    
    kt  = k*temp
    lambda1 = (wvl - (dwvl / 2._f)) * 1e-2
    lambda2 = (wvl + (dwvl / 2._f)) * 1e-2
    ll1 = log(lambda1)
    ll2 = log(lambda2)
  
!    print *,'number of quadrature points', inumber
!    print *,'wavelength range',lambda1,lambda2
  
    ! step through partition
    do i = 1,inumber
    
     ! find bounds, exponentially distributed on freqency 
     l1 = lambda1 * exp( (i-1)*(ll2-ll1)/fnumber) !lambda1 + (i-1)*(lambda2-lambda1)/50._f
     l2 = lambda1 * exp( i*(ll2-ll1)/fnumber) !l1 + (lambda2-lambda1)/50._f
     fr1 = c/l2
     fr2 = c/l1
     !print*, 'wavelength interval',l1,l2
    
     ! integrate
    
     ! constants
     dfr = (fr2-fr1)/2._f  ! delta freq
     m = (fr1+fr2)/2._f    ! mean freq
     e = exp(h*m/kt)        ! exponential
     d = 1._f-1._f/e      ! delta
     a = h/kt               ! alpha
    
     ! integrals
     o = fr2-fr1                     ! int 1 deps
     tt = 2._f*(dfr*dfr*dfr)/3._f   ! int eps^2 deps
     !print*,'dfr',dfr,'m',m,'e',e,'d',d,'a',a,'o',o,'t',t
    
     ! frontpiece
     coeff = 2._f*h*m*m*m/c/c/(e-1._f)
     ! term and 3'rd order correction
     t1 = o
     t3 = 3._f/m/m - 3._f*a/d/m + a*a/d/d - a*a/2._f/d
     
     ! sum it up.  Total is 3rd order, total2 is 1st order
     total  = total + coeff*(t1+tt*t3)
!       total2 = total2 + coeff*t1
  
     end do
     
     ! Convert to erg/cm2/s/sr/cm
     planckBandIntensityConley2011 = total * 1e7 / 1e4 / dwvl
     
     return
  end function
end
