! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!!  This routine evaluates supersaturations <supsatl> and <supsati> for all gases.
!!
!! @author Andy Ackerman, Chuck Bardeen
!! @version Dec-1995, Aug-2010
subroutine supersat(carma, cstate, iz, igas, rc)

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
  integer, intent(in)                  :: iz      !! z index
  integer, intent(in)                  :: igas    !! gas index
  integer, intent(inout)               :: rc      !! return code, negative indicates failure

  ! Local declarations
  real(kind=f)  :: rvap
  real(kind=f)  :: gc_cgs
  real(kind=f)  :: alpha

  ! Calculate vapor pressures.
  call vaporp(carma, cstate, iz, igas, rc)

  ! Define gas constant for this gas
  rvap = RGAS / gwtmol(igas)

  gc_cgs = gc(iz,igas) / (zmet(iz)*xmet(iz)*ymet(iz))

  supsatl(iz,igas) = (gc_cgs * rvap * t(iz) - pvapl(iz,igas)) / pvapl(iz,igas)
  supsati(iz,igas) = (gc_cgs * rvap * t(iz) - pvapi(iz,igas)) / pvapi(iz,igas)

  return
end
