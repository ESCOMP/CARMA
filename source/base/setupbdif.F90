! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine evaluates particle vertical diffusion coefficients,
!! dkz(k,i,j) [cm^2 s^-1].
!!
!! Method: Uses equation 8.73 from Seinfeld and Pandis [1998] along
!! with the slip correction factor (bpm) calculated in the fall
!! velocity setup.
!!
!! This routine requires that vertical profiles of temperature <t>,
!! air density <rhoa>, viscosity <rmu>, and slip correction <bpm> are
!! defined (i.e., initatm.f and setupvf.f must be called before this).
!!
!! NOTE: Eddy diffusion is carried out by the parent model, so the only
!! diffusion that CARMA does is Brownian diffusion.
!!
!! @author Chuck Bardeen
!! @version Aug-2010
subroutine setupbdif(carma, cstate, rc)

  !     types
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod

  implicit none

  type(carma_type), intent(in)         :: carma   !! the carma object
  type(carmastate_type), intent(inout) :: cstate  !! the carma state object
  integer, intent(inout)               :: rc      !! return code, negative indicates failure

  ! Local declarations  
  integer                 :: igroup, i, j, k, k1, k2, ibin, iz, nzm1

  !  Define formats
  2 format(/,'Brownian diffusion coefficient (prior to interpolation)')
  3 format(/,'Particle group ',i3,/,' bin   lev  p [dyne/cm2]  T [K]       r [cm]    wet r [cm]    dkz [cm2/s]',/)
  4 format(i3,4x,i3,5(1pe11.3,4x)) 

  
  !  Loop over all groups.
  do igroup = 1, NGROUP
      
    ! Loop over particle size bins.
    do i = 1,NBIN

      ! Loop over all atltitudes.
      do k = 1, NZ

        ! Vertical brownian diffusion coefficient
        dkz(k,i,j) = (BK*t(k)*bpm(k,i,j)) / (6._f*PI*rmu(k)*r_wet(k,i,j))

      enddo
    enddo
  enddo
    
  ! Print out diffusivities.
!#ifdef DEBUG
  if (do_print_init) then
   
    write(LUNOPRT,2)

    do j = 1, NGROUP
        
      write(LUNOPRT,3) j
      
      do i = 1,NBIN
    
        do k = NZ, 1, -1
          write(LUNOPRT,4) i,k,p(k),t(k),r(i,j),r_wet(k,i,j),dkz(k,i,j)
        end do
      enddo
    enddo
  
    write(LUNOPRT,*) ""
  end if
!#endif

  !  Interpolate <dkz> from layer mid-pts to layer boundaries.
  !  <dkz(k)> is the diffusion coefficient at the lower edge of the layer
  nzm1 = max(1, NZ-1)

  ! Set upper boundary before averaging
  dkz(NZP1,:,:) = dkz(NZ,:,:)

  if (NZ .gt. 1) then
    dkz(NZ,:,:) = sqrt(dkz(nzm1,:,:) * dkz(NZ,:,:))

    if (NZ .gt. 2) then
      do iz = NZ-1, 2, -1
        dkz(iz,:,:) = sqrt(dkz(iz-1,:,:) * dkz(iz,:,:))
      enddo
    endif
  endif
  
  ! Scale cartesian diffusivities to the appropriate vertical coordinate system.
  ! Non--cartesion coordinates are assumed to be positive downward, but
  ! vertical velocities in this model are always assumed to be positive upward. 
  if( igridv .ne. I_CART )then
    do igroup=1,NGROUP
      do ibin=1,NBIN
        dkz(:,ibin,igroup) = dkz(:,ibin,igroup) / (zmetl(:)**2)
      enddo
    enddo
  endif

  ! Return to caller with fall velocities evaluated.
  return
end
