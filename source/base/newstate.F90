! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine manages the calculations that update state variables
!! of the model with new values at the current simulation time.
!!
!! @author McKie
!! @version Oct-1995
subroutine newstate(carma, cstate, rc)

  ! types
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod

  implicit none

  type(carma_type), intent(inout)      :: carma   !! the carma object
  type(carmastate_type), intent(inout) :: cstate  !! the carma state object
  integer, intent(inout)               :: rc      !! return code, negative indicates failure
  
  real(kind=f)                    :: sedlayer(NBIN,NELEM)
  real(kind=f)                    :: pcd_last(NBIN,NELEM)
  real(kind=f)                    :: pc_orig(NZ,NBIN,NELEM)
  real(kind=f)                    :: gc_orig(NZ,NGAS)
  real(kind=f)                    :: t_orig(NZ)
  real(kind=f)                    :: scale_cldfrc(NZ)
  integer                         :: kb
  integer                         :: ke
  integer                         :: idk
  integer                         :: iz
  integer                         :: isubstep
  integer                         :: igroup
  integer                         :: igas
  integer                         :: ielem
  integer                         :: ibin
  integer                         :: ntsubsteps
  logical                         :: takeSteps
  real(kind=f)                    :: fraction            ! Fraction of dT, dgc and pdc to be added in a substep.

  1 format(/,'newstate::ERROR - Substep failed, maximum retries execeed. : iz=',i4,',isubstep=',i12, &
             ',ntsubsteps=',i12,',nretries=',F9.0)

  ! Calculate changes due to vertical transport
  if (do_vtran) then
  
    call vertical(carma, cstate, rc)
    if (rc < RC_OK) return
  endif
  
  ! If doing doing incloud processing, then scale the parameters for incloud concentrations.
  ! 
  ! NOTE: Don't want to do this before sedimentation, since sedimentation doesn't take into
  ! account the varying cloud fractions, and thus a particle scaled at one level and cloud
  ! fraction would be scaled inappropriately at another level and cloud fraction.
  if (do_incloud) then
  
    ! Convert "cloud" particles to in-cloud values.
    !
    ! NOTE: If a particle is a "cloud" particle, it means that the entire mass of the
    ! particle is in the incloud portion of the grid box. Particle that are not "cloud
    ! particles" have their mass spread throughout the grid box.
    pc_orig(:,:,:) = pc(:,:,:)
    gc_orig(:,:)   = gc(:,:)
    t_orig(:)      = t(:)
    
    ! If the cloud fraction gets too small it causes the microphysics to require a
    ! lot of substeps. Enofrce a minimum cloud fraction for the purposes of scaling
    ! to incloud values.
    scale_cldfrc(:) = max(CLDFRC_MIN, cldfrc(:))
    
    do ielem = 1, NELEM
      igroup = igelem(ielem)
      
      if (is_grp_cloud(igroup)) then
        do ibin = 1, NBIN
          pc(:, ibin, ielem)  = pc(:, ibin, ielem)  / scale_cldfrc(:)
          pcd(:, ibin, ielem) = pcd(:, ibin, ielem) / scale_cldfrc(:)
        end do
      end if
    end do
  end if
  
  ! Redetermine the maximum particle values.
  if ((do_vtran) .or. do_incloud) then
    do iz = 1, NZ
      call maxconc(carma, cstate, iz, rc)
      if (rc < RC_OK) return
    end do
  end if
  
  ! Calculate changes in particle concentrations due to microphysical
  ! processes, part 1.  (potentially slower microphysical calcs)
  ! All spatial points are handled by one call to this routine.
  if (do_coag) then
    call microslow(carma, cstate, rc)
    if (rc < RC_OK) return
  endif
  
  ! If there is any microsphysics that happens on a faster time scale,
  ! then check to see if the time step needs to be subdivided and then
  ! perform the fast microphysical calculations.
  if (do_grow) then
  
    ! Set vertical loop index to increment downwards
    ! (for substepping of sedimentation)
    if (igridv .eq. I_CART) then
      kb  = NZ
      ke  = 1
      idk = -1
    else
      kb  = 1
      ke  = NZ
      idk = 1
    endif
    
    ! Initialize sedimentation source to zero at top of model
    dpc_sed(:,:) = 0._f

    ! Save the results from the slow operations, since we might need to retry the
    ! fast operations
    pcl(:,:,:) = pc(:,:,:)
    
    if (do_substep) then
      do igas = 1,NGAS
        gcl(:,igas) = gc(:,igas)
      end do
      told(:) = t(:)
    endif
    
      
    do iz = kb,ke,idk

      ! Compute or specify number of sub-timestep intervals for current spatial point
      ! (Could be same for all spatial pts, or could vary as a function of location)
      ntsubsteps = minsubsteps
      
      call nsubsteps(carma, cstate, iz, dtime_orig, ntsubsteps, rc)
      if (rc <  RC_OK) return
      
      ! Grab sedimentation source for entire step for this layer
      ! and set accumlated source for underlying layer to zero
      sedlayer(:,:) = dpc_sed(:,:)
      
      ! Do sub-timestepping for current spatial grid point, and allow for
      ! retrying should this level of substepping not be enough to keep the
      ! gas concentration from going negative.
      nretries = 0._f
      takeSteps = .true.
      
      do while (takeSteps)
      
        ! Compute sub-timestep time interval for current spatial grid point
        dtime = dtime_orig / ntsubsteps
        
        ! Don't retry unless requested.
        takeSteps = .false.

        ! Reset the amount that has been collected to sedimented down to the
        ! layer below.
        dpc_sed(:,:) = 0._f
        
        ! Reset the total nucleation for the step.
        pc_nucl(iz,:,:) = 0._f

        ! Remember the amount of detrained particles.
        if (do_detrain) then
          pcd_last(:,:) = pcd(iz,:,:)
        end if
        
        ! Reset average heating rates.
        rlheat(iz)     = 0._f
        partheat(iz)   = 0._f
        
        do isubstep = 1,ntsubsteps
          
          ! If substepping, then increment the gas concentration and the temperature by
          ! an amount for one substep.
          if (do_substep) then

            ! Since we don't really know how the gas and temperature changes arrived during the
            ! step, we can try different assumptions for how the gas and temperature are add to
            ! the values from the previous substep.
  
            ! Linear increment for substepping.
            fraction     = 1._f / ntsubsteps

            do igas = 1,NGAS
              gc(iz,igas) = gc(iz,igas) + d_gc(iz,igas) * fraction
            enddo
            
            t(iz) = t(iz) + d_t(iz) * fraction
          
 
            ! Detrainment puts the full gridbox amount into the incloud portion.
            if (do_detrain) then
              pc(iz,:,:)  = pc(iz,:,:)  + pcd_last(:,:) * fraction
              pcd(iz,:,:) = pcd(iz,:,:) - pcd_last(:,:) * fraction
            end if
          endif
            
  
          ! Redetermine maximum particle concentrations.
          call maxconc(carma, cstate, iz, rc)
          if (rc < RC_OK) return
          
          ! Calculate changes in particle concentrations for current spatial point
          ! due to microphysical processes, part 2.  (faster microphysical calcs)
          call microfast(carma, cstate, iz, rc)
          if (rc < RC_OK) return

  
          ! If there was a retry warning message and substepping is enabled, then retry
          ! the operation with more substepping.
          if (rc == RC_WARNING_RETRY) then
            if (do_substep) then
          
              ! Only retry for so long ...
              nretries = nretries + 1
              
              if (nretries > maxretries) then
                if (do_print) write(LUNOPRT,1) iz, isubstep, ntsubsteps, nretries - 1._f
                rc = RC_ERROR
                exit
              end if
            
              ! Try twice the substeps
              !
              ! NOTE: We are going to rely upon retries, so don't clutter the log
              ! with retry print statements. They slow down the run.
              ntsubsteps = ntsubsteps * 2
              
!              if (do_print) write(LUNOPRT,*) "newstate::WARNING - Substep failed, retrying with ", ntsubsteps, " substeps."
  
              ! Reset the state to the beginning of the step
              pc(iz,:,:) = pcl(iz,:,:)
              pcd(iz,:,:) = pcd_last(:,:)
              t(iz) = told(iz)
              do igas = 1,NGAS
                gc(iz,igas) = gcl(iz,igas)

                ! Now that we have reset the gas concentration, we need to recalculate the supersaturation.  
                call supersat(carma, cstate, iz, igas, rc)
                if (rc < RC_OK) return
              end do
              
              rc = RC_OK
              takeSteps = .true.
              exit
              
              
            ! If substepping is not enabled, than the retry warning should be treated as an error.
            else
            
              if (do_print) write(LUNOPRT,*) "newstate::ERROR - Step failed, suggest enabling substepping."
              rc = RC_ERROR
              exit
            end if            
          end if
        end do
      end do

      ! Keep track of substepping and retry statistics for performance tuning.
      max_nsubstep = max(max_nsubstep, ntsubsteps)
      max_nretry   = max(max_nretry, nretries)

      nstep    = nstep    + 1._f
      nsubstep = nsubstep + ntsubsteps
      nretry   = nretry   + nretries
      
      if (do_substep) zsubsteps(iz) = ntsubsteps
    end do
  
    ! Restore normal timestep
    dtime = dtime_orig
    
    ! Calculate average heating rates.
    rlheat(:)    = rlheat(:)   / dtime
    partheat(:)  = partheat(:) / dtime
    
  else
  
    ! If there is no reason to substep, but substepping was enabled, get the gas and
    ! temperature back to their final states.
    if (do_substep) then
  
      do igas = 1,NGAS
        gc(:,igas) = gc(:,igas) + d_gc(:,igas)
      enddo

      t(:) = t(:) + d_t(:)
    end if
    
    ! Do the detrainment, if it was being done in the growth loop.
    if (do_detrain) then
      pc(:,:,:)    = pc(:,:,:) + pcd(:,:,:)
      
      ! Remove the ice from the detrained ice, so that total ice will be conserved.
      pcd(:,:,:)   = 0._f
    end if
  end if
  
  
  ! Convert particles, gas and temperature to gridbox average values
  !
  ! NOTE: For particles that are not in the cloud, the unchanged value outside of the
  ! cloud needs to be merged with the changes in the part of the grid box that is
  ! occupied by cloud. The values used for the rest of the gridbox need to be the
  ! values from the original state.
  !
  ! For particles at the surface, assume a maximum cloud overlap.
  if (do_incloud) then
    do ielem = 1, NELEM
      igroup = igelem(ielem)
      
      if (is_grp_cloud(igroup)) then
        do ibin = 1, NBIN
          pc(:, ibin, ielem)   = pc(:, ibin, ielem) * scale_cldfrc(:)
        end do
      else
        do ibin = 1, NBIN
          pc(:, ibin, ielem)   = (1._f - scale_cldfrc(:)) * pc_orig(:, ibin, ielem) + pc(:, ibin, ielem) * scale_cldfrc(:)
        end do
      end if
    end do
        
    t(:) = (1._f - scale_cldfrc(:)) * t_orig(:) + scale_cldfrc(:) * t(:)
    
    if (do_substep) then
      t(:) = t(:) + (1._f - scale_cldfrc(:)) * d_t(:)
    end if
    
    do igas = 1, NGAS
      gc(:, igas) = (1._f - scale_cldfrc(:)) * gc_orig(:, igas) + gc(:, igas) * scale_cldfrc(:)
    
      if (do_substep) then
        gc(:, igas) = gc(:, igas) + (1._f - scale_cldfrc(:)) * d_gc(:, igas)
      end if

      ! Recalculate gridbox average supersaturation.
      do iz = 1, NZ
        call supersat(carma, cstate, iz, igas, rc)
        if (rc < RC_OK) return
      end do
    end do
  end if

    
  ! Return to caller with new state computed 
  return
end
