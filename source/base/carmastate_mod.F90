!! The CARMA state module contains the atmospheric data for use with the CARMA
!! module. This implementation has been customized to work within other model 
!! frameworks. CARMA adds a lot of extra state information (atmospheric
!! properties, fall velocities, coagulation kernels, growth kernels, ...) and
!! thus has a large memory footprint. Because only one column will be operated
!! upon at a time per thread, only one cstate object needs to be instantiated
!! at a time and each cstate object only represents one column. This keeps
!! the memory requirements of CARMA to a minimum.
!!
!! @version Feb-2009 
!! @author  Chuck Bardeen, Pete Colarco, Jamie Smith 
!
! NOTE: Documentation for this code can be generated automatically using f90doc,
! which is freely available from:
!   http://erikdemaine.org/software/f90doc/
! Comment lines with double comment characters are processed by f90doc, and there are
! some special characters added to the comments to control the documentation process.
! In addition to the special characters mentioned in the f990doc documentation, html
! formatting tags (e.g. <i></i>, <sup></sup>, ...) can also be added to the f90doc
! comments.
module carmastate_mod

  ! This module maps the parents models constants into the constants need by CARMA.
  ! NOTE: CARMA constants are in CGS units, while the parent models are typically in
  ! MKS units.
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod

  ! cstate explicitly declares all variables. 
  implicit none

  ! All cstate variables and procedures are private except those explicitly
  ! declared to be public.
  private
  
  ! Declare the public methods.
  public CARMASTATE_Create
  public CARMASTATE_CreateFromReference
  public CARMASTATE_Destroy
  public CARMASTATE_Get
  public CARMASTATE_GetBin
  public CARMASTATE_GetDetrain
  public CARMASTATE_GetGas
  public CARMASTATE_GetState
  public CARMASTATE_SetBin
  public CARMASTATE_SetDetrain
  public CARMASTATE_SetGas
  public CARMASTATE_SetState
  public CARMASTATE_Step
  
contains
  
  ! These are the methods that provide the interface between the parent model and
  ! the atmospheric state data of the CARMA microphysical model. There are many other
  ! methods that are not in this file that are used to implement the microphysical
  ! calculations needed by the CARMA model. These other methods are in effect private
  ! methods of the CARMA module, but are in individual files since that is the way that
  ! CARMA has traditionally been structured and where users may want to extend or
  ! replace code to affect the microphysics.

  !! Create the CARMASTATE object, which contains information about the
  !! atmospheric state. Internally, CARMA uses CGS units, but this interface uses
  !! MKS units which are more commonly used in parent models. The units and grid
  !! orientation depend on the grid type:
  !!
  !!  - igridh
  !!    - I_CART   : Cartesian coordinates, units in [m]
  !!    - I_LL     : Lat/Lon coordinates, units in [degrees]
  !!
  !!  - igridv
  !!    - I_CART   : Cartesian coordinates, units in [m], bottom at NZ=1
  !!    - I_SIG    : Sigma coordinates, unitless [P/P0], top at NZ=1
  !!    - I_HYBRID : Hybrid coordinates, unitless [~P/P0], top at NZ=1
  !!
  !! NOTE: The supplied CARMA object should already have been created, configured,
  !! and initialized.
  !!
  !! NOTE: The relative humidity is optional, but needs to be supplied if particles
  !! are subject to swelling based upon relative humidity. The specific humdity can
  !! can be specified instead. If both are specified, then the realtive humidity is
  !! used.
  !!
  !! @author Chuck Bardeen
  !! @version Feb-2009
  !! @see CARMA_Create
  !! @see CARMA_Initialize
  !! @see CARMASTATE_Destroy
  subroutine CARMASTATE_Create(cstate, carma_ptr, time, dtime, NZ, igridv, igridh,  &
      lat, lon, xc, dx, yc, dy, zc, zl, p, pl, t, rc, qh2o, relhum, told)
    type(carmastate_type), intent(inout)    :: cstate      !! the carma state object
    type(carma_type), pointer, intent(in)   :: carma_ptr   !! (in) the carma object
    real(kind=f), intent(in)                :: time        !! the model time [s]
    real(kind=f), intent(in)                :: dtime       !! the timestep size [s]
    integer, intent(in)                     :: NZ          !! the number of vertical grid points
    integer, intent(in)                     :: igridv      !! vertical grid type
    integer, intent(in)                     :: igridh      !! horizontal grid type
    real(kind=f), intent(in)                :: lat         !! latitude at center [degrees north]
    real(kind=f), intent(in)                :: lon         !! longitude at center [degrees east]
    real(kind=f), intent(in)                :: xc(NZ)      !! x at center
    real(kind=f), intent(in)                :: dx(NZ)      !! ix width
    real(kind=f), intent(in)                :: yc(NZ)      !! y at center
    real(kind=f), intent(in)                :: dy(NZ)      !! y width
    real(kind=f), intent(in)                :: zc(NZ)      !! z at center
    real(kind=f), intent(in)                :: zl(NZ+1)    !! z at edge
    real(kind=f), intent(in)                :: p(NZ)       !! pressure at center [Pa]
    real(kind=f), intent(in)                :: pl(NZ+1)    !! presssure at edge [Pa]
    real(kind=f), intent(in)                :: t(NZ)       !! temperature at center [K]
    integer, intent(out)                    :: rc          !! return code, negative indicates failure
    real(kind=f), intent(in) , optional     :: qh2o(NZ)    !! specific humidity at center [mmr]
    real(kind=f), intent(in) , optional     :: relhum(NZ)  !! relative humidity at center [fraction]
    real(kind=f), intent(in) , optional     :: told(NZ)    !! previous temperature at center [K]

    integer                                 :: iz
    real(kind=f)                            :: rvap
    real(kind=f)                            :: pvap_liq
    real(kind=f)                            :: pvap_ice
    real(kind=f)                            :: gc_cgs  

    ! Assume success.
    rc = RC_OK

    ! Save the defintion of the number of comonents involved in the microphysics.
    cstate%carma => carma_ptr

    ! Save the model timing.
    cstate%time       = time
    cstate%dtime_orig = dtime
    cstate%dtime      = dtime
    cstate%nretries   = 0
    
    ! Save the grid dimensions.
    cstate%NZ   = NZ
    cstate%NZP1 = NZ+1
    
    ! Save the grid definition.
    cstate%igridv = igridv
    cstate%igridh = igridh
    
    ! Store away the grid location information.
    cstate%lat  = lat
    cstate%lon  = lon
    
    ! Allocate all the dynamic variables related to state.
    call CARMASTATE_Allocate(cstate, rc)
    if (rc < 0) return
    
    cstate%xc(:)  = xc(:)
    cstate%dx(:)  = dx(:)
    cstate%yc(:)  = yc(:)
    cstate%dy(:)  = dy(:)        
    cstate%zc(:)  = zc(:)
    cstate%zl(:)  = zl(:)

    ! Store away the grid state, doing any necessary unit conversions from MKS to CGS.
    cstate%p(:)  = p(:)  * RPA2CGS    
    cstate%pl(:) = pl(:) * RPA2CGS    
    cstate%t(:)  = t(:)
    
    cstate%pcd(:,:,:)     = 0._f
    
    if (carma_ptr%do_substep) then
      if (present(told)) then
        cstate%told(:) = told
      else
        if (carma_ptr%do_print) write(carma_ptr%LUNOPRT,*) "CARMASTATE_Create: Error - Need to specify told when substepping."
        rc = RC_ERROR
        
        return
      end if
    end if 
    
    ! Calculate the metrics, ...
    ! if Cartesian coordinates were specifed, then the units need to be converted
    ! from MKS to CGS.
    if (cstate%igridh == I_CART) then
      cstate%xc = cstate%xc * RM2CGS
      cstate%dx = cstate%dx * RM2CGS
      cstate%yc = cstate%yc * RM2CGS
      cstate%dy = cstate%dy * RM2CGS
    end if
    
    if (cstate%igridv == I_CART) then
      cstate%zc = cstate%zc * RM2CGS
      cstate%zl = cstate%zl * RM2CGS
    end if
    
    ! Initialize the state of the atmosphere.
    call setupatm(carma_ptr, cstate, carma_ptr%do_fixedinit, rc)
    if (rc < 0) return
    
    ! Set the realtive humidity. If necessary, it will be calculated from
    ! the specific humidity.
    if (present(relhum)) then
      cstate%relhum(:) = relhum(:)
    else if (present(qh2o)) then
    
      ! Define gas constant for this gas
      rvap = RGAS/WTMOL_H2O

      ! Calculate relative humidity
      do iz = 1, NZ
        call vaporp_h2o_murphy2005(carma_ptr, cstate, iz, rc, pvap_liq, pvap_ice)
        if (rc < 0) return

        gc_cgs = qh2o(iz)*cstate%rhoa_wet(iz) / (cstate%zmet(iz)*cstate%xmet(iz)*cstate%ymet(iz))
        cstate%relhum(iz) = ( gc_cgs * rvap * t(iz)) / pvap_liq
      enddo
    end if
    
    ! Need for vertical transport.
    !
    ! NOTE: How should these be set? Optional parameters?
    if (carma_ptr%do_vtran) then
      cstate%ftoppart(:,:) = 0._f
      cstate%fbotpart(:,:) = 0._f
      cstate%pc_topbnd(:,:) = 0._f
      cstate%pc_botbnd(:,:) = 0._f
    end if
        
    return
  end subroutine CARMASTATE_Create


  !! Create the CARMASTATE object, which contains information about the
  !! atmospheric state.
  !! 
  !! This call is similar to CARMASTATE_Create, but differs in that all the
  !! initialization happens here based on the the fixed state inofrmation provided rather
  !! than occurring in CARMASTATE_Step.
  !!
  !! This call should be done before CARMASTATE_Create when do_fixedinit has been
  !! specified. The temperatures and pressures specified here should be the reference
  !! state used for all columns, not an actual column from the model. This approach
  !! should not be used if particle swelling occurs, since the initialization needs to
  !! be recalculated based upon the wet radius.
  !!
  !! CARMASTATE_Create should still be called again after this call with the actual
  !! column of state information from the model. The initialization will be done once 
  !! from the reference state, but the microphysical calculations will be done on the
  !! model state. Multiple CARMASTATE_Create ... CARMASTATE_Step calls can be done
  !! before a CARMASTATE_Destroy. This reduces the amount of memory allocations and
  !! when used with do_fixedinit, reduces the amount of time spent initializing.
  !!
  !! @author Chuck Bardeen
  !! @version June-2010
  !! @see CARMA_Create
  !! @see CARMA_Initialize
  !! @see CARMASTATE_Destroy
  subroutine CARMASTATE_CreateFromReference(cstate, carma_ptr, time, dtime, NZ, igridv, igridh,  &
      lat, lon, xc, dx, yc, dy, zc, zl, p, pl, t, rc)
    type(carmastate_type), intent(inout)    :: cstate      !! the carma state object
    type(carma_type), pointer, intent(in)   :: carma_ptr   !! (in) the carma object
    real(kind=f), intent(in)                :: time        !! the model time [s]
    real(kind=f), intent(in)                :: dtime       !! the timestep size [s]
    integer, intent(in)                     :: NZ          !! the number of vertical grid points
    integer, intent(in)                     :: igridv      !! vertical grid type
    integer, intent(in)                     :: igridh      !! horizontal grid type
    real(kind=f), intent(in)                :: lat         !! latitude at center [degrees north]
    real(kind=f), intent(in)                :: lon         !! longitude at center [degrees east]
    real(kind=f), intent(in)                :: xc(NZ)      !! x at center
    real(kind=f), intent(in)                :: dx(NZ)      !! ix width
    real(kind=f), intent(in)                :: yc(NZ)      !! y at center
    real(kind=f), intent(in)                :: dy(NZ)      !! y width
    real(kind=f), intent(in)                :: zc(NZ)      !! z at center
    real(kind=f), intent(in)                :: zl(NZ+1)    !! z at edge
    real(kind=f), intent(in)                :: p(NZ)       !! pressure at center [Pa]
    real(kind=f), intent(in)                :: pl(NZ+1)    !! presssure at edge [Pa]
    real(kind=f), intent(in)                :: t(NZ)       !! temperature at center [K]
    integer, intent(out)                    :: rc          !! return code, negative indicates failure

    integer                                 :: iz
    real(kind=f)                            :: rvap
    real(kind=f)                            :: pvap_liq
    real(kind=f)                            :: pvap_ice
    real(kind=f)                            :: gc_cgs  

    ! Assume success.
    rc = RC_OK

    ! Save the defintion of the number of comonents involved in the microphysics.
    cstate%carma => carma_ptr

    ! Save the model timing.
    cstate%time       = time
    cstate%dtime_orig = dtime
    cstate%dtime      = dtime
    cstate%nretries   = 0
    
    ! Save the grid dimensions.
    cstate%NZ   = NZ
    cstate%NZP1 = NZ+1
    
    ! Save the grid definition.
    cstate%igridv = igridv
    cstate%igridh = igridh
    
    ! Store away the grid location information.
    cstate%lat  = lat
    cstate%lon  = lon
    
    ! Allocate all the dynamic variables related to state.
    call CARMASTATE_Allocate(cstate, rc)
    if (rc < 0) return
    
    cstate%xc(:)  = xc(:)
    cstate%dx(:)  = dx(:)
    cstate%yc(:)  = yc(:)
    cstate%dy(:)  = dy(:)        
    cstate%zc(:)  = zc(:)
    cstate%zl(:)  = zl(:)

    ! Store away the grid state, doing any necessary unit conversions from MKS to CGS.
    cstate%p(:)  = p(:)  * RPA2CGS    
    cstate%pl(:) = pl(:) * RPA2CGS    
    cstate%t(:)  = t(:)
    
    cstate%pcd(:,:,:)     = 0._f
    
    ! Calculate the metrics, ...
    ! if Cartesian coordinates were specifed, then the units need to be converted
    ! from MKS to CGS.
    if (cstate%igridh == I_CART) then
      cstate%xc = cstate%xc * RM2CGS
      cstate%dx = cstate%dx * RM2CGS
      cstate%yc = cstate%yc * RM2CGS
      cstate%dy = cstate%dy * RM2CGS
    end if
    
    if (cstate%igridv == I_CART) then
      cstate%zc = cstate%zc * RM2CGS
      cstate%zl = cstate%zl * RM2CGS
    end if
    
    ! Initialize the state of the atmosphere.
    call setupatm(carma_ptr, cstate, .false., rc)
    if (rc < 0) return

    ! Need for vertical transport.
    !
    ! NOTE: How should these be set? Optional parameters?
    if (carma_ptr%do_vtran) then
      cstate%ftoppart(:,:) = 0._f
      cstate%fbotpart(:,:) = 0._f
      cstate%pc_topbnd(:,:) = 0._f
      cstate%pc_botbnd(:,:) = 0._f
    end if
    
    
    ! Now do the initialization that is normally done in CARMASTATE_Step. However
    ! here it is done using the reference atmosphere.
    
    ! Determine the particle densities.
    call rhopart(cstate%carma, cstate, rc)
    if (rc < 0) return

    ! If configured for fixed initialization, then we will lose some accuracy
    ! in the calculation of the fall velocities, growth kernels, ... and in return
    ! will gain a significant performance by not having to initialize as often.
  
    ! Initialize the vertical transport.
    if (cstate%carma%do_vtran .or. cstate%carma%do_coag .or. cstate%carma%do_grow) then
      call setupvf(cstate%carma, cstate, rc)
      
      if (cstate%carma%do_vdiff) then
        call setupbdif(cstate%carma, cstate, rc)
      end if
    end if

    ! Intialize the nucleation, growth and evaporation.      
    if (cstate%carma%do_grow)  then
      call setupgrow(cstate%carma, cstate, rc)
      if (rc < 0) return

      call setupgkern(cstate%carma, cstate, rc)
      if (rc < 0) return
      
       call setupnuc(cstate%carma, cstate, rc)
      if (rc < 0) return
    end if
    
    ! Initialize the coagulation.
    if (cstate%carma%do_coag) then
      call setupckern(cstate%carma, cstate, rc)
      if (rc < 0) return
    end if
    
    return
  end subroutine CARMASTATE_CreateFromReference


  subroutine CARMASTATE_Allocate(cstate, rc)
    type(carmastate_type), intent(inout)  :: cstate
    integer, intent(out)                  :: rc
    
    ! Local Variables
    integer                               :: ier
    integer                               :: NZ
    integer                               :: NZP1
    integer                               :: NGROUP
    integer                               :: NELEM
    integer                               :: NBIN
    integer                               :: NGAS
    
    ! Assume success.
    rc = RC_OK

    ! Check to see if the arrays are already allocated. If so, just reuse the
    ! existing allocations.
    
    ! Allocate the variables needed for setupatm.
    if (.not. (allocated(cstate%xmet))) then
    
      NZ      = cstate%NZ
      NZP1    = cstate%NZP1
      NGROUP  = cstate%carma%NGROUP
      NELEM   = cstate%carma%NELEM
      NBIN    = cstate%carma%NBIN
      NGAS    = cstate%carma%NGAS
    
      allocate( &
        cstate%xmet(NZ), &
        cstate%ymet(NZ), &
        cstate%zmet(NZ), &
        cstate%zmetl(NZP1), &
        cstate%xc(NZ), &
        cstate%yc(NZ), &
        cstate%zc(NZ), &
        cstate%dx(NZ), &
        cstate%dy(NZ), &
        cstate%dz(NZ), &
        cstate%zl(NZP1), &
        cstate%pc(NZ,NBIN,NELEM), &
        cstate%pcd(NZ,NBIN,NELEM), &
        cstate%pc_surf(NBIN,NELEM), &
        cstate%gc(NZ,NGAS), &
        cstate%cldfrc(NZ), &
        cstate%rhcrit(NZ), &
        cstate%rhop(NZ,NBIN,NGROUP), &
        cstate%r_wet(NZ,NBIN,NGROUP), &
        cstate%rhop_wet(NZ,NBIN,NGROUP), &
        cstate%rhoa(NZ), &
        cstate%rhoa_wet(NZ), &
        cstate%t(NZ), &
        cstate%p(NZ), &
        cstate%pl(NZP1), &
        cstate%relhum(NZ), &
        cstate%rmu(NZ), &
        cstate%thcond(NZ), &
        cstate%dpc_sed(NBIN,NELEM), &
        cstate%pconmax(NZ,NGROUP), &
        cstate%pcl(NZ,NBIN,NELEM), &
        stat=ier)
      if (ier /= 0) then
        if (cstate%carma%do_print) write(cstate%carma%LUNOPRT, *) "CARMASTATE_Allocate::ERROR allocating atmosphere arrays, status=", ier
        rc = RC_ERROR
        return
      end if
      
      cstate%relhum(:)      = 0._f
      cstate%pc(:,:,:)      = 0._f
      cstate%pcd(:,:,:)     = 0._f
      cstate%pc_surf(:,:)   = 0._f
      cstate%cldfrc(:)      = 1._f
      cstate%rhcrit(:)      = 1._f
      
      ! Allocate the last fields if they are needed for substepping.
      if (cstate%carma%do_substep) then
        allocate( &
          cstate%gcl(NZ,NGAS), &
          cstate%d_gc(NZ,NGAS), &
          cstate%told(NZ), &
          cstate%d_t(NZ), &
          cstate%zsubsteps(NZ), &
          stat=ier)
        if (ier /= 0) then
          if (cstate%carma%do_print) write(cstate%carma%LUNOPRT, *) "CARMASTATE_Allocate::ERROR allocating stepping arrays, status=", ier
          rc = RC_ERROR
          return
        endif
      
        ! Initialize
        cstate%gcl(:,:)     = 0._f
        cstate%d_gc(:,:)    = 0._f
        cstate%told(:)      = 0._f
        cstate%d_t(:)       = 0._f
        cstate%zsubsteps(:) = 0._f

        ! When substepping is enabled, we want to initialize these statistics once for
        ! the life of the object.
        cstate%max_nsubstep = 0
        cstate%max_nretry   = 0._f
        cstate%nstep        = 0._f
        cstate%nsubstep     = 0
        cstate%nretry       = 0._f
      endif

      
      ! Allocate the variables needed for setupvf.
      !
      ! NOTE: Coagulation and dry deposition also need bpm, vf and re.
      if (cstate%carma%do_vtran .or. cstate%carma%do_coag .or. cstate%carma%do_grow .or. cstate%carma%do_drydep) then
        allocate( &
          cstate%bpm(NZ,NBIN,NGROUP), &
          cstate%vf(NZP1,NBIN,NGROUP), &
          cstate%re(NZ,NBIN,NGROUP), &
          cstate%dkz(NZP1,NBIN,NGROUP), &
          cstate%ftoppart(NBIN,NELEM), &
          cstate%fbotpart(NBIN,NELEM), &
          cstate%pc_topbnd(NBIN,NELEM), &
          cstate%pc_botbnd(NBIN,NELEM), &
          cstate%vd(NBIN, NGROUP), &
          stat=ier)
        if (ier /= 0) then
          if (cstate%carma%do_print) write(cstate%carma%LUNOPRT, *) "CARMASTATE_Allocate::ERROR allocating vertical transport arrays, status=", ier
          rc = RC_ERROR
          return
        endif

        ! Initialize
        cstate%bpm(:,:,:) = 0._f
        cstate%vf(:,:,:) = 0._f
        cstate%re(:,:,:) = 0._f
        cstate%dkz(:,:,:) = 0._f
        cstate%ftoppart(:,:) = 0._f
        cstate%fbotpart(:,:) = 0._f
        cstate%pc_topbnd(:,:) = 0._f
        cstate%pc_botbnd(:,:) = 0._f
        cstate%vd(:, :) = 0._f
      end if
      
      
      
      if (cstate%carma%NGAS > 0) then
        allocate( &
          cstate%pvapl(NZ,NGAS), &
          cstate%pvapi(NZ,NGAS), &
          cstate%supsatl(NZ,NGAS), &
          cstate%supsati(NZ,NGAS), &
          cstate%supsatlold(NZ,NGAS), &
          cstate%supsatiold(NZ,NGAS), &
          stat=ier)
        if (ier /= 0) then
          if (cstate%carma%do_print) write(cstate%carma%LUNOPRT, *) "CARMASTATE_Allocate::ERROR allocating gas arrays, status=", ier
          rc = RC_ERROR
          return
        endif
      end if
 
      
      if (cstate%carma%do_grow) then
        allocate( &
          cstate%diffus(NZ,NGAS), &
          cstate%rlhe(NZ,NGAS), &
          cstate%rlhm(NZ,NGAS), &
          cstate%surfctwa(NZ), &
          cstate%surfctiw(NZ), &
          cstate%surfctia(NZ), &
          cstate%akelvin(NZ,NGAS), &
          cstate%akelvini(NZ,NGAS), &
          cstate%ft(NZ,NBIN,NGROUP), &
          cstate%gro(NZ,NBIN,NGROUP),  &
          cstate%gro1(NZ,NBIN,NGROUP),  &
          cstate%gro2(NZ,NGROUP),  &
          cstate%scrit(NZ,NBIN,NGROUP), &
          cstate%rnuclg(NBIN,NGROUP,NGROUP),&
          cstate%rnucpe(NBIN,NELEM), &
          cstate%pc_nucl(NZ,NBIN,NELEM), &
          cstate%growpe(NBIN,NELEM), &
          cstate%evappe(NBIN,NELEM), &
          cstate%evcore(NELEM), &
          cstate%growlg(NBIN,NGROUP), &
          cstate%evaplg(NBIN,NGROUP), &
          cstate%gasprod(NGAS), &
          cstate%qrad(NZ,NBIN,NGROUP), &
          cstate%cmf(NBIN,NGROUP), &
          cstate%totevap(NBIN,NGROUP), &
          stat=ier)
        if (ier /= 0) then
          if (cstate%carma%do_print) write(cstate%carma%LUNOPRT, *) "CARMASTATE_Allocate::ERROR allocating growth arrays, status=", ier
          rc = RC_ERROR
          return
        endif
        
        cstate%qrad(:,:,:) = 0._f
      end if
      
      if (cstate%carma%do_coag) then
        allocate( &
          cstate%coaglg(NZ,NBIN,NGROUP), &
          cstate%coagpe(NZ,NBIN,NELEM), &
          cstate%ckernel(NZ,NBIN,NBIN,NGROUP,NGROUP), &
          cstate%pkernel(NZ,NBIN,NBIN,NGROUP,NGROUP,NGROUP,6), &
          stat = ier)
        if (ier /= 0) then
          if (cstate%carma%do_print) write(cstate%carma%LUNOPRT, *) "CARMASTATE_Allocate::ERROR allocating coag arrays, status=", ier
          rc = RC_ERROR
          return
        end if

        ! Initialize
        cstate%coaglg(:,:,:) = 0._f
        cstate%coagpe(:,:,:) = 0._f
        cstate%ckernel(:,:,:,:,:) = 0._f
        cstate%pkernel(:,:,:,:,:,:,:) = 0._f
      end if
    end if
    
    return
  end subroutine CARMASTATE_Allocate
    

  !! The routine should be called when the carma state object is no longer needed.
  !! It deallocates any memory allocations made by CARMA during CARMASTATE_Create(), 
  !! and failure to call this routine could result in memory leaks.
  !!
  !! @author Chuck Bardeen
  !! @version Feb-2009
  !! @see CARMASTATE_Create
  subroutine CARMASTATE_Destroy(cstate, rc)
    type(carmastate_type), intent(inout)    :: cstate
    integer, intent(out)                    :: rc
    
    ! Local variables
    integer   :: ier
    
    ! Assume success.
    rc = RC_OK

    ! Check to see if the arrays are already allocated. If so, deallocate them.

    ! Allocate the variables needed for setupatm.
    if (allocated(cstate%xmet)) then
    
      deallocate( &
        cstate%xmet, &
        cstate%ymet, &
        cstate%zmet, &
        cstate%zmetl, &
        cstate%xc, &
        cstate%yc, &
        cstate%zc, &
        cstate%dx, &
        cstate%dy, &
        cstate%dz, &
        cstate%zl, &
        cstate%pc, &
        cstate%pcd, &
        cstate%pc_surf, &
        cstate%gc, &
        cstate%cldfrc, &
        cstate%rhcrit, &
        cstate%rhop, &
        cstate%r_wet, &
        cstate%rhop_wet, &
        cstate%rhoa, &
        cstate%rhoa_wet, &
        cstate%t, &
        cstate%p, &
        cstate%pl, &
        cstate%relhum, &
        cstate%rmu, &
        cstate%thcond, &
        cstate%dpc_sed, &
        cstate%pconmax, &
        cstate%pcl, &
        stat=ier)
      if (ier /= 0) then
        if (cstate%carma%do_print) write(cstate%carma%LUNOPRT, *) "CARMASTATE_Destroy::ERROR deallocating atmosphere arrays, status=", ier
        rc = RC_ERROR
        return
      end if
      
      ! Allocate the last fields if they are needed for substepping stepping.
      if (allocated(cstate%gcl)) then
        deallocate( &
          cstate%gcl, &
          cstate%d_gc, &
          cstate%told, &
          cstate%d_t, &
          cstate%zsubsteps, &
          stat=ier)
        if (ier /= 0) then
          if (cstate%carma%do_print) write(cstate%carma%LUNOPRT, *) "CARMASTATE_Destroy::ERROR deallocating stepping arrays, status=", ier
          rc = RC_ERROR
          return
        endif
      endif
      
      ! Allocate the variables needed for setupvf.
      !
      ! NOTE: Coagulation also needs bpm, vf and re.
      if (allocated(cstate%bpm)) then
        deallocate( &
          cstate%bpm, &
          cstate%vf, &
          cstate%re, &
          cstate%dkz, &
          cstate%ftoppart, &
          cstate%fbotpart, &
          cstate%pc_topbnd, &
          cstate%pc_botbnd, &
          stat=ier)
        if (ier /= 0) then
          if (cstate%carma%do_print) write(cstate%carma%LUNOPRT, *) "CARMASTATE_Destroy::ERROR deallocating vertical transport arrays, status=", ier
          rc = RC_ERROR
          return
        endif
      end if
      
      if (allocated(cstate%diffus)) then
        deallocate( &
          cstate%diffus, &
          cstate%rlhe, &
          cstate%rlhm, &
          cstate%surfctwa, &
          cstate%surfctiw, &
          cstate%surfctia, &
          cstate%akelvin, &
          cstate%akelvini, &
          cstate%ft, &
          cstate%gro, &
          cstate%gro1, &
          cstate%gro2, &
          cstate%scrit, &
          cstate%rnuclg,&
          cstate%rnucpe, &
          cstate%pc_nucl, &
          cstate%growpe, &
          cstate%evappe, &
          cstate%evcore, &
          cstate%growlg, &
          cstate%evaplg, &
          cstate%gasprod, &
          cstate%qrad, &
          cstate%cmf, &
          cstate%totevap, &
          stat=ier)
        if (ier /= 0) then
          if (cstate%carma%do_print) write(cstate%carma%LUNOPRT, *) "CARMASTATE_Destroy::ERROR deallocating growth arrays, status=", ier
          rc = RC_ERROR
          return
        endif
      end if
      
      if (allocated(cstate%pvapl)) then
        deallocate( &
          cstate%pvapl, &
          cstate%pvapi, &
          cstate%supsatl, &
          cstate%supsati, &
          cstate%supsatlold, &
          cstate%supsatiold, &
          stat=ier)
        if (ier /= 0) then
          if (cstate%carma%do_print) write(cstate%carma%LUNOPRT, *) "CARMASTATE_Destroy::ERROR deallocating gas arrays, status=", ier
          rc = RC_ERROR
          return
        endif
      end if
      
      if (allocated(cstate%coaglg)) then
        deallocate( &
          cstate%coaglg, &
          cstate%coagpe, &
          cstate%ckernel, &
          cstate%pkernel, &
          stat = ier)
        if (ier /= 0) then
          if (cstate%carma%do_print) write(cstate%carma%LUNOPRT, *) "CARMASTATE_Destroy::ERROR deallocating coag arrays, status=", ier
          rc = RC_ERROR
          return
        end if
      end if
    end if
    
    return
  end subroutine CARMASTATE_Destroy


  !! The routine performs the main CARMA processing for one timestep of
  !! the parent model. The state variables should have all been set before
  !! calling CARMASTATE_Step(). When this routine returns, the state will
  !! have been updated to reflect the changes from the CARMA microphysics.
  !! If tendencies are desired, then the difference between the final and
  !! initial state will need to be computed by the caller.
  !!
  !! @author Chuck Bardeen
  !! @version Feb-2009
  subroutine CARMASTATE_Step(cstate, rc, cldfrc, rhcrit, surfric, ram, landfrac, ocnfrac, icefrac)
    type(carmastate_type), intent(inout)  :: cstate
    integer, intent(out)                  :: rc
    real(kind=f), intent(in), optional    :: cldfrc(cstate%NZ)   !! cloud fraction [fraction]
    real(kind=f), intent(in), optional    :: rhcrit(cstate%NZ)   !! relative humidity for onset of liquid clouds [fraction]
    real(kind=f), intent(in), optional    :: surfric             !! surface friction velocity [m/s]
    real(kind=f), intent(in), optional    :: ram                 !! aerodynamic resistance [s/m]
    real(kind=f), intent(in), optional    :: landfrac            !! land fraction
    real(kind=f), intent(in), optional    :: ocnfrac             !! ocn fraction
    real(kind=f), intent(in), optional    :: icefrac             !! ice fraction

    
    integer                               :: iz     ! vertical index
    integer                               :: igas   ! gas index
    integer                               :: ielem
    integer                               :: ibin
    integer                               :: igroup
  
    ! Assume success.
    rc = RC_OK
    
    ! Store the cloud fraction if specified
    cstate%cldfrc(:) = 1._f
    cstate%rhcrit(:) = 1._f
    
    if (present(cldfrc)) cstate%cldfrc(:) = cldfrc(:)
    if (present(rhcrit)) cstate%rhcrit(:) = rhcrit(:)
    
    ! Determine the gas supersaturations.
    do iz = 1, cstate%NZ
      do igas = 1, cstate%carma%NGAS
        call supersat(cstate%carma, cstate, iz, igas, rc)
        if (rc < 0) return
      end do
    end do

    ! Determine the particle densities.
    call rhopart(cstate%carma, cstate, rc)
    if (rc < 0) return
    

    ! We have to hold off initialization until now, because the particle density
    ! (rhop) can not be determined until the particle masses are known (i.e. after
    ! CARMASTATE_SetBin), because rhop is used to determine the fall velocity.
    !
    ! NOTE: If configured for fixed initialization, then we will lose some accuracy
    ! in the calculation of the fall velocities, growth kernels, ... and in return
    ! will gain a significant performance by not having to initialize as often.
    if (.not. cstate%carma%do_fixedinit) then

      ! Initialize the vertical transport.
      if (cstate%carma%do_vtran .or. cstate%carma%do_coag .or. cstate%carma%do_grow) then
        call setupvf(cstate%carma, cstate, rc)

        if (cstate%carma%do_vdiff) then
          call setupbdif(cstate%carma, cstate, rc)
        end if
      end if
      
      ! intialize the dry deposition
      if (cstate%carma%do_drydep) then
        if (present(surfric) .and. present(ram) .and. present(landfrac) .and. present(ocnfrac) .and. present(icefrac)) then
        
          ! NOTE: Need to convert surfric and ram from mks to cgs units.
          call setupvdry(cstate%carma, cstate, surfric * 100._f, ram / 100._f, landfrac, ocnfrac, icefrac, rc)
          if (rc < RC_OK) return
        else
          write(cstate%carma%LUNOPRT, *) "CARMASTATE_Step: do_drydep requires that the optional inputs surfric, ram, landfrac, ocnfrac and icefrac be provided."
          rc = RC_ERROR
        end if
      end if
       
      ! Intialize the nucleation, growth and evaporation.      
      if (cstate%carma%do_grow)  then
        call setupgrow(cstate%carma, cstate, rc)
        if (rc < RC_OK) return
  
        call setupgkern(cstate%carma, cstate, rc)
        if (rc < RC_OK) return
        
         call setupnuc(cstate%carma, cstate, rc)
        if (rc < RC_OK) return
      end if
      
      ! Initialize the coagulation.
      if (cstate%carma%do_coag) then
        call setupckern(cstate%carma, cstate, rc)
        if (rc < RC_OK) return
      end if
    end if
    
    ! Calculate the impact of microphysics upon the state.
    call step(cstate%carma, cstate, rc)

    return
  end subroutine CARMASTATE_Step


  ! Query, Control and State I/O

  !! Gets the mass mixing ratio for the gas (igas). After a call to CARMA_Step(),
  !! the new mass mixing ratio of the gas can be retrieved.
  !!
  !! @author Chuck Bardeen
  !! @version Feb-2009
  !! @see CARMA_AddGas
  !! @see CARMA_GetGas
  !! @see CARMA_Step 
  !! @see CARMASTATE_SetGas
  subroutine CARMASTATE_Get(cstate, rc, max_nsubstep, max_nretry, nstep, nsubstep, nretry, zsubsteps)
    type(carmastate_type), intent(in)     :: cstate            !! the carma state object
    integer, intent(out)                  :: rc                !! return code, negative indicates failure
    integer, optional, intent(out)        :: max_nsubstep      !! maximum number of substeps in a step
    real(kind=f), optional, intent(out)   :: max_nretry        !! maximum number of retries in a step
    real(kind=f), optional, intent(out)   :: nstep             !! total number of steps taken
    integer, optional, intent(out)        :: nsubstep          !! total number of substeps taken
    real(kind=f), optional, intent(out)   :: nretry            !! total number of retries taken
    real(kind=f), optional, intent(out)   :: zsubsteps(cstate%NZ) !! number of substeps taken per vertical grid point
    
    ! Assume success.
    rc = RC_OK

    if (present(max_nsubstep)) max_nsubstep = cstate%max_nsubstep
    if (present(max_nretry))   max_nretry   = cstate%max_nretry
    if (present(nstep))        nstep        = cstate%nstep
    if (present(nsubstep))     nsubstep     = cstate%nsubstep
    if (present(nretry))       nretry       = cstate%nretry
    if (present(zsubsteps))    zsubsteps    = cstate%zsubsteps
    
    return
  end subroutine CARMASTATE_Get
  

  !! Gets the mass of the bins (ibin) for each particle element (ielem). After the
  !! CARMA_Step() call, new particle concentrations are determined. The number density
  !! and the nucleation rate are only calculated if the element is the number density
  !! element for the group.
  !!
  !! @author Chuck Bardeen
  !! @version Feb-2009
  !! @see CARMA_AddElement
  !! @see CARMA_AddGroup
  !! @see CARMA_Step 
  !! @see CARMASTATE_SetBin
  subroutine CARMASTATE_GetBin(cstate, ielem, ibin, mmr, rc, nmr, numberDensity, nucleationRate, r_wet, rhop_wet, surface)
    type(carmastate_type), intent(in)     :: cstate         !! the carma state object
    integer, intent(in)                   :: ielem          !! the element index
    integer, intent(in)                   :: ibin           !! the bin index
    real(kind=f), intent(out)             :: mmr(cstate%NZ) !! the bin mass mixing ratio [kg/kg]
    integer, intent(out)                  :: rc             !! return code negative indicates failure
    real(kind=f), optional, intent(out)   :: nmr(cstate%NZ) !! number mixing ratio [#/kg]
    real(kind=f), optional, intent(out)   :: numberDensity(cstate%NZ)  !! number density [#/cm3]
    real(kind=f), optional, intent(out)   :: nucleationRate(cstate%NZ) !! nucleation rate [1/cm3/s]
    real(kind=f), optional, intent(out)   :: r_wet(cstate%NZ)          !! wet particle radius [cm]
    real(kind=f), optional, intent(out)   :: rhop_wet(cstate%NZ)       !! wet particle density [g/cm3]
    real(kind=f), optional, intent(out)   :: surface        !! particle mass on the surface [kg/m2]
    
    integer                               :: ienconc        !! index of element that is the particle concentration for the group
    integer                               :: igroup         ! Group containing this bin

    ! Assume success.
    rc = RC_OK
    
    ! Determine the particle group for the bin.    
    igroup = cstate%carma%element(ielem)%igroup

    ! Make sure there are enough elements allocated.
    if (ielem > cstate%carma%NELEM) then
      if (cstate%carma%do_print) write(cstate%carma%LUNOPRT, *) "CARMASTATE_SetBin:: ERROR - The specifed element (", &
        ielem, ") is larger than the number of elements (", cstate%carma%NELEM, ")."
      rc = RC_ERROR
      return
    end if
    
    ! Make sure there are enough bins allocated.
    if (ibin > cstate%carma%NBIN) then
      if (cstate%carma%do_print) write(cstate%carma%LUNOPRT, *) "CARMA_SetBin:: ERROR - The specifed bin (", &
        ibin, ") is larger than the number of bins (", cstate%carma%NBIN, ")."
      rc = RC_ERROR
      return
    end if

    
    ! Use the specified mass mixing ratio and the air density to determine the mass
    ! of the particles in g/x/y/z.
    mmr(:) = cstate%pc(:, ibin, ielem) / cstate%rhoa_wet(:)


    ! Handle the special cases for different types of elements ...
    if ((cstate%carma%element(ielem)%itype == I_INVOLATILE) .or. (cstate%carma%element(ielem)%itype == I_VOLATILE)) then
      mmr(:) = mmr(:) * cstate%carma%group(igroup)%rmass(ibin)
    else if (cstate%carma%element(ielem)%itype == I_CORE2MOM) then
      mmr(:) = mmr(:) / cstate%carma%group(igroup)%rmass(ibin)
    end if
    
    ! If the number of particles in the group is less than the minimum value represented
    ! by CARMA, then return and mmr of 0.0 for all elements.
    ienconc = cstate%carma%group(igroup)%ienconc
!    where (cstate%pc(:, ibin, ienconc) <= SMALL_PC) mmr(:) = 0.0_f


    ! Do they also want the mass flux of particles that sedimented to the surface?
    if (present(surface)) then
      
      ! Convert from g/cm2 to kg/m2
      surface = cstate%pc_surf(ibin, ielem) * 1e4_f / 1e3_f

      ! Handle the special cases for different types of elements ...
      if ((cstate%carma%element(ielem)%itype == I_INVOLATILE) .or. (cstate%carma%element(ielem)%itype == I_VOLATILE)) then
        surface = surface * cstate%carma%group(igroup)%rmass(ibin)
      else if (cstate%carma%element(ielem)%itype == I_CORE2MOM) then
        surface = surface / cstate%carma%group(igroup)%rmass(ibin)
      end if
    end if
    
    ! If this is the partcile # element, then determine some other statistics.
    if (ienconc == ielem) then
      if (present(nmr))           nmr(:)             = (cstate%pc(:, ibin, ielem) / cstate%rhoa_wet(:)) * 1000._f
      if (present(numberDensity)) numberDensity(:)   = cstate%pc(:, ibin, ielem) / (cstate%xmet(:)*cstate%ymet(:)*cstate%zmet(:))
      if (present(r_wet))         r_wet(:)           = cstate%r_wet(:, ibin, igroup)
      if (present(rhop_wet))      rhop_wet(:)        = cstate%rhop_wet(:, ibin, igroup)

      if (cstate%carma%do_grow) then
        if (present(nucleationRate)) nucleationRate(:) = cstate%pc_nucl(:, ibin, ielem) / (cstate%xmet(:)*cstate%ymet(:)*cstate%zmet(:)) / cstate%dtime
      else
        if (present(nucleationRate)) nucleationRate(:) = CAM_FILL
      end if
    else
      if (present(nmr))            nmr(:)             = CAM_FILL
      if (present(numberDensity))  numberDensity(:)   = CAM_FILL
      if (present(nucleationRate)) nucleationRate(:)  = CAM_FILL
    end if
   
    return
  end subroutine CARMASTATE_GetBin
  

  !! Gets the mass of the detrained condensate for the bins (ibin) for each particle
  !! element (ielem) in the grid.
  !!
  !!
  !! @author Chuck Bardeen
  !! @version Feb-2009
  !! @see CARMA_AddElement
  !! @see CARMA_AddGroup
  !! @see CARMA_Step 
  !! @see CARMASTATE_SetDetrain
  subroutine CARMASTATE_GetDetrain(cstate, ielem, ibin, mmr, rc, nmr, numberDensity, r_wet, rhop_wet)
    type(carmastate_type), intent(in)     :: cstate         !! the carma state object
    integer, intent(in)                   :: ielem          !! the element index
    integer, intent(in)                   :: ibin           !! the bin index
    real(kind=f), intent(out)             :: mmr(cstate%NZ) !! the bin mass mixing ratio [kg/kg]
    integer, intent(out)                  :: rc             !! return code negative indicates failure
    real(kind=f), optional, intent(out)   :: nmr(cstate%NZ) !! number mixing ratio [#/kg]
    real(kind=f), optional, intent(out)   :: numberDensity(cstate%NZ)  !! number density [#/cm3]
    real(kind=f), optional, intent(out)   :: r_wet(cstate%NZ)          !! wet particle radius [cm]
    real(kind=f), optional, intent(out)   :: rhop_wet(cstate%NZ)       !! wet particle density [g/cm3]
    
    integer                               :: ienconc        !! index of element that is the particle concentration for the group
    integer                               :: igroup         ! Group containing this bin

    ! Assume success.
    rc = RC_OK
    
    ! Determine the particle group for the bin.    
    igroup = cstate%carma%element(ielem)%igroup

    ! Make sure there are enough elements allocated.
    if (ielem > cstate%carma%NELEM) then
      if (cstate%carma%do_print) write(cstate%carma%LUNOPRT, *) "CARMASTATE_SetDetrain:: ERROR - The specifed element (", &
        ielem, ") is larger than the number of elements (", cstate%carma%NELEM, ")."
      rc = RC_ERROR
      return
    end if
    
    ! Make sure there are enough bins allocated.
    if (ibin > cstate%carma%NBIN) then
      if (cstate%carma%do_print) write(cstate%carma%LUNOPRT, *) "CARMA_SetDetrainin:: ERROR - The specifed bin (", &
        ibin, ") is larger than the number of bins (", cstate%carma%NBIN, ")."
      rc = RC_ERROR
      return
    end if

    
    ! Use the specified mass mixing ratio and the air density to determine the mass
    ! of the particles in g/x/y/z.
    mmr(:) = cstate%pcd(:, ibin, ielem) / cstate%rhoa_wet(:)


    ! Handle the special cases for different types of elements ...
    if ((cstate%carma%element(ielem)%itype == I_INVOLATILE) .or. (cstate%carma%element(ielem)%itype == I_VOLATILE)) then
      mmr(:) = mmr(:) * cstate%carma%group(igroup)%rmass(ibin)
    else if (cstate%carma%element(ielem)%itype == I_CORE2MOM) then
      mmr(:) = mmr(:) / cstate%carma%group(igroup)%rmass(ibin)
    end if
       
    ! If this is the partcile # element, then determine some other statistics.
    ienconc = cstate%carma%group(igroup)%ienconc
    if (ienconc == ielem) then
      if (present(nmr))           nmr(:)             = (cstate%pcd(:, ibin, ielem) / cstate%rhoa_wet(:)) * 1000._f
      if (present(numberDensity)) numberDensity(:)   = cstate%pcd(:, ibin, ielem) / (cstate%xmet(:)*cstate%ymet(:)*cstate%zmet(:))
      if (present(r_wet))         r_wet(:)           = cstate%r_wet(:, ibin, igroup)
      if (present(rhop_wet))      rhop_wet(:)        = cstate%rhop_wet(:, ibin, igroup)
    else
      if (present(nmr))            nmr(:)             = CAM_FILL
      if (present(numberDensity))  numberDensity(:)   = CAM_FILL
    end if
    
   return
  end subroutine CARMASTATE_GetDetrain
  

  !! Gets the mass mixing ratio for the gas (igas). After a call to CARMA_Step(),
  !! the new mass mixing ratio of the gas can be retrieved.
  !!
  !! @author Chuck Bardeen
  !! @version Feb-2009
  !! @see CARMA_AddGas
  !! @see CARMA_GetGas
  !! @see CARMA_Step 
  !! @see CARMASTATE_SetGas
  subroutine CARMASTATE_GetGas(cstate, igas, mmr, rc, satice, satliq)
    type(carmastate_type), intent(in)     :: cstate            !! the carma state object
    integer, intent(in)                   :: igas              !! the gas index
    real(kind=f), intent(out)             :: mmr(cstate%NZ)    !! the gas mass mixing ratio [kg/kg]
    integer, intent(out)                  :: rc                !! return code, negative indicates failure
    real(kind=f), optional, intent(out)   :: satice(cstate%NZ) !! the gas saturation wrt ice
    real(kind=f), optional, intent(out)   :: satliq(cstate%NZ) !! the gas saturation wrt liquid

    ! Assume success.
    rc = RC_OK

    ! Make sure there are enough gases allocated.
    if (igas > cstate%carma%NGAS) then
      if (cstate%carma%do_print) write(cstate%carma%LUNOPRT, *) "CARMASTATE_GetGas:: ERROR - The specifed gas (", &
        igas, ") is larger than the number of gases (", cstate%carma%NGAS, ")."
      rc = RC_ERROR
      return
    end if
    
    ! Use the specified mass mixing ratio and the air density to determine the mass
    ! of the gas in g/x/y/z.
    mmr(:) = cstate%gc(:, igas) / cstate%rhoa_wet(:)

    if (present(satice)) satice(:) = cstate%supsati(:, igas) + 1._f
    if (present(satice)) satliq(:) = cstate%supsatl(:, igas) + 1._f
    
    return
  end subroutine CARMASTATE_GetGas
  

  !! Gets information about the state of the atmosphere. After the CARMA_Step() call,
  !! a new atmospheric state is determined.
  !!
  !! @author Chuck Bardeen
  !! @version Feb-2009
  !! @see CARMA_Step 
  !! @see CARMASTATE_Create
  subroutine CARMASTATE_GetState(cstate, rc, t, p, rhoa_wet)
    type(carmastate_type), intent(in)     :: cstate              !! the carma state object
    integer, intent(out)                  :: rc                  !! return code, negative indicates failure
    real(kind=f), optional, intent(out)   :: t(cstate%NZ)        !! the air temperature [K]
    real(kind=f), optional, intent(out)   :: p(cstate%NZ)        !! the air pressure [Pa]
    real(kind=f), optional, intent(out)   :: rhoa_wet(cstate%NZ) !! air density [kg m-3]
    
    ! Assume success.
    rc = RC_OK

    ! Return the temperature, pressure, and/or density.
    if (present(t))         t(:) = cstate%t(:)
    
    ! DYNE -> Pa
    if (present(p))         p(:) = cstate%p(:) / RPA2CGS
    
    ! Convert rhoa from the scaled units to mks.
    if (present(rhoa_wet))  rhoa_wet(:) = (cstate%rhoa_wet(:) / (cstate%zmet(:)*cstate%xmet(:)*cstate%ymet(:))) * 1e6_f / 1e3_f
    
    return
  end subroutine CARMASTATE_GetState
  

  !! Sets the mass of the bins (ibin) for each particle element (ielem) in the grid.
  !! This call should be made after CARMASTATE_Create() and before CARMA_Step().
  !!
  !! @author Chuck Bardeen
  !! @version Feb-2009
  !! @see CARMA_AddBin
  !! @see CARMA_Step 
  !! @see CARMASTATE_GetBin
  subroutine CARMASTATE_SetBin(cstate, ielem, ibin, mmr, rc, surface)
    type(carmastate_type), intent(inout)  :: cstate         !! the carma state object
    integer, intent(in)                   :: ielem          !! the element index
    integer, intent(in)                   :: ibin           !! the bin index
    real(kind=f), intent(in)              :: mmr(cstate%NZ) !! the bin mass mixing ratio [kg/kg]
    integer, intent(out)                  :: rc             !! return code, negative indicates failure
    real(kind=f), optional, intent(in)    :: surface        !! particles mass on the surface [kg/m2]
    
    integer                               :: igroup         ! Group containing this bin

    ! Assume success.
    rc = RC_OK
    
    ! Determine the particle group for the bin.    
    igroup = cstate%carma%element(ielem)%igroup

    ! Make sure there are enough elements allocated.
    if (ielem > cstate%carma%NELEM) then
      if (cstate%carma%do_print) write(cstate%carma%LUNOPRT, *) "CARMASTATE_SetBin:: ERROR - The specifed element (", &
        ielem, ") is larger than the number of elements (", cstate%carma%NELEM, ")."
      rc = RC_ERROR
      return
    end if
    
    ! Make sure there are enough bins allocated.
    if (ibin > cstate%carma%NBIN) then
      if (cstate%carma%do_print) write(cstate%carma%LUNOPRT, *) "CARMASTATE_SetBin:: ERROR - The specifed bin (", &
        ibin, ") is larger than the number of bins (", cstate%carma%NBIN, ")."
      rc = RC_ERROR
      return
    end if
    
    ! Use the specified mass mixing ratio and the air density to determine the mass
    ! of the particles in g/x/y/z.
    cstate%pc(:, ibin, ielem) = mmr(:) * cstate%rhoa_wet(:)
    
    ! Handle the special cases for different types of elements ...
    if ((cstate%carma%element(ielem)%itype == I_INVOLATILE) .or. (cstate%carma%element(ielem)%itype == I_VOLATILE)) then
      cstate%pc(:, ibin, ielem) = cstate%pc(:, ibin, ielem) / cstate%carma%group(igroup)%rmass(ibin)
    else if (cstate%carma%element(ielem)%itype == I_CORE2MOM) then
      cstate%pc(:, ibin, ielem) = cstate%pc(:, ibin, ielem) * cstate%carma%group(igroup)%rmass(ibin)
    end if
    
    ! If they specified an initial mass of particles on the surface, then use that
    ! value.
    if (present(surface)) then
      
      ! Convert from g/cm2 to kg/m2
      cstate%pc_surf(ibin, ielem) = surface / 1e4_f * 1e3_f

      ! Handle the special cases for different types of elements ...
      if ((cstate%carma%element(ielem)%itype == I_INVOLATILE) .or. (cstate%carma%element(ielem)%itype == I_VOLATILE)) then
        cstate%pc_surf(ibin, ielem) = cstate%pc_surf(ibin, ielem) / cstate%carma%group(igroup)%rmass(ibin)
      else if (cstate%carma%element(ielem)%itype == I_CORE2MOM) then
        cstate%pc_surf(ibin, ielem) = cstate%pc_surf(ibin, ielem) * cstate%carma%group(igroup)%rmass(ibin)
      end if
    else
      cstate%pc_surf(ibin, ielem) = 0.0_f
    end if
        
    return
  end subroutine CARMASTATE_SetBin
  

  !! Sets the mass of the detrained condensate for the bins (ibin) for each particle
  !! element (ielem) in the grid. This call should be made after CARMASTATE_Create()
  !! and before CARMA_Step().
  !!
  !! @author Chuck Bardeen
  !! @version May-2010
  !! @see CARMA_AddBin
  !! @see CARMA_Step 
  !! @see CARMASTATE_GetDetrain
  subroutine CARMASTATE_SetDetrain(cstate, ielem, ibin, mmr, rc)
    type(carmastate_type), intent(inout)  :: cstate         !! the carma state object
    integer, intent(in)                   :: ielem          !! the element index
    integer, intent(in)                   :: ibin           !! the bin index
    real(kind=f), intent(in)              :: mmr(cstate%NZ) !! the bin mass mixing ratio [kg/kg]
    integer, intent(out)                  :: rc             !! return code, negative indicates failure
    
    integer                               :: igroup         ! Group containing this bin

    ! Assume success.
    rc = RC_OK
    
    ! Determine the particle group for the bin.    
    igroup = cstate%carma%element(ielem)%igroup

    ! Make sure there are enough elements allocated.
    if (ielem > cstate%carma%NELEM) then
      if (cstate%carma%do_print) write(cstate%carma%LUNOPRT, *) "CARMASTATE_SetDetrain:: ERROR - The specifed element (", &
        ielem, ") is larger than the number of elements (", cstate%carma%NELEM, ")."
      rc = RC_ERROR
      return
    end if
    
    ! Make sure there are enough bins allocated.
    if (ibin > cstate%carma%NBIN) then
      if (cstate%carma%do_print) write(cstate%carma%LUNOPRT, *) "CARMASTATE_SetDetrain:: ERROR - The specifed bin (", &
        ibin, ") is larger than the number of bins (", cstate%carma%NBIN, ")."
      rc = RC_ERROR
      return
    end if
    
    ! Use the specified mass mixing ratio and the air density to determine the mass
    ! of the particles in g/x/y/z.
    cstate%pcd(:, ibin, ielem) = mmr(:) * cstate%rhoa_wet(:)
    
    ! Handle the special cases for different types of elements ...
    if ((cstate%carma%element(ielem)%itype == I_INVOLATILE) .or. (cstate%carma%element(ielem)%itype == I_VOLATILE)) then
      cstate%pcd(:, ibin, ielem) = cstate%pcd(:, ibin, ielem) / cstate%carma%group(igroup)%rmass(ibin)
    else if (cstate%carma%element(ielem)%itype == I_CORE2MOM) then
      cstate%pcd(:, ibin, ielem) = cstate%pcd(:, ibin, ielem) * cstate%carma%group(igroup)%rmass(ibin)
    end if
        
    return
  end subroutine CARMASTATE_SetDetrain
  


  !! Sets the mass of the gas (igas) in the grid. This call should be made after
  !! CARMASTATE_Create() and before CARMA_Step().
  !!
  !! @author Chuck Bardeen
  !! @version Feb-2009
  !! @see CARMA_AddGas
  !! @see CARMA_GetGas
  !! @see CARMA_InitializeStep 
  !! @see CARMA_Step 
  subroutine CARMASTATE_SetGas(cstate, igas, mmr, rc, mmr_old, satice_old, satliq_old)
    type(carmastate_type), intent(inout)  :: cstate         !! the carma object
    integer, intent(in)                   :: igas           !! the gas index
    real(kind=f), intent(in)              :: mmr(cstate%NZ) !! the gas mass mixing ratio [kg/kg]
    integer, intent(out)                  :: rc             !! return code, negative indicates failure
    real(kind=f), intent(in), optional    :: mmr_old(cstate%NZ) !! the previous gas mass mixing ratio [kg/kg]
    real(kind=f), intent(inout), optional :: satice_old(cstate%NZ) !! the previous gas saturation wrt ice, calculates if -1
    real(kind=f), intent(inout), optional :: satliq_old(cstate%NZ) !! the previous gas saturation wrt liquid, calculates if -1
    
    real(kind=f)                          :: tnew(cstate%NZ)
    integer                               :: iz
    logical                               :: calculateOld
    
    ! Assume success.
    rc = RC_OK

    ! Make sure there are enough gases allocated.
    if (igas > cstate%carma%NGAS) then
      if (cstate%carma%do_print) write(cstate%carma%LUNOPRT, *) "CARMASTATE_SetGas:: ERROR - The specifed gas (", &
        igas, ") is larger than the number of gases (", cstate%carma%NGAS, ")."
      rc = RC_ERROR
      return
    end if
    
    if (cstate%carma%do_substep) then
      if (.not. present(mmr_old)) then
        if (cstate%carma%do_print) write(cstate%carma%LUNOPRT,*) "CARMASTATE_SetGas: Error - Need to specify mmr_old, satic_old, satliq_old when substepping."
        rc = RC_ERROR
        
        return
        
      else
        cstate%gcl(:, igas) = mmr_old(:) * cstate%rhoa_wet(:) * cstate%t(:) / cstate%told(:)
      
        ! A value of -1 for the saturation ratio means that it needs to be calculated from the old temperature
        ! and the old gc.
        !
        ! NOTE: This is typically just a problem for the first step, so we just need to get close.
        calculateOld = .false.
        if (present(satice_old) .and. present(satliq_old)) then
          if (any(satice_old(:) == -1._f) .or. any(satliq_old(:) == -1._f)) calculateOld = .true.
        else 
          calculateOld = .true.
        end if
        
        if (calculateOld) then
          
          ! This is a bit of a hack, because of the way CARMA has the vapor pressure and saturation
          ! routines implemented.
          
          ! Temporarily set the temperature and gc of to the old state
          
          tnew(:)      = cstate%t(:)
          cstate%t(:)  = cstate%told(:)
       
          cstate%gc(:, igas) = mmr_old(:) * cstate%rhoa_wet(:)
          
          do iz = 1, cstate%NZ
            call supersat(cstate%carma, cstate, iz, igas, rc)
            if (rc /= RC_OK) return
          
            if (present(satice_old)) then
              if (satice_old(iz) == -1._f) then
                cstate%supsatiold(iz, igas) = cstate%supsati(iz, igas)
              else
                cstate%supsatiold(iz, igas) = satice_old(iz) - 1._f
              endif
            else
              cstate%supsatiold(iz, igas) = cstate%supsati(iz, igas)
            end if
            
            if (present(satliq_old)) then
              if (satliq_old(iz) == -1._f) then
                cstate%supsatlold(iz, igas) = cstate%supsatl(iz, igas)
              else
                cstate%supsatlold(iz, igas) = satliq_old(iz) - 1._f
              endif
            else
              cstate%supsatlold(iz, igas) = cstate%supsatl(iz, igas)
            end if
          end do
          
          cstate%t(:) = tnew(:)
        
        else
          cstate%supsatiold(:, igas) = satice_old(:) - 1._f
          cstate%supsatlold(:, igas) = satliq_old(:) - 1._f
        end if
      end if
    end if

    ! Use the specified mass mixing ratio and the air density to determine the mass
    ! of the gas in g/x/y/z.
    cstate%gc(:, igas)  = mmr(:) * cstate%rhoa_wet(:)
    
    return
  end subroutine CARMASTATE_SetGas
  
  
  !! Sets information about the state of the atmosphere.
  !!
  !! @author Chuck Bardeen
  !! @version Feb-2009
  !! @see CARMA_Step 
  !! @see CARMASTATE_Create
  subroutine CARMASTATE_SetState(cstate, rc, t, rhoa_wet)
    type(carmastate_type), intent(inout)  :: cstate              !! the carma state object
    integer, intent(out)                  :: rc                  !! return code, negative indicates failure
    real(kind=f), optional, intent(in)    :: t(cstate%NZ)        !! the air temperature [K]
    real(kind=f), optional, intent(in)    :: rhoa_wet(cstate%NZ) !! air density [kg m-3]
    
    ! Assume success.
    rc = RC_OK

    ! Return the temperature or density.
    if (present(t))         cstate%t(:) = t(:)
    
    ! Convert rhoa from mks to the scaled units.
    if (present(rhoa_wet))  cstate%rhoa_wet(:) = (rhoa_wet(:) * (cstate%zmet(:)*cstate%xmet(:)*cstate%ymet(:))) / 1e6_f * 1e3_f
    
    return
  end subroutine CARMASTATE_SetState
end module
