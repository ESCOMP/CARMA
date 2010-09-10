!! The CARMA module contains an interface to the Community Aerosol and Radiation
!! Model for Atmospheres (CARMA) bin microphysical model [Turco et al. 1979; 
!! Toon et al. 1988]. This implementation has been customized to work within
!! other model frameworks, so although it can be provided with an array of
!! columns, it does not do horizontal transport and just does independent 1-D
!! calculations upon each column.
!!
!! The typical usage for the CARMA and CARMASTATE objects within a model would be:
!!>
!!   ! This first section of code is done during the parent model's initialzation,
!!   ! and there should be a unique CARMA object created for each thread of
!!   ! execution.
!!
!!   ! Create the CARMA object.   
!!   call CARMA_Create(carma, ...)
!!
!!   ! Define the microphysical components. 
!!   call CARMAGROUP_Create(carma, ...)      ! One or more calls
!!
!!   call CARMAELEMENT_Create(carma, ...)  ! One or more calls
!!
!!   call CARMASOLUTE_Create(carma, ...)    ! Zero or more calls
!!
!!   call CARMAGAS_Create(carma, ...)          ! Zero or more calls
!!
!!   ! Define the relationships for the microphysical processes. 
!!   call CARMA_AddCoagulation(carma, ...)    ! Zero or more calls
!!   call CARMA_AddGrowth(carma, ...)         ! Zero or more calls
!!   call CARMA_AddNucleation(carma, ...)     ! Zero or more calls
!!
!!   ! Initialize things that are state and timestep independent.
!!   call CARMA_Initialize(carma, ...)
!!
!!   ...
!!   
!!   ! This section of code is within the parent model's timing loop.
!!   !
!!   ! NOTE: If using OPEN/MP, then each thread will execute one of
!!   ! of these loops per column of data. To avoid having to destroy
!!   ! the CARMASTATE object, a pool of CARMASTATE objects could be
!!   ! created so that there is one per thread and then the
!!   ! CARMA_Destroy() could be called after all columns have been
!!   ! processed.
!!
!!   ! Initialize CARMA for this model state and timestep.
!!   call CARMASTATE_Create(cstate, carma, ...)
!!
!!   ! Set the model state for each bin and gas.
!!   call CARMASTATE_SetBin(cstate, ...)          ! One call for each bin
!!   call CARMASTATE_SetGas(cstate, ...)          ! One call for each gas
!!
!!   ! Calculate the new state
!!   call CARMASTATE_Step(cstate, ...)
!!
!!   ! Get the results to return back to the parent model.
!!   call CARMASTATE_GetBin(cstate, ...)      ! One call for each Bin
!!   call CARMASTATE_GetGas(cstate, ...)      ! One call for each gas
!!   call CARMASTATE_GetState(cstate, ...)    ! Zero or one calls
!!
!!   ! (optional) Deallocate arrays that are not needed beyond this timestep.
!!   call CARMASTATE_Destroy(cstate)
!!
!!   ...
!!
!!   ! This section of code is done during the parent model's cleanup.
!!
!!   ! Deallocate all arrays.
!!   call CARMA_Destroy(carma)
!!<
!!
!!  @version Feb-2009 
!!  @author  Chuck Bardeen, Pete Colarco, Jamie Smith 
!
! NOTE: Documentation for this code can be generated automatically using f90doc,
! which is freely available from:
!   http://erikdemaine.org/software/f90doc/
! Comment lines with double comment characters are processed by f90doc, and there are
! some special characters added to the comments to control the documentation process.
! In addition to the special characters mentioned in the f990doc documentation, html
! formatting tags (e.g. <i></i>, <sup></sup>, ...) can also be added to the f90doc
! comments.
module carma_mod

  ! This module maps the parents models constants into the constants need by CARMA. NOTE: CARMA
  ! constants are in CGS units, while the parent models are typically in MKS units.
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod

  ! CARMA explicitly declares all variables. 
  implicit none

  ! All CARMA variables and procedures are private except those explicitly declared to be public.
  private

  ! Declare the public methods.
  public CARMA_AddCoagulation
  public CARMA_AddGrowth
  public CARMA_AddNucleation
  public CARMA_Create
  public CARMA_Destroy
  public CARMA_Get
  public CARMA_Initialize

contains

  ! These are the methods that provide the interface between the parent model and the CARMA
  ! microphysical model. There are many other methods that are not in this file that are
  ! used to implement the microphysical calculations needed by the CARMA model. These other
  ! methods are in effect private methods of the CARMA module, but are in individual files
  ! since that is the way that CARMA has traditionally been structured and where users may
  ! want to extend or replace code to affect the microphysics.

  !! Creates the CARMA object and allocates arrays to store configuration information
  !! that will follow from the CARMA_AddXXX() methods. When the CARMA object is no longer
  !! needed, the CARMA_Destroy() method should be used to clean up any allocations
  !! that have happened. If LUNOPRT is specified, then the logical unit should be open and
  !! ready for output. The caller is responsible for closing the LUNOPRT logical unit
  !! after the CARMA object has been destroyed.
  !!
  !!  @version Feb-2009 
  !!  @author  Chuck Bardeen 
  subroutine CARMA_Create(carma, NBIN, NELEM, NGROUP, NSOLUTE, NGAS, NWAVE, rc, LUNOPRT, wave)
    type(carma_type), intent(out)      :: carma     !! the carma object
    integer, intent(in)                :: NBIN      !! number of radius bins per group
    integer, intent(in)                :: NELEM     !! total number of elements
    integer, intent(in)                :: NGROUP    !! total number of groups
    integer, intent(in)                :: NSOLUTE   !! total number of solutes
    integer, intent(in)                :: NGAS      !! total number of gases
    integer, intent(in)                :: NWAVE     !! number of wavelengths
    integer, intent(out)               :: rc        !! return code, negative indicates failure
    integer, intent(in), optional      :: LUNOPRT   !! logical unit number for output
    real(kind=f), intent(in), optional :: wave(NWAVE)  !! wavelength (cm)
    
    ! Local Varaibles      
    integer                            :: ier
    
    ! Assume success.
    rc = RC_OK
    
    ! Save off the logic unit used for output if one was provided. If one was provided,
    ! then assume that CARMA can print output.
    if (present(LUNOPRT)) then 
      carma%LUNOPRT = LUNOPRT
      carma%do_print = .TRUE.
    end if
    
    ! Save the defintion of the number of comonents involved in the microphysics.
    carma%NGROUP  = NGROUP 
    carma%NELEM   = NELEM
    carma%NBIN    = NBIN
    carma%NGAS    = NGAS
    carma%NSOLUTE = NSOLUTE      
    carma%NWAVE   = NWAVE


    ! Allocate tables for the groups.
    allocate( &
      carma%group(NGROUP), &
      carma%icoag(NGROUP, NGROUP), &
      carma%inucgas(NGROUP), &
      stat=ier) 
    if(ier /= 0) then
      if (carma%do_print) write(carma%LUNOPRT, *) "CARMA_Create: ERROR allocating groups, NGROUP=", &
        carma%NGROUP, ", status=", ier
      rc = RC_ERROR
      return
    endif
    
    ! Initialize
    carma%icoag(:, :) = 0
    carma%inucgas(:) = 0
    
    
    ! Allocate tables for the elements.
    allocate( &
      carma%element(NELEM), &
      carma%igrowgas(NELEM), &
      carma%inuc2elem(NELEM, NELEM), &
      carma%inucproc(NELEM, NELEM), &
      carma%ievp2elem(NELEM), &
      carma%nnuc2elem(NELEM), &
      carma%nnucelem(NELEM), &
      carma%inucelem(NELEM,NELEM*NGROUP), &
      carma%if_nuc(NELEM,NELEM), &
      carma%rlh_nuc(NELEM, NELEM), &
      carma%icoagelem(NELEM, NGROUP), &
      carma%icoagelem_cm(NELEM, NGROUP), &
      stat=ier) 
    if(ier /= 0) then
      if (carma%do_print) write(carma%LUNOPRT, *) "CARMA_Create: ERROR allocating elements, NELEM=", &
        carma%NELEM, ", status=", ier
      rc = RC_ERROR
      return
    endif
    
    ! Initialize
    carma%igrowgas(:) = 0
    carma%inuc2elem(:,:) = 0
    carma%inucproc(:,:) = 0
    carma%ievp2elem(:) = 0
    carma%nnuc2elem(:) = 0
    carma%nnucelem(:) = 0
    carma%inucelem(:,:) = 0
    carma%if_nuc(:,:) = .FALSE.
    carma%rlh_nuc(:,:) = 0._f
    carma%icoagelem(:,:) = 0
    carma%icoagelem_cm(:,:) = 0
    
    
    ! Allocate talbes for the bins.
    allocate( &
      carma%inuc2bin(NBIN,NGROUP,NGROUP), &
      carma%ievp2bin(NBIN,NGROUP,NGROUP), &
      carma%nnucbin(NGROUP,NBIN,NGROUP), &
      carma%inucbin(NBIN*NGROUP,NGROUP,NBIN,NGROUP), &
      carma%diffmass(NBIN, NGROUP, NBIN, NGROUP), &
      carma%volx(NGROUP,NGROUP,NGROUP,NBIN,NBIN), &
      carma%ilow(NGROUP,NBIN,NBIN*NBIN), &
      carma%jlow(NGROUP,NBIN,NBIN*NBIN), &
      carma%iup(NGROUP,NBIN,NBIN*NBIN), &
      carma%jup(NGROUP,NBIN,NBIN*NBIN), &
      carma%npairl(NGROUP,NBIN), &
      carma%npairu(NGROUP,NBIN), &
      carma%iglow(NGROUP,NBIN,NBIN*NBIN), &
      carma%jglow(NGROUP,NBIN,NBIN*NBIN), &
      carma%igup(NGROUP,NBIN,NBIN*NBIN), &
      carma%jgup(NGROUP,NBIN,NBIN*NBIN), &
      carma%kbin(NGROUP,NGROUP,NGROUP,NBIN,NBIN), &
      carma%pratt(3,NBIN,NGROUP), &
      carma%prat(4,NBIN,NGROUP), &
      carma%pden1(NBIN,NGROUP), &
      carma%palr(4,NGROUP), &
      stat=ier) 
    if(ier /= 0) then
      if (carma%do_print) write(carma%LUNOPRT, *) "CARMA_Create: ERROR allocating bins, NBIN=", &
        carma%NBIN, ", status=", ier
      rc = RC_ERROR
      return
    endif
    
    ! Initialize
    carma%inuc2bin(:,:,:) = 0
    carma%ievp2bin(:,:,:) = 0
    carma%nnucbin(:,:,:) = 0
    carma%inucbin(:,:,:,:) = 0
    carma%diffmass(:, :, :, :) = 0._f
    carma%volx(:,:,:,:,:) = 0._f
    carma%ilow(:,:,:) = 0
    carma%jlow(:,:,:) = 0
    carma%iup(:,:,:) = 0
    carma%jup(:,:,:) = 0
    carma%npairl(:,:) = 0
    carma%npairu(:,:) = 0
    carma%iglow(:,:,:) = 0
    carma%jglow(:,:,:) = 0
    carma%igup(:,:,:) = 0
    carma%jgup(:,:,:) = 0
    carma%kbin(:,:,:,:,:) = 0._f
    carma%pratt(:,:,:) = 0._f
    carma%prat(:,:,:) = 0._f
    carma%pden1(:,:) = 0._f
    carma%palr(:,:) = 0._f
      

    ! Allocate tables for solutes, if any are needed.
    if (NSOLUTE > 0) then
      allocate( &
        carma%solute(NSOLUTE), &
        stat=ier) 
      if(ier /= 0) then
        if (carma%do_print) write(carma%LUNOPRT, *) "CARMA_Create: ERROR allocating solutes, NSOLUTE=", &
          carma%NSOLUTE, ", status=", ier 
        rc = RC_ERROR
        return
      endif
    end if
    
   
    ! Allocate tables for gases, if any are needed.
    if (NGAS > 0) then
      allocate( &
        carma%gas(NGAS), &
        stat=ier) 
      if(ier /= 0) then
        if (carma%do_print) write(carma%LUNOPRT, *) "CARMA_Create: ERROR allocating gases, NGAS=", &
          carma%NGAS, ", status=", ier
        rc = RC_ERROR
        return
      endif
    end if
    
    
    ! Allocate tables for optical properties, if any are needed.
    if (NWAVE > 0) then
      allocate( &
        carma%wave(NWAVE), &
        stat=ier) 
      if(ier /= 0) then
        if (carma%do_print) write(carma%LUNOPRT, *) "CARMA_Create: ERROR allocating wavelengths, NWAVE=", &
          carma%NWAVE, ", status=", ier
        rc = RC_ERROR
        return
      endif

      ! Initialize
      carma%wave(:) = wave(:)
    end if
    
    return
 end subroutine CARMA_Create

  !! Called after the CARMA object has been created and the microphysics description has been
  !! configured. The optional flags control which microphysical processes are enabled and all of
  !! them default to FALSE. For a microphysical process to be active it must have been both
  !! configured (using a CARMA_AddXXX() method) and enabled here.
  !!
  !! NOTE: After initialization, the structure of the particle size bins is determined, and
  !! the resulting r, dr, rmass and dm can be retrieved with the CARMA_GetGroup() method.
  !!
  !!  @version Feb-2009 
  !!  @author  Chuck Bardeen 
  subroutine CARMA_Initialize(carma, rc, do_cnst_rlh, do_coag, do_detrain, do_fixedinit, do_grow, do_incloud, do_explised, do_print_init, do_substep, &
      do_thermo, do_vdiff, do_vtran, vf_const, minsubsteps, maxsubsteps, maxretries, conmax)
    type(carma_type), intent(inout)     :: carma         !! the carma object
    integer, intent(out)                :: rc            !! return code, negative indicates failure
    logical, intent(in), optional       :: do_cnst_rlh   !! use constant values for latent heats (instead of varying with temperature)?
    logical, intent(in), optional       :: do_coag       !! do coagulation?
    logical, intent(in), optional       :: do_detrain    !! do detrainement?
    logical, intent(in), optional       :: do_fixedinit  !! do initialization from reference atm?
    logical, intent(in), optional       :: do_grow       !! do nucleation, growth and evaporation?
    logical, intent(in), optional       :: do_incloud    !! do incloud growth and coagulation?
    logical, intent(in), optional       :: do_explised   !! do sedimentation with substepping
    logical, intent(in), optional       :: do_substep    !! do substepping
    logical, intent(in), optional       :: do_print_init !! do prinit initializtion information
    logical, intent(in), optional       :: do_thermo     !! do thermodynamics
    logical, intent(in), optional       :: do_vdiff      !! do Brownian diffusion
    logical, intent(in), optional       :: do_vtran      !! do sedimentation
    real(kind=f), intent(in), optional  :: vf_const      !! if specified and non-zero, constant fall velocity for all particles [cm/s]
    integer, intent(in), optional       :: minsubsteps   !! minimum number of substeps, default = 1
    integer, intent(in), optional       :: maxsubsteps   !! maximum number of substeps, default = 1
    integer, intent(in), optional       :: maxretries    !! maximum number of substep retries, default = 5
    real(kind=f), intent(in), optional  :: conmax        !! minimum relative concentration to consider, default = 1e-1
    
    ! Assume success.
    rc = RC_OK

    ! Set default values for control flags.
    carma%do_cnst_rlh   = .FALSE.
    carma%do_coag       = .FALSE.
    carma%do_detrain    = .FALSE.
    carma%do_fixedinit  = .FALSE.
    carma%do_grow       = .FALSE.
    carma%do_incloud    = .FALSE.
    carma%do_explised   = .FALSE.
    carma%do_print_init = .FALSE.
    carma%do_substep    = .FALSE.
    carma%do_thermo     = .FALSE.
    carma%do_vdiff      = .FALSE.
    carma%do_vtran      = .FALSE.
    
    ! Store off any control flag values that have been supplied.
    if (present(do_cnst_rlh))   carma%do_cnst_rlh   = do_cnst_rlh
    if (present(do_coag))       carma%do_coag       = do_coag
    if (present(do_detrain))    carma%do_detrain    = do_detrain
    if (present(do_fixedinit))  carma%do_fixedinit  = do_fixedinit
    if (present(do_grow))       carma%do_grow       = do_grow
    if (present(do_incloud))     carma%do_incloud   = do_incloud
    if (present(do_explised))   carma%do_explised   = do_explised
    if (present(do_print_init)) carma%do_print_init = (do_print_init .and. carma%do_print)
    if (present(do_substep))    carma%do_substep    = do_substep
    if (present(do_thermo))     carma%do_thermo     = do_thermo
    if (present(do_vdiff))      carma%do_vdiff      = do_vdiff
    if (present(do_vtran))      carma%do_vtran      = do_vtran 
    
    ! Setup the bin structure.
    call setupbins(carma, rc)
    if (rc < 0) return
    
    ! Substepping
    carma%minsubsteps = 1         ! minimum number of substeps
    carma%maxsubsteps = 1         ! maximum number of substeps
    carma%maxretries  = 1         ! maximum number of retries
    carma%conmax      = 1.e-1_f
    
    if (present(minsubsteps)) carma%minsubsteps = minsubsteps
    if (present(maxsubsteps)) carma%maxsubsteps = maxsubsteps
    if (present(maxretries))  carma%maxretries  = maxretries
    if (present(conmax))      carma%conmax      = conmax

    carma%do_step = .TRUE.
    
    ! Calculate the Optical Properties
    call CARMA_InitializeOptics(carma, rc)
    if (rc < 0) return

    ! If any of the processes have initialization that can be done without the state
    ! information, then perform that now. This will mostly be checking the configuration
    ! and setting up any tables based upon the configuration.
    if (carma%do_vtran .or. carma%do_coag)  then
      call CARMA_InitializeVertical(carma, rc, vf_const)
      if (rc < 0) return
    end if
    
    if (carma%do_coag) then
      call setupcoag(carma, rc)
      if (rc < 0) return
    end if
      
    if (carma%do_grow) then
      call CARMA_InitializeGrowth(carma, rc)
      if (rc < 0) return
    end if
      
    if (carma%do_thermo) then
      call CARMA_InitializeThermo(carma, rc)
       if (rc < 0) return
    end if
   
    return
  end subroutine CARMA_Initialize


  subroutine CARMA_InitializeGrowth(carma, rc)
    type(carma_type), intent(inout)    :: carma
    integer, intent(out)             :: rc
      
    ! Local Variables
    integer                            :: i
    logical                            :: bad_grid
    integer                            :: igroup   ! group index
    integer                            :: igas     ! gas index
    integer                            :: isol     ! solute index
    integer                            :: ielem    ! element index
    integer                            :: ibin     ! bin index
    integer                            :: igfrom
    integer                            :: igto
    integer                            :: ibto
    integer                            :: ieto
    integer                            :: ifrom
    integer                            :: iefrom
    integer                            :: jefrom
    integer                            :: ip
    integer                            :: jcore
    integer                            :: iecore
    integer                            :: im
    integer                            :: jnucelem
    integer                            :: inuc2
    integer                            :: neto
    integer                            :: jfrom
    integer                            :: j
    integer                            :: nnucb

    ! Define formats
    1 format(a,':  ',12i6)
    2 format(/,a,':  ',i6)
    3 format(a,a)
    4 format(a,':  ',1pe12.3)
    5 format(/,'Particle nucleation mapping arrays (setupnuc):')
    7 format(/,'Warning: nucleation cannot occur from group',i3, &
               '   bin',i3,'   into group',i3,'   (<inuc2bin> is zero)')
    
    
    ! Assume success.
    rc = RC_OK

    ! Compute radius-dependent terms used in PPM advection scheme
    do igroup = 1, carma%NGROUP
      do i = 2,carma%NBIN-1
        carma%pratt(1,i,igroup) = carma%group(igroup)%dm(i) / &
              ( carma%group(igroup)%dm(i-1) + carma%group(igroup)%dm(i) + carma%group(igroup)%dm(i+1) )
        carma%pratt(2,i,igroup) = ( 2._f*carma%group(igroup)%dm(i-1) + carma%group(igroup)%dm(i) ) / &
              ( carma%group(igroup)%dm(i+1) + carma%group(igroup)%dm(i) )
        carma%pratt(3,i,igroup) = ( 2._f*carma%group(igroup)%dm(i+1) + carma%group(igroup)%dm(i) ) / &
              ( carma%group(igroup)%dm(i-1) + carma%group(igroup)%dm(i) )
      enddo

      do i = 2,carma%NBIN-2
        carma%prat(1,i,igroup) = carma%group(igroup)%dm(i) / &
                ( carma%group(igroup)%dm(i) + carma%group(igroup)%dm(i+1) )
        carma%prat(2,i,igroup) = 2._f * carma%group(igroup)%dm(i+1) * carma%group(igroup)%dm(i) / &
               ( carma%group(igroup)%dm(i) + carma%group(igroup)%dm(i+1) )
        carma%prat(3,i,igroup) = ( carma%group(igroup)%dm(i-1) + carma%group(igroup)%dm(i) ) / &
               ( 2._f*carma%group(igroup)%dm(i) + carma%group(igroup)%dm(i+1) )
        carma%prat(4,i,igroup) = ( carma%group(igroup)%dm(i+2) + carma%group(igroup)%dm(i+1) ) / &
               ( 2._f*carma%group(igroup)%dm(i+1) + carma%group(igroup)%dm(i) )
        carma%pden1(i,igroup) = carma%group(igroup)%dm(i-1) + carma%group(igroup)%dm(i) + &
               carma%group(igroup)%dm(i+1) + carma%group(igroup)%dm(i+2)
      enddo

      if( carma%NBIN .gt. 1 )then
        carma%palr(1,igroup) = &
             (carma%group(igroup)%rmassup(1)-carma%group(igroup)%rmass(1)) / &
             (carma%group(igroup)%rmass(2)-carma%group(igroup)%rmass(1))
        carma%palr(2,igroup) = &
             (carma%group(igroup)%rmassup(1)/carma%group(igroup)%rmrat-carma%group(igroup)%rmass(1)) / &
             (carma%group(igroup)%rmass(2)-carma%group(igroup)%rmass(1))
        carma%palr(3,igroup) = &
             (carma%group(igroup)%rmassup(carma%NBIN-1)-carma%group(igroup)%rmass(carma%NBIN-1)) &
             / (carma%group(igroup)%rmass(carma%NBIN)-carma%group(igroup)%rmass(carma%NBIN-1))
        carma%palr(4,igroup) = &
             (carma%group(igroup)%rmassup(carma%NBIN)-carma%group(igroup)%rmass(carma%NBIN-1)) &
             / (carma%group(igroup)%rmass(carma%NBIN)-carma%group(igroup)%rmass(carma%NBIN-1))
      endif
    end do
    
    
    ! Check the nucleation mapping.
    !
    ! NOTE: This code was moved from setupnuc, because it is not dependent on the model's
    ! state. A small part of setupnuc which deals with scrit is state specific, and that was
    ! left in setupnuc.

    ! Bin mapping for nucleation : nucleation would transfer mass from particles
    ! in <ifrom,igfrom> into target bin <inuc2bin(ifrom,igfrom,igto)> in group
    ! <igto>.  The target bin is the smallest bin in the target size grid with
    ! mass exceeding that of nucleated particle.
    do igfrom = 1,carma%NGROUP    ! nucleation source group
      do igto = 1,carma%NGROUP        ! nucleation target group
        do ifrom = 1,carma%NBIN   ! nucleation source bin
  
          carma%inuc2bin(ifrom,igfrom,igto) = 0
  
          do ibto = carma%NBIN,1,-1        ! nucleation target bin
  
            if( carma%group(igto)%rmass(ibto) .ge. carma%group(igfrom)%rmass(ifrom) )then
              carma%inuc2bin(ifrom,igfrom,igto) = ibto
            endif
          enddo
        enddo
      enddo
    enddo

    ! Mappings for nucleation sources: 
    !
    !  <nnucelem(ielem)> is the number of particle elements that nucleate to
    !   particle element <ielem>.
    !
    !  <inuc2elem(jefrom,ielem)> are the particle elements that
    !   nucleate to particle element <ielem>, where 
    !   jefrom = 1,nnucelem(ielem).
    !
    !  <if_nuc(iefrom,ieto)> is true if nucleation transfers mass from element
    !   <iefrom> to element <ieto>.
    !
    !  <nnucbin(igfrom,ibin,igroup)> is the number of particle bins that nucleate
    !   to particles in bin <ibin,igroup> from group <igfrom>.
    !
    !  <inucbin(jfrom,igfrom,ibin,igto)> are the particle bins 
    !   that nucleate to particles in bin <ibin,igto>, where
    !   jfrom = 1,nnucbin(igfrom,ibin,igto).
    !
    !
    ! First, calculate <nnucelem(ielem)> and <if_nuc(iefrom,ieto)>
    ! based on <inucelem(jefrom,ielem)>
    do iefrom = 1,carma%NELEM
      do ieto = 1,carma%NELEM
        carma%if_nuc(iefrom,ieto) = .false.
      enddo
    enddo
    
    do ielem = 1,carma%NELEM
      carma%nnuc2elem(ielem) = 0
      
      do jefrom = 1,carma%NGROUP
        if( carma%inuc2elem(jefrom,ielem) .ne. 0 ) then
          carma%nnuc2elem(ielem) = carma%nnuc2elem(ielem) + 1
          carma%if_nuc(ielem,carma%inuc2elem(jefrom,ielem)) = .true.

      
          ! Also check for cases where neither the source or destinaton don't have cores (e.g.
          ! melting ice to water drops).
          if ((carma%group(carma%element(ielem)%igroup)%ncore .eq. 0) .and. &
              (carma%group(carma%element(carma%inuc2elem(jefrom,ielem))%igroup)%ncore .eq. 0)) then
      
            ! For particle concentration target elements, only count source elements
            ! that are also particle concentrations.
            carma%nnucelem(carma%inuc2elem(jefrom,ielem)) = carma%nnucelem(carma%inuc2elem(jefrom,ielem)) + 1
            carma%inucelem(carma%nnucelem(carma%inuc2elem(jefrom,ielem)),carma%inuc2elem(jefrom,ielem)) = ielem
          end if
        endif
      enddo
    enddo
    
    ! Next, enumerate and count elements that nucleate to cores.
    do igroup = 1,carma%NGROUP

      ip = carma%group(igroup)%ienconc    ! target particle number concentration element

      do jcore = 1,carma%group(igroup)%ncore

        iecore = carma%group(igroup)%icorelem(jcore)    ! target core element 
!        carma%nnucelem(iecore) = 0

        do iefrom = 1,carma%NELEM

          if( carma%if_nuc(iefrom,iecore) ) then
            carma%nnucelem(iecore) = carma%nnucelem(iecore) + 1
            carma%inucelem(carma%nnucelem(iecore),iecore) = iefrom
          endif
        enddo      ! iefrom=1,NELEM
      enddo        ! jcore=1,ncore
    enddo          ! igroup=1,NGROUP
    

    ! Now enumerate and count elements nucleating to particle concentration
    ! (itype=I_INVOLATILE and itype=I_VOLATILE) and core second moment
    ! (itype=I_COREMASS).  Elements with itype = I_VOLATILE are special because all
    ! nucleation sources for core elements in same group are also sources
    ! for the itype = I_VOLATILE element.
    do igroup = 1,carma%NGROUP
    
      ip = carma%group(igroup)%ienconc    ! target particle number concentration element
      im = carma%group(igroup)%imomelem   ! target core second moment element

!      carma%nnucelem(ip) = 0
!      if( im .ne. 0 )then
!        carma%nnucelem(im) = 0
!      endif

      do jcore = 1,carma%group(igroup)%ncore

        iecore = carma%group(igroup)%icorelem(jcore)       ! target core mass element

        do jnucelem = 1,carma%nnucelem(iecore)  ! elements nucleating to cores

          iefrom = carma%inucelem(jnucelem,iecore)  ! source
          
          ! For particle concentration target elements, only count source elements
          ! that are also particle concentrations.
          carma%nnucelem(ip) = carma%nnucelem(ip) + 1
          carma%inucelem(carma%nnucelem(ip),ip) = carma%group(carma%element(iefrom)%igroup)%ienconc

          if( im .ne. 0 )then
            carma%nnucelem(im) = carma%nnucelem(im) + 1
            carma%inucelem(carma%nnucelem(im),im) = iefrom
          endif
        enddo
      enddo       ! jcore=1,ncore
    enddo         ! igroup=1,NGROUP


    ! Now enumerate and count nucleating bins.
    do igroup = 1,carma%NGROUP    ! target group
      do ibin = 1,carma%NBIN    ! target bin
        do igfrom = 1,carma%NGROUP    ! source group

          carma%nnucbin(igfrom,ibin,igroup) = 0

          do ifrom = 1,carma%NBIN   ! source bin

            if( carma%inuc2bin(ifrom,igfrom,igroup) .eq. ibin ) then
              carma%nnucbin(igfrom,ibin,igroup) = carma%nnucbin(igfrom,ibin,igroup) + 1
              carma%inucbin(carma%nnucbin(igfrom,ibin,igroup),igfrom,ibin,igroup) = ifrom
            endif
          enddo
        enddo   ! igfrom=1,NGROUP
      enddo   ! ibin=1,NBIN=1,NGROUP
    enddo   ! igroup=1,NGROUP

    if (carma%do_print_init) then
      
      !  Report nucleation mapping arrays (should be 'write' stmts, of course)

      write(carma%LUNOPRT,*) ' '
      write(carma%LUNOPRT,*) 'Nucleation mapping arrays (setupnuc):'
      write(carma%LUNOPRT,*) ' '
      write(carma%LUNOPRT,*) 'Elements mapping:'
      
      do ielem = 1,carma%NELEM
        write(carma%LUNOPRT,*) 'ielem,nnucelem=',ielem,carma%nnucelem(ielem)
       
        if(carma%nnucelem(ielem) .gt. 0) then
          do jfrom = 1,carma%nnucelem(ielem)
            write(carma%LUNOPRT,*) 'jfrom,inucelem=  ',jfrom,carma%inucelem(jfrom,ielem)
          enddo
        endif
      enddo
      
      write(carma%LUNOPRT,*) ' '
      write(carma%LUNOPRT,*) 'Bin mapping:'
      
      do igfrom = 1,carma%NGROUP
        do igroup = 1,carma%NGROUP
          write(carma%LUNOPRT,*) ' '
          write(carma%LUNOPRT,*) 'Groups (from, to) = ', igfrom, igroup
          
          do ibin = 1,carma%NBIN
            nnucb = carma%nnucbin(igfrom,ibin,igroup)
            if(nnucb .eq. 0) write(carma%LUNOPRT,*) '  None for bin ',ibin
            if(nnucb .gt. 0) then
              write(carma%LUNOPRT,*) '  ibin,nnucbin=',ibin,nnucb
              write(carma%LUNOPRT,*) '   inucbin=',(carma%inucbin(j,igfrom,ibin,igroup),j=1,nnucb)
            endif
          enddo
        enddo
      enddo
    endif


    ! Check that values are valid.
    do ielem = 1, carma%NELEM

      if( carma%element(ielem)%isolute .gt. carma%NSOLUTE )then
        if (carma%do_print) write(carma%LUNOPRT,*) 'CARMA_InitializeGrowth::ERROR - component of isolute > NSOLUTE'
        rc = RC_ERROR
        return
      endif

      if( carma%ievp2elem(ielem) .gt. carma%NELEM )then
        if (carma%do_print) write(carma%LUNOPRT,*) 'CARMA_InitializeGrowth::ERROR - component of ievp2elem > NELEM'
        rc = RC_ERROR
        return
      endif

      ! Check that <isolute> is consistent with <ievp2elem>.
      if( carma%ievp2elem(ielem) .ne. 0 .and. carma%element(ielem)%itype .eq. I_COREMASS )then
        if( carma%element(ielem)%isolute .ne. carma%element(carma%ievp2elem(ielem))%isolute)then
          if (carma%do_print) write(carma%LUNOPRT,*) 'CARMA_InitializeGrowth::ERROR - isolute and ievp2elem are inconsistent'
          rc = RC_ERROR
          return
        endif
      endif

      ! Check that <isolute> is consistent with <inucgas>.
      igas = carma%inucgas( carma%element(ielem)%igroup )
      if( igas .ne. 0 )then
        if( carma%element(ielem)%itype .eq. I_COREMASS .and. carma%element(ielem)%isolute .eq. 0 )then
          if (carma%do_print) write(carma%LUNOPRT,*) 'CARMA_InitializeGrowth::ERROR - inucgas ne 0 but isolute eq 0'
          rc = RC_ERROR
          return
        endif
      endif
    enddo

    do ielem = 1, carma%NELEM
      if( carma%nnuc2elem(ielem) .gt. 0 ) then
        do inuc2 = 1, carma%nnuc2elem(ielem)
          if( carma%inuc2elem(inuc2,ielem) .gt. carma%NELEM )then
            if (carma%do_print) write(carma%LUNOPRT,*) 'CARMA_InitializeGrowth::ERROR - component of inuc2elem > NELEM'
            rc = RC_ERROR
            return
          endif
        enddo
      endif
    enddo

    ! Particle grids are incompatible if there is no target bin with enough
    ! mass to accomodate nucleated particle.
    bad_grid = .false.

    do iefrom = 1,carma%NELEM   ! source element

      igfrom = carma%element(iefrom)%igroup
      neto   = carma%nnuc2elem(iefrom)

      if( neto .gt. 0 )then

        do inuc2 = 1,neto
          ieto = carma%inuc2elem(inuc2,iefrom)
          igto = carma%element(ieto)%igroup

          do ifrom = 1,carma%NBIN   ! source bin
            if( carma%inuc2bin(ifrom,igfrom,igto) .eq. 0 )then
              if (carma%do_print) write(carma%LUNOPRT,7) igfrom,ifrom,igto
              bad_grid = .true.
            endif
          enddo
        enddo
      endif
    enddo

    if( bad_grid )then
      if (carma%do_print) write(carma%LUNOPRT,*) 'CARMA_InitializeGrowth::ERROR - incompatible grids for nucleation'
      rc = RC_ERROR
      return
    endif
      
    if (carma%do_print_init) then
    
      ! Report some initialization values!
      write(carma%LUNOPRT,5)
      write(carma%LUNOPRT,1) 'inucgas  ',(carma%inucgas(i),i=1,carma%NGROUP)
      write(carma%LUNOPRT,1) 'inuc2elem',(carma%inuc2elem(1,i),i=1,carma%NELEM)
      write(carma%LUNOPRT,1) 'ievp2elem',(carma%ievp2elem(i),i=1,carma%NELEM)
      write(carma%LUNOPRT,1) 'isolute ',(carma%element(i)%isolute,i=1,carma%NELEM)
    
      do isol = 1,carma%NSOLUTE
        write(carma%LUNOPRT,2) 'solute number   ',isol
        write(carma%LUNOPRT,3) 'solute name:    ',carma%solute(isol)%name
        write(carma%LUNOPRT,4) 'molecular weight',carma%solute(isol)%wtmol
        write(carma%LUNOPRT,4) 'mass density    ',carma%solute(isol)%rho    
      enddo
    endif

    return
  end subroutine CARMA_InitializeGrowth         

  !! Calculate the optical properties for each particle bin at each of
  !! the specified wavelengths. The optical properties include the
  !! extinction efficiency, the single scattering albedo and the
  !! asymmetry factor.
  !!
  !! NOTE: For these calculations, the particles are assumed to be spheres and
  !! Mie code is used to calculate the optical properties.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  subroutine CARMA_InitializeOptics(carma, rc)
    type(carma_type), intent(inout)    :: carma
    integer, intent(out)               :: rc
    
    integer                            :: igroup      !! group index
    integer                            :: iwave       !! wavelength index
    integer                            :: ibin        !! bin index
    real(kind=f)                       :: theta(IT)
    real(kind=f)                       :: Qext
    real(kind=f)                       :: Qsca
    real(kind=f)                       :: Qbs
    real(kind=f)                       :: ctbrqs
    real(kind=f)                       :: wvno
    real(kind=f)                       :: rfr
    real(kind=f)                       :: rfi
    
    ! Assume success.
    rc = RC_OK    
    
    ! We only care about the forward direction.
    theta(:) = 0.0_f
    
    ! Were any wavelengths specified?
    do iwave = 1, carma%NWAVE
    
      ! Calcualte the wave number.
      wvno = 2._f * PI / (carma%wave(iwave))
   
      do igroup = 1, carma%NGROUP
     
        ! Should we calculate mie properties for this group?
        if (carma%group(igroup)%do_mie) then 
       
          rfr = real(carma%group(igroup)%refidx(iwave))
          rfi = imag(carma%group(igroup)%refidx(iwave))
        
          do ibin = 1, carma%NBIN
          
            ! Assume the particle is homogeneous (no core).
            call miess(carma, &
                       carma%group(igroup)%r(ibin), &
                       rfr, &
                       rfi, &
                       theta, &
                       1, &
                       Qext, &
                       Qsca, &
                       Qbs,&
                       ctbrqs, &
                       0.0_f, &
                       rfr, &
                       rfi, &
                       wvno, &
                       rc)
            if (rc < RC_OK) return
  
            carma%group(igroup)%qext(iwave, ibin) = Qext
            carma%group(igroup)%ssa(iwave, ibin)  = Qsca / Qext
            carma%group(igroup)%asym(iwave, ibin) = ctbrqs / Qsca
          end do
        end if
      end do
    end do   
    
    return
  end subroutine CARMA_InitializeOptics         

  !! Perform initialization of variables related to thermodynamical calculations that
  !! are not dependent on the model state.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  subroutine CARMA_InitializeThermo(carma, rc)
    type(carma_type), intent(inout)    :: carma
    integer, intent(out)               :: rc
        
    ! Assume success.
    rc = RC_OK

    return
  end subroutine CARMA_InitializeThermo         

  !! Perform initialization of variables related to vertical transport that are not dependent
  !! on the model state.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  subroutine CARMA_InitializeVertical(carma, rc, vf_const)
    type(carma_type), intent(inout)    :: carma
    integer, intent(out)               :: rc
    real(kind=f), intent(in), optional :: vf_const
        
    ! Assume success.
    rc = RC_OK

    ! Was a constant vertical velocity specified?
    carma%ifall = 1
    carma%vf_const = 0._f
    
    if (present(vf_const)) then
      if (vf_const /= 0._f) then
        carma%ifall = 0
        carma%vf_const = vf_const
      end if
    end if
    
    ! Specify the boundary conditions for vertical transport.
    carma%itbnd_pc  = I_FIXED_CONC
    carma%ibbnd_pc  = I_FIXED_CONC
    
    return
  end subroutine CARMA_InitializeVertical         

  !! The routine should be called when the carma object is no longer needed. It deallocates
  !! any memory allocations made by CARMA (during CARMA_Create()), and failure to call this
  !!routine could result in memory leaks.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  !!
  !! @see CARMA_Create
  subroutine CARMA_Destroy(carma, rc)
    use carmaelement_mod
    use carmagas_mod
    use carmagroup_mod
    use carmasolute_mod

    type(carma_type), intent(inout)    :: carma
    integer, intent(out)               :: rc
    
    ! Local variables
    integer   :: ier
    integer   :: igroup
    integer   :: ielem
    integer   :: isolute
    integer   :: igas
    
    ! Assume success.
    rc = RC_OK
    
    ! If allocated, deallocate all the variables that were allocated in the Create() method.
    if (allocated(carma%group)) then
      do igroup = 1, carma%NGROUP
        call CARMAGROUP_Destroy(carma, igroup, rc)
        if (rc < 0) return
      end do
      
      deallocate( &
        carma%group, &
        carma%icoag, &
        carma%inucgas, &
        stat=ier) 
      if(ier /= 0) then
        if (carma%do_print) write(carma%LUNOPRT, *) "CARMA_Destroy: ERROR deallocating groups, status=", ier
        rc = RC_ERROR
      endif
    endif
    
    if (allocated(carma%element)) then
      do ielem = 1, carma%NELEM
        call CARMAELEMENT_Destroy(carma, ielem, rc)
        if (rc < RC_OK) return
      end do
      
      deallocate( &
        carma%element, &
        carma%igrowgas, &
        carma%inuc2elem, &
        carma%inucproc, &
        carma%ievp2elem, &
        carma%nnuc2elem, &
        carma%nnucelem, &
        carma%inucelem, &
        carma%if_nuc, &
        carma%rlh_nuc, &
        carma%icoagelem, &
        carma%icoagelem_cm, &
        stat=ier) 
      if(ier /= 0) then
        if (carma%do_print) write(carma%LUNOPRT, *) "CARMA_Destroy: ERROR deallocating elements, status=", ier
        rc = RC_ERROR
      endif
    endif
    
    if (allocated(carma%inuc2bin)) then
      deallocate( &
        carma%inuc2bin, &
        carma%ievp2bin, &
        carma%nnucbin, &
        carma%inucbin, &
        carma%diffmass, &
        carma%volx, &
        carma%ilow, &
        carma%jlow, &
        carma%iup, &
        carma%jup, &
        carma%npairl, &
        carma%npairu, &
        carma%iglow, &
        carma%jglow, &
        carma%igup, &
        carma%jgup, &
        carma%kbin, &
        stat=ier) 
      if(ier /= 0) then
        if (carma%do_print) write(carma%LUNOPRT, *) "CARMA_Destroy: ERROR deallocating bins, status=", ier
        rc = RC_ERROR
      endif
    endif
    
    if (carma%NSOLUTE > 0) then
      do isolute = 1, carma%NSOLUTE
        call CARMASOLUTE_Destroy(carma, isolute, rc)
        if (rc < RC_OK) return
      end do
      
      if (allocated(carma%solute)) then
        deallocate( &
          carma%solute, &
          stat=ier) 
        if(ier /= 0) then
          if (carma%do_print) write(carma%LUNOPRT, *) "CARMA_Destroy: ERROR deallocating solutes, status=", ier 
          rc = RC_ERROR
        endif
      endif
    end if
     
    if (carma%NGAS > 0) then
      do igas = 1, carma%NGAS
        call CARMAGAS_Destroy(carma, igas, rc)
        if (rc < RC_OK) return
      end do
      
      if (allocated(carma%gas)) then
        deallocate( &
          carma%gas, &
          stat=ier) 
        if(ier /= 0) then
          if (carma%do_print) write(carma%LUNOPRT, *) "CARMA_Destroy: ERROR deallocating gases, status=", ier
          rc = RC_ERROR
        endif
      endif
    end if     
    
    if (carma%NWAVE > 0) then
      if (allocated(carma%wave)) then
        deallocate( &
          carma%wave, &
          stat=ier) 
        if(ier /= 0) then
          if (carma%do_print) write(carma%LUNOPRT, *) "CARMA_Destroy: ERROR deallocating wavelengths, status=", ier
          rc = RC_ERROR
          return
        endif
      endif
    endif
    
    return
  end subroutine CARMA_Destroy

  ! Configuration    
        
  !! Add a coagulation process between two groups (<i>igroup1</i> and <i>igroup2</i>), with the resulting
  !! particle being in the destination group (<i>igroup3</i>). If <i>ck0</i> is specifed, then a constant
  !! coagulation kernel will be used.
  subroutine CARMA_AddCoagulation(carma, igroup1, igroup2, igroup3, icollec, rc, ck0, grav_e_coll0)
    type(carma_type), intent(inout)    :: carma         !! the carma object
    integer, intent(in)                :: igroup1       !! first source group
    integer, intent(in)                :: igroup2       !! second source group
    integer, intent(in)                :: igroup3       !! destination group
    integer, intent(in)                :: icollec       !! collection technique [I_COLLEC_CONST | I_COLLEC_FUCHS | I_COLLEC_DATA] 
    integer, intent(out)               :: rc            !! return code, negative indicates failure
    real(kind=f), intent(in), optional :: ck0           !! if specified, forces a constant coagulation kernel
    real(kind=f), intent(in), optional :: grav_e_coll0  !! if <i>icollec</i> is I_COLLEC_CONST, the constant gravitational collection efficiency      
    
    ! Assume success.
    rc = RC_OK

    ! Make sure the groups exists.
    if (igroup1 > carma%NGROUP) then
      if (carma%do_print) write(carma%LUNOPRT, '(a,i3,a,i3,a)') "CARMA_AddCoagulation:: ERROR - The specifed group (", &
        igroup1, ") is larger than the number of groups (", carma%NGROUP, ")."
      rc = RC_ERROR
      return
    end if

    if (igroup2 > carma%NGROUP) then
      if (carma%do_print) write(carma%LUNOPRT, '(a,i3,a,i3,a)') "CARMA_AddCoagulation:: ERROR - The specifed group (", &
        igroup2, ") is larger than the number of groups (", carma%NGROUP, ")."
      rc = RC_ERROR
      return
    end if
    
    if (igroup3 > carma%NGROUP) then
      if (carma%do_print) write(carma%LUNOPRT, '(a,i3,a,i3,a)') "CARMA_AddCoagulation:: ERROR - The specifed group (", &
        igroup3, ") is larger than the number of groups (", carma%NGROUP, ")."
      rc = RC_ERROR
      return
    end if
    
    ! Indicate that the groups coagulate together.
    carma%icoag(igroup1, igroup2) = igroup3
    
    ! If ck0 was specified, then we use a fixed coagulation rate of ck0.
    if (present(ck0)) then
      carma%ck0 = ck0
      carma%icoagop = I_COAGOP_CONST
    else
      carma%icoagop = I_COAGOP_CALC
    end if
    
    ! What collection technique is specified.
    if (icollec > I_COLLEC_DATA) then
      if (carma%do_print) write(carma%LUNOPRT, '(a,i3,a)') "CARMA_AddCoagulation:: ERROR - The specifed collection method (", &
        icollec, ") is unknown."
      rc = RC_ERROR
      return
    end if
    
    if (icollec == I_COLLEC_CONST) then
      if (present(grav_e_coll0)) then
        carma%grav_e_coll0 = grav_e_coll0
      else
        if (carma%do_print) write(carma%LUNOPRT, *) "CARMA_AddCoagulation:: ERROR - A constant gravitational collection was requests, but grav_e_coll0 was not provided."
        rc = RC_ERROR
        return
      end if
    end if
    
    carma%icollec = icollec
    
    return
  end subroutine CARMA_AddCoagulation
    
  !! Add a growth process between the element (<i>ielem</i>) and gas (<i>igas</i>) specifed. The element
  !! and gas should have already been defined using <i>CARMA_AddElement()</i> and <i>CARMA_AddGas()</i>.
  !!
  !! NOTE: Each element can only have one volatile component.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  !!
  !! @see CARMA_AddElement
  !! @see CARMA_AddGas
  subroutine CARMA_AddGrowth(carma, ielem, igas, rc)
    type(carma_type), intent(inout)    :: carma    !! the carma object
    integer, intent(in)                :: ielem    !! the element index
    integer, intent(in)                :: igas     !! the gas index
    integer, intent(out)               :: rc       !! return code, negative indicates failure
    
    ! Assume success.
    rc = RC_OK
    
    ! Make sure the element exists.
    if (ielem > carma%NELEM) then
      if (carma%do_print) write(carma%LUNOPRT, *) "CARMA_AddGrowth:: ERROR - The specifed element (", &
        ielem, ") is larger than the number of elements (", carma%NELEM, ")."
      rc = RC_ERROR
      return
    end if

    ! Make sure there are enough gases allocated.
    if (igas > carma%NGAS) then
      if (carma%do_print) write(carma%LUNOPRT, *) "CARMA_AddGrowth:: ERROR - The specifed gas (", &
        igas, ") is larger than the number of gases (", carma%NGAS, ")."
      rc = RC_ERROR
      return
    end if
    
    ! If not already defined, indicate that the element can grow with the specified gas.
    if (carma%igrowgas(ielem) /= 0) then
      if (carma%do_print) write(carma%LUNOPRT, *) "CARMA_AddGrowth:: ERROR - The specifed element (", &
        ielem, ") already has gas (", carma%igrowgas(ielem), ") condensing on it."
      rc = RC_ERROR
      return
    else 
      carma%igrowgas(ielem) = igas
    end if
    
    return
  end subroutine CARMA_AddGrowth
    
  !! Add a nucleation process that nucleates one element (<i>elemfrom</i>) to another element (<i>elemto</i>)
  !! using the specified gas (<i>igas</i>). The elements and gas should have already been defined
  !! using <i>CARMA_AddElement()</i> and <i>CARMA_AddGas()</i>. The nucleation scheme is indicated by
  !! inucproc, and can be one of:
  !!
  !!   - <i>I_DROPACT</i>
  !!   - <i>I_AERFREEZE</i>
  !!   - <i>I_DROPFREEZE</i>
  !!   - <i>I_ICEMELT</i>
  !!   - <i>I_HETNUC</i>
  !!   - <i>I_GLFREEZE</i>
  !!   - <i>I_GLAERFREEZE</i>
  !! 
  !! Total evaporation transfers particle mass from the destination element back to the
  !! element indicated by <i>ievp2elem</i>. This relationship is not automatically generated,
  !! because multiple elements can nucleate to a particular element and therefore the 
  !! reverse mapping is not unique.
  !!
  !! NOTE: The gas used for nucleation must be the same for all nucleation defined from
  !! elements of the same group.
  !!
  !! @author Chuck Bardeen
  !! @version Feb-2009
  !! @see I_DROPACT
  !! @see I_AERFREEZE
  !! @see I_DROPFREEZE
  !! @see I_ICEMELT
  !! @see I_HETNUC
  !! @see I_GLFREEZE
  !! @see I_GLAERFREEZE
  !! @see CARMA_AddElement
  !! @see CARMA_AddGas
  subroutine CARMA_AddNucleation(carma, ielemfrom, ielemto, inucproc, &
      rlh_nuc, rc, igas, ievp2elem)
      
    use carmaelement_mod, only         : CARMAELEMENT_Get
    
    type(carma_type), intent(inout)    :: carma       !! the carma object
    integer, intent(in)                :: ielemfrom   !! the source element
    integer, intent(in)                :: ielemto     !! the destination element
    integer, intent(in)                :: inucproc    !! the nucleation process [I_DROPACT | I_AERFREEZE | I_ICEMELT | I_HETNUC | I_GLFREEZE | I_GLAERFREEZE]
    real(kind=f), intent(in)           :: rlh_nuc     !! the latent heat of nucleation [cm<sup>2</sup>/s<sup>2</sup>]
    integer, intent(out)               :: rc          !! return code, negative indicated failure
    integer, optional, intent(in)      :: igas        !! the gas
    integer, optional, intent(in)      :: ievp2elem   !! the element created upon evaporation
    
    integer                            :: igroup      !! group for source element
    
    ! Assume success.
    rc = RC_OK
    
    ! Make sure the elements exist.
    if (ielemfrom > carma%NELEM) then
      if (carma%do_print) write(carma%LUNOPRT, *) "CARMA_AddNucleation:: ERROR - The specifed element (", &
        ielemfrom, ") is larger than the number of elements (", carma%NELEM, ")."
      rc = RC_ERROR
      return
    end if

    if (ielemto > carma%NELEM) then
      if (carma%do_print) write(carma%LUNOPRT, *) "CARMA_AddNucleation:: ERROR - The specifed element (", &
        ielemto, ") is larger than the number of elements (", carma%NELEM, ")."
      rc = RC_ERROR
      return
    end if

    if (present(ievp2elem)) then
      if (ievp2elem > carma%NELEM) then
        if (carma%do_print) write(carma%LUNOPRT, *) "CARMA_AddNucleation:: ERROR - The specifed element (", &
          ievp2elem, ") is larger than the number of elements (", carma%NELEM, ")."
        rc = RC_ERROR
        return
      end if
    end if


    ! Make sure there are enough gases allocated.
    if (present(igas)) then
      if (igas > carma%NGAS) then
        if (carma%do_print) write(carma%LUNOPRT, *) "CARMA_AddGrowth:: ERROR - The specifed gas (", &
          igas, ") is larger than the number of gases (", carma%NGAS, ")."
        rc = RC_ERROR
        return
      end if
    end if
    
    
    ! Array <inucgas> maps a particle group to its associated gas for nucleation:
    ! Nucleation from group <igroup> is associated with gas <inucgas(igroup)>
    ! Set to zero if particles are not subject to nucleation.
    if (present(igas)) then
      call CARMAELEMENT_Get(carma, ielemfrom, rc, igroup=igroup)
      
      if (rc >= RC_OK) then
        carma%inucgas(igroup) = igas
      end if
    end if


    ! Nucleation transfers particle mass from element <ielem> to element
    ! <inuc2elem(i,ielem)>, where <i> ranges from 0 to the number of elements
    !  nucleating from <ielem>.
!    carma%nnucelem(ielemto) = carma%nnucelem(ielemto) + 1
!    carma%inucelem(carma%nnucelem(ielemto), ielemto) = ielemfrom
    carma%nnuc2elem(ielemfrom) = carma%nnuc2elem(ielemfrom) + 1
    carma%inuc2elem(carma%nnuc2elem(ielemfrom), ielemfrom) = ielemto
!    carma%if_nuc(ielemfrom,carma%inuc2elem(carma%nnuc2elem(ielemfrom), ielemfrom)) = .true.

    ! <inucproc(iefrom,ieto)> specifies what nucleation process nucleates
    ! particles from element <ielem> to element <ieto>:
    !   I_DROPACT:  Aerosol activation to droplets 
    !   I_AERFREEZE: Aerosol homogeneous freezing
    !   I_DROPFREEZE: Droplet homogeneous freezing
    !   I_GLFREEZE: Glassy Aerosol heteroogeneous freezing
    !   I_GLAERFREEZE: Glassy & Aerosol freezing
    carma%inucproc(ielemfrom, ielemto) = inucproc


    ! Total evaporation mapping: total evaporation transfers particle mass from
    ! element <ielem> to element <ievp2elem(ielem)>.
    !
    ! NOTE: This array is not automatically derived from <inuc2elem> because multiple
    ! elements can nucleate to a particular element (reverse mapping is not
    ! unique).
    if (present(ievp2elem)) carma%ievp2elem(ielemto) = ievp2elem


    ! <rlh_nuc(iefrom,ieto)> is the latent heat released by nucleation
    ! from element <iefrom> to element <ieto> [cm^2/s^2].
    carma%rlh_nuc(ielemfrom,ielemto) = rlh_nuc

    return
  end subroutine


  ! Query, Control and State I/O

  !! Gets the information about the carma object.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  !!
  !! @see CARMA_Create
  subroutine CARMA_Get(carma, rc, wave)
    type(carma_type), intent(in)        :: carma              !! the carma object
    integer, intent(out)                :: rc                 !! return code, negative indicates failure
    real(kind=f), optional, intent(out) :: wave(carma%NWAVE)  !! the wavelengths
    
    ! Assume success.
    rc = RC_OK
    
    if (present(wave)) wave(:) = carma%wave(:)
    
    return
  end subroutine CARMA_Get

end module
