!! The CARMAGROUP module contains configuration information about a CARMA partcile.
!!
!! NOTE: Because of the way Fortran handles pointers and allocations, it is much
!! simpiler to have these methods directly access the group array that is in the
!! CARMA object rather than having this as its own objects. Some compilers (like
!! IBM on AIX do not by default automatically deallocate automatically created
!! derived types that contain allocations. This can result in memory leaks that
!! are difficult to find.
!!
!! These calls are written like they are part of CARMA, but they are called
!! CARMAGROUP and kept by themselves in their own file to make it easier to keep
!! track of what is required when adding an attribute to a group.
!!
!!  @version July-2009 
!!  @author  Chuck Bardeen 
module carmagroup_mod

  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod

  ! CARMA explicitly declares all variables. 
  implicit none

  ! All CARMA variables and procedures are private except those explicitly declared to be public.
  private

  ! Declare the public methods.
  public CARMAGROUP_Create
  public CARMAGROUP_Destroy
  public CARMAGROUP_Get
  public CARMAGROUP_Print

contains

  subroutine CARMAGROUP_Create(carma, igroup, name, rmin, rmrat, ishape, eshape, is_ice, &
      rc, irhswell, irhswcomp, refidx, do_mie, do_wetdep, do_drydep, do_vtran, solfac, scavcoef, shortname, &
      cnsttype, maxbin, ifallrtn, is_cloud, rmassmin)
    type(carma_type), intent(inout)             :: carma               !! the carma object
    integer, intent(in)                         :: igroup              !! the group index
    character(*), intent(in)                    :: name                !! the group name, maximum of 255 characters
    real(kind=f), intent(in)                    :: rmin                !! the minimum radius, can be specified [cm]
    real(kind=f), intent(in)                    :: rmrat               !! the volume ratio between bins
    integer, intent(in)                         :: ishape              !! the type of the particle shape [I_SPHERE | I_HEXAGON | I_CYLINDER]
    real(kind=f), intent(in)                    :: eshape              !! the aspect ratio of the particle shape (length/diameter)
    logical, intent(in)                         :: is_ice              !! is this an ice particle?
    integer, intent(out)                        :: rc                  !! return code, negative indicates failure
    integer, optional, intent(in)               :: irhswell            !! the parameterization for particle swelling from relative humidity [I_FITZGERALD | I_GERBER]
    integer, optional, intent(in)               :: irhswcomp           !! the composition for particle swelling from relative humidity [I_FITZGERALD | I_GERBER]
    complex(kind=f), optional, intent(in)       :: refidx(carma%NWAVE) !! refractive index for the particle
    logical, optional, intent(in)               :: do_mie              !! do mie calculations?
    logical, optional, intent(in)               :: do_wetdep           !! do wet deposition for this particle?
    logical, optional, intent(in)               :: do_drydep           !! do dry deposition for this particle?
    logical, optional, intent(in)               :: do_vtran            !! do sedimentation for this particle?
    real(kind=f), intent(in), optional          :: solfac              !! the solubility factor for wet deposition
    real(kind=f), intent(in), optional          :: scavcoef            !! the scavenging coefficient for wet deposition
    character(*), optional, intent(in)          :: shortname           !! the group shortname, maximum of 6 characters
    integer, optional, intent(in)               :: cnsttype            !! constituent type in parent model [I_CNSTTYPE_PROGNOSTIC | I_CNSTTYPE_DIAGNOSTIC]
    integer, optional, intent(in)               :: maxbin              !! bin number of the last prognostic bin, the remaining bins are diagnostic
    integer, optional, intent(in)               :: ifallrtn            !! fall velocity routine [I_FALLRTN_STD | I_FALLRTN_STD_SHAPE | I_FALLRTN_HEYMSFIELD2010 | I_FALLRTN_ACKERMAN_DROP | I_FALLRTN_ACKERMAN_ICE]
    logical, optional, intent(in)               :: is_cloud            !! is this a cloud particle?
    real(kind=f), optional, intent(in)          :: rmassmin            !! the minimum mass, when used overrides rmin[g]

    ! Local variables
    integer                               :: ier
    
    ! Assume success.
    rc = RC_OK
    
    ! Make sure there are enough groups allocated.
    if (igroup > carma%NGROUP) then
      if (carma%do_print) write(carma%LUNOPRT, *) "CARMAGROUP_Add:: ERROR - The specifed group (", &
        igroup, ") is larger than the number of groups (", carma%NGROUP, ")."
      rc = RC_ERROR
      return
    end if
    
    allocate( &
      carma%group(igroup)%r(carma%NBIN), &
      carma%group(igroup)%rmass(carma%NBIN), &
      carma%group(igroup)%vol(carma%NBIN), &
      carma%group(igroup)%dr(carma%NBIN), &
      carma%group(igroup)%dm(carma%NBIN), &
      carma%group(igroup)%rmassup(carma%NBIN), &
      carma%group(igroup)%rup(carma%NBIN), &
      carma%group(igroup)%rlow(carma%NBIN), &
      carma%group(igroup)%icorelem(carma%NELEM), &
      stat=ier) 
    if(ier /= 0) then
        if (carma%do_print) write(carma%LUNOPRT, *) "CARMAGROUP_Add: ERROR allocating, status=", ier
      rc = RC_ERROR
      return
    end if
    
    ! Initialize
    carma%group(igroup)%r(:)        = 0._f
    carma%group(igroup)%rmass(:)    = 0._f
    carma%group(igroup)%vol(:)      = 0._f
    carma%group(igroup)%dr(:)       = 0._f
    carma%group(igroup)%dm(:)       = 0._f
    carma%group(igroup)%rmassup(:)  = 0._f
    carma%group(igroup)%rup(:)      = 0._f
    carma%group(igroup)%rlow(:)     = 0._f
    carma%group(igroup)%icorelem(:) = 0
    carma%group(igroup)%ifallrtn    = I_FALLRTN_STD
    carma%group(igroup)%is_cloud    = .false.
    

    ! Any optical properties?
    if (carma%NWAVE > 0) then
      allocate( &
        carma%group(igroup)%refidx(carma%NWAVE), &
        carma%group(igroup)%qext(carma%NWAVE,carma%NBIN), &
        carma%group(igroup)%ssa(carma%NWAVE,carma%NBIN), &
        carma%group(igroup)%asym(carma%NWAVE,carma%NBIN), &
        stat=ier) 
      if(ier /= 0) then
        if (carma%do_print) write(carma%LUNOPRT, *) "CARMAGROUP_Add: ERROR allocating, status=", ier
        rc = RC_ERROR
        return
      endif

      ! Initialize
      carma%group(igroup)%refidx(:) = (0._f, 0._f)
      carma%group(igroup)%qext(:,:) = 0._f
      carma%group(igroup)%ssa(:,:)  = 0._f
      carma%group(igroup)%asym(:,:) = 0._f
    end if
    

    ! Save off the settings.
    carma%group(igroup)%name        = name
    carma%group(igroup)%rmin        = rmin
    carma%group(igroup)%rmrat       = rmrat
    carma%group(igroup)%ishape      = ishape
    carma%group(igroup)%eshape      = eshape
    carma%group(igroup)%is_ice      = is_ice
    
    
    ! Defaults for optional parameters
    carma%group(igroup)%irhswell    = 0
    carma%group(igroup)%do_mie      = .false.
    carma%group(igroup)%do_wetdep   = .false.
    carma%group(igroup)%do_drydep   = .false.
    carma%group(igroup)%grp_do_vtran  = .true.
    carma%group(igroup)%solfac      = 0.3_f
    carma%group(igroup)%scavcoef    = 0.1_f
    carma%group(igroup)%shortname   = ""
    carma%group(igroup)%cnsttype    = I_CNSTTYPE_PROGNOSTIC
    carma%group(igroup)%maxbin      = carma%NBIN
    carma%group(igroup)%rmassmin    = 0.0_f
    
    ! Set optional parameters.
    if (present(irhswell))   carma%group(igroup)%irhswell     = irhswell
    if (present(irhswcomp))  carma%group(igroup)%irhswcomp    = irhswcomp
    if (present(refidx))     carma%group(igroup)%refidx(:)    = refidx(:)
    if (present(do_mie))     carma%group(igroup)%do_mie       = do_mie
    if (present(do_wetdep))  carma%group(igroup)%do_wetdep    = do_wetdep
    if (present(do_drydep))  carma%group(igroup)%do_drydep    = do_drydep
    if (present(do_vtran))   carma%group(igroup)%grp_do_vtran = do_vtran
    if (present(solfac))     carma%group(igroup)%solfac       = solfac
    if (present(scavcoef))   carma%group(igroup)%scavcoef     = scavcoef
    if (present(shortname))  carma%group(igroup)%shortname    = shortname
    if (present(cnsttype))   carma%group(igroup)%cnsttype     = cnsttype
    if (present(maxbin))     carma%group(igroup)%maxbin       = maxbin
    if (present(ifallrtn))   carma%group(igroup)%ifallrtn     = ifallrtn
    if (present(is_cloud))   carma%group(igroup)%is_cloud     = is_cloud
    if (present(rmassmin))   carma%group(igroup)%rmassmin     = rmassmin

    
    ! Initialize other properties.
    carma%group(igroup)%nelem         = 0
    carma%group(igroup)%if_sec_mom    = .FALSE.
    carma%group(igroup)%ncore         = 0
    carma%group(igroup)%ienconc       = 0
    carma%group(igroup)%imomelem      = 0
    
    return
  end subroutine CARMAGROUP_Create
    

  !! Deallocates the memory associated with a CARMAGROUP object.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  !!
  !! @see CARMAGROUP_Create
  subroutine CARMAGROUP_Destroy(carma, igroup, rc)
    type(carma_type), intent(inout)      :: carma         !! the carma object
    integer, intent(in)                  :: igroup        !! the group index
    integer, intent(out)                 :: rc            !! return code, negative indicates failure

    ! Local variables
    integer                              :: ier
    
    ! Assume success.
    rc = RC_OK
    
    ! Make sure there are enough groups allocated.
    if (igroup > carma%NGROUP) then
      if (carma%do_print) write(carma%LUNOPRT, *) "CARMAGROUP_Destroy:: ERROR - The specifed group (", &
        igroup, ") is larger than the number of groups (", carma%NGROUP, ")."
      rc = RC_ERROR
      return
    end if
    
    if (allocated(carma%group(igroup)%refidx)) then
      deallocate( &
        carma%group(igroup)%refidx, &
        carma%group(igroup)%qext, &
        carma%group(igroup)%ssa, &
        carma%group(igroup)%asym, &
        stat=ier) 
      if(ier /= 0) then
        if (carma%do_print) write(carma%LUNOPRT, *) "CARMAGROUP_Destroy: ERROR deallocating, status=", ier
        rc = RC_ERROR
        return
      endif
    endif
    
    ! Allocate dynamic data.
    if (allocated(carma%group(igroup)%r)) then
      deallocate( &
        carma%group(igroup)%r, &
        carma%group(igroup)%rmass, &
        carma%group(igroup)%vol, &
        carma%group(igroup)%dr, &
        carma%group(igroup)%dm, &
        carma%group(igroup)%rmassup, &
        carma%group(igroup)%rup, &
        carma%group(igroup)%rlow, &
        carma%group(igroup)%icorelem, &
        stat=ier) 
      if(ier /= 0) then
        if (carma%do_print) write(carma%LUNOPRT, *) "CARMAGROUP_Destroy: ERROR deallocating, status=", ier
        rc = RC_ERROR
        return
      endif
    endif

    return
  end subroutine CARMAGROUP_Destroy


  !! Gets information about a group.
  !!
  !! The group name and most other properties are available after a call to
  !! CARMAGROUP_Create(). After a call to CARMA_Initialize(), the bin
  !! dimensions and optical properties can be retrieved.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  !!
  !! @see CARMAGROUP_Create
  !! @see CARMA_GetGroup
  !! @see CARMA_Initialize 
  subroutine CARMAGROUP_Get(carma, igroup, rc, name, shortname, rmin, rmrat, ishape, eshape, is_ice, &
      irhswell, irhswcomp, cnsttype, r, rlow, rup, dr, rmass, dm, vol, qext, ssa, asym, do_mie, &
      do_wetdep, do_drydep, do_vtran, solfac, scavcoef, ienconc, refidx, ncore, icorelem, maxbin, ifallrtn, is_cloud, rmassmin)
      
    type(carma_type), intent(in)              :: carma                        !! the carma object
    integer, intent(in)                       :: igroup                       !! the group index
    integer, intent(out)                      :: rc                           !! return code, negative indicates failure
    character(len=*), optional, intent(out)   :: name                         !! the group name
    character(len=*), optional, intent(out)   :: shortname                    !! the group short name
    real(kind=f), optional, intent(out)       :: rmin                         !! the minimum radius [cm]
    real(kind=f), optional, intent(out)       :: rmrat                        !! the volume ratio between bins
    integer, optional, intent(out)            :: ishape                       !! the type of the particle shape
    real(kind=f), optional, intent(out)       :: eshape                       !! the aspect ratio of the particle shape
    logical, optional, intent(out)            :: is_ice                       !! is this an ice particle?
    integer, optional, intent(out)            :: irhswell                     !! the parameterization for particle swelling from relative humidity
    integer, optional, intent(out)            :: irhswcomp                    !! the composition for particle swelling from relative humidity
    integer, optional, intent(out)            :: cnsttype                     !! constituent type in the parent model
    real(kind=f), intent(out), optional       :: r(carma%NBIN)                !! the bin radius [cm]
    real(kind=f), intent(out), optional       :: rlow(carma%NBIN)             !! the bin radius lower bound [cm]
    real(kind=f), intent(out), optional       :: rup(carma%NBIN)              !! the bin radius upper bound [cm]
    real(kind=f), intent(out), optional       :: dr(carma%NBIN)               !! the bin width in radius space [cm]
    real(kind=f), intent(out), optional       :: rmass(carma%NBIN)            !! the bin mass [g]
    real(kind=f), intent(out), optional       :: dm(carma%NBIN)               !! the bin width in mass space [g]
    real(kind=f), intent(out), optional       :: vol(carma%NBIN)              !! the bin volume [cm<sup>3</sup>]
    complex(kind=f), intent(out), optional    :: refidx(carma%NWAVE)          !! the refractive index at each wavelength
    real(kind=f), intent(out), optional       :: qext(carma%NWAVE,carma%NBIN) !! extinction efficiency
    real(kind=f), intent(out), optional       :: ssa(carma%NWAVE,carma%NBIN)  !! single scattering albedo
    real(kind=f), intent(out), optional       :: asym(carma%NWAVE,carma%NBIN) !! asymmetry factor
    logical, optional, intent(out)            :: do_mie                       !! do mie calculations?
    logical, optional, intent(out)            :: do_wetdep                    !! do wet deposition for this particle?
    logical, optional, intent(out)            :: do_drydep                    !! do dry deposition for this particle?
    logical, optional, intent(out)            :: do_vtran                     !! do sedimentation for this particle?
    real(kind=f), intent(out), optional       :: solfac                       !! the solubility factor for wet deposition
    real(kind=f), intent(out), optional       :: scavcoef                     !! the scavenging coefficient for wet deposition
    integer, intent(out), optional            :: ienconc                      !! Particle number conc. element for group
    integer, intent(out), optional            :: ncore                        !! Number of core mass elements for group
    integer, intent(out), optional            :: icorelem(carma%NELEM)        !! Element index of core mass elements for group
    integer, optional, intent(out)            :: maxbin                       !! the last prognostic bin in the group
    integer, optional, intent(out)            :: ifallrtn                     !! fall velocity routine [I_FALLRTN_STD | I_FALLRTN_STD_SHAPE | I_FALLRTN_HEYMSFIELD2010 | I_FALLRTN_ACKERMAN_DROP | I_FALLRTN_ACKERMAN_ICE]
    logical, optional, intent(out)            :: is_cloud                     !! is this a cloud particle?
    real(kind=f), optional, intent(out)       :: rmassmin                     !! the minimum mass [g]

    ! Assume success.
    rc = RC_OK

    ! Make sure there are enough groups allocated.
    if (igroup > carma%NGROUP) then
      if (carma%do_print) write(carma%LUNOPRT, *) "CARMAGROUP_Get:: ERROR - The specifed group (", &
        igroup, ") is larger than the number of groups (", carma%NGROUP, ")."
      rc = RC_ERROR
      return
    end if
      
    ! Return any requested properties of the group.
    if (present(name))         name         = carma%group(igroup)%name
    if (present(shortname))    shortname    = carma%group(igroup)%shortname
    if (present(rmin))         rmin         = carma%group(igroup)%rmin
    if (present(rmrat))        rmrat        = carma%group(igroup)%rmrat
    if (present(ishape))       ishape       = carma%group(igroup)%ishape
    if (present(eshape))       eshape       = carma%group(igroup)%eshape
    if (present(is_ice))       is_ice       = carma%group(igroup)%is_ice
    if (present(irhswell))     irhswell     = carma%group(igroup)%irhswell
    if (present(irhswcomp))    irhswcomp    = carma%group(igroup)%irhswcomp
    if (present(cnsttype))     cnsttype     = carma%group(igroup)%cnsttype
    if (present(r))            r(:)         = carma%group(igroup)%r(:)
    if (present(rlow))         rlow(:)      = carma%group(igroup)%rlow(:)
    if (present(rup))          rup(:)       = carma%group(igroup)%rup(:)
    if (present(dr))           dr(:)        = carma%group(igroup)%dr(:)
    if (present(rmass))        rmass(:)     = carma%group(igroup)%rmass(:)
    if (present(dm))           dm(:)        = carma%group(igroup)%dm(:)
    if (present(vol))          vol(:)       = carma%group(igroup)%vol(:)
    if (present(do_mie))       do_mie       = carma%group(igroup)%do_mie
    if (present(do_wetdep))    do_wetdep    = carma%group(igroup)%do_wetdep
    if (present(do_drydep))    do_drydep    = carma%group(igroup)%do_drydep
    if (present(do_vtran))     do_vtran     = carma%group(igroup)%grp_do_vtran
    if (present(solfac))       solfac       = carma%group(igroup)%solfac
    if (present(scavcoef))     scavcoef     = carma%group(igroup)%scavcoef
    if (present(ienconc))      ienconc      = carma%group(igroup)%ienconc
    if (present(ncore))        ncore        = carma%group(igroup)%ncore
    if (present(icorelem))     icorelem     = carma%group(igroup)%icorelem(:)
    if (present(maxbin))       maxbin       = carma%group(igroup)%maxbin
    if (present(ifallrtn))     ifallrtn     = carma%group(igroup)%ifallrtn
    if (present(is_cloud))     is_cloud     = carma%group(igroup)%is_cloud
    if (present(rmassmin))     rmassmin     = carma%group(igroup)%rmassmin
    
    if (carma%NWAVE == 0) then
      if (present(refidx) .or. present(qext) .or. present(ssa) .or. present(asym)) then
        if (carma%do_print) write(carma%LUNOPRT, *) "CARMAGROUP_Get: ERROR no optical properties defined."
        rc = RC_ERROR
        return
      end if
    else
      if (present(refidx))     refidx(:)    = carma%group(igroup)%refidx(:)
      if (present(qext))       qext(:,:)    = carma%group(igroup)%qext(:,:)
      if (present(ssa))        ssa(:,:)     = carma%group(igroup)%ssa(:,:)
      if (present(asym))       asym(:,:)    = carma%group(igroup)%asym(:,:)
    end if
    
    return
  end subroutine CARMAGROUP_Get
  
  
  
  !! Prints information about a group.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  !!
  !! @see CARMAGROUP_Get
  subroutine CARMAGROUP_Print(carma, igroup, rc)
    type(carma_type), intent(in)              :: carma              !! the carma object
    integer, intent(in)                       :: igroup             !! the group index
    integer, intent(out)                      :: rc                 !! return code, negative indicates failure
    
    ! Local variables
    integer                                   :: i
    character(len=CARMA_NAME_LEN)             :: name               ! name
    character(len=CARMA_SHORT_NAME_LEN)       :: shortname          ! shortname
    real(kind=f)                              :: rmin               ! the minimum radius [cm]
    real(kind=f)                              :: rmrat              ! the volume ratio between bins
    integer                                   :: ishape             ! the type of the particle shape
    real(kind=f)                              :: eshape             ! the aspect ratio of the particle shape
    logical                                   :: is_ice             ! is this an ice particle?
    integer                                   :: irhswell           ! the parameterization for particle swelling from relative humidity
    integer                                   :: irhswcomp          ! the composition for particle swelling from relative humidity
    integer                                   :: cnsttype           ! constituent type in the parent model
    real(kind=f)                              :: r(carma%NBIN)      ! the bin radius [m]
    real(kind=f)                              :: dr(carma%NBIN)     ! the bin width in radius space [m]
    real(kind=f)                              :: rmass(carma%NBIN)  ! the bin mass [kg]
    real(kind=f)                              :: dm(carma%NBIN)     ! the bin width in mass space [kg]
    real(kind=f)                              :: vol(carma%NBIN)    ! the bin volume [m<sup>3</sup>]
    integer                                   :: ifallrtn           ! fall velocity routine [I_FALLRTN_STD | I_FALLRTN_STD_SHAPE | I_FALLRTN_HEYMSFIELD2010 | I_FALLRTN_ACKERMAN_DROP | I_FALLRTN_ACKERMAN_ICE]
    logical                                   :: is_cloud           ! is this a cloud particle?
    real(kind=f)                              :: rmassmin           ! the minimum mass [g]

    ! Assume success.
    rc = RC_OK

    ! Test out the Get method.
    if (carma%do_print) then
      call CARMAGROUP_Get(carma, igroup, rc, name=name, shortname=shortname, rmin=rmin, rmrat=rmrat, ishape=ishape, eshape=eshape, &
                        is_ice=is_ice, is_cloud=is_cloud, irhswell=irhswell, irhswcomp=irhswcomp, cnsttype=cnsttype, r=r, dr=dr, &
                        rmass=rmass, dm=dm, vol=vol, ifallrtn=ifallrtn, rmassmin=rmassmin)
      if (rc < 0) return

    
      write(carma%LUNOPRT,*) "    name          : ", trim(name)
      write(carma%LUNOPRT,*) "    shortname     : ", trim(shortname)
      write(carma%LUNOPRT,*) "    rmin          : ", rmin, " (cm)"
      write(carma%LUNOPRT,*) "    rmassmin      : ", rmassmin, " (g)"
      write(carma%LUNOPRT,*) "    rmrat         : ", rmrat

      select case(ishape)
        case (I_SPHERE)
          write(carma%LUNOPRT,*) "    ishape        :    spherical"
        case (I_HEXAGON)
          write(carma%LUNOPRT,*) "    ishape        :    hexagonal"
        case (I_CYLINDER)
          write(carma%LUNOPRT,*) "    ishape        :    cylindrical"
        case default
          write(carma%LUNOPRT,*) "    ishape        :    unknown, ", ishape
      end select

      write(carma%LUNOPRT,*) "    eshape        : ", eshape
      write(carma%LUNOPRT,*) "    is_ice        : ", is_ice
      write(carma%LUNOPRT,*) "    is_cloud      : ", is_cloud
      
      select case(irhswell)
        case (0)
          write(carma%LUNOPRT,*) "    irhswell      :    none"
        case (I_FITZGERALD)
          write(carma%LUNOPRT,*) "    irhswell      :    Fitzgerald"
        case (I_GERBER)
          write(carma%LUNOPRT,*) "    irhswell      :    Gerber"
        case default
          write(carma%LUNOPRT,*) "    irhswell      :    unknown, ", irhswell
      end select

      select case(irhswcomp)
        case (0)
          write(carma%LUNOPRT,*) "    irhswcomp     :    none"

        case (I_SWF_NH42SO4)
          write(carma%LUNOPRT,*) "    irhswcomp     :    (NH4)2SO4 (Fitzgerald)"
        case (I_SWF_NH4NO3)
          write(carma%LUNOPRT,*) "    irhswcomp     :    NH4NO3 (Fitzgerald)"
        case (I_SWF_NANO3)
          write(carma%LUNOPRT,*) "    irhswcomp     :    NaNO3 (Fitzgerald)"
        case (I_SWF_NH4CL)
          write(carma%LUNOPRT,*) "    irhswcomp     :    NH4Cl (Fitzgerald)"
        case (I_SWF_CACL2)
          write(carma%LUNOPRT,*) "    irhswcomp     :    CaCl2 (Fitzgerald)"
        case (I_SWF_NABR)
          write(carma%LUNOPRT,*) "    irhswcomp     :    NaBr (Fitzgerald)"
        case (I_SWF_NACL)
          write(carma%LUNOPRT,*) "    irhswcomp     :    NaCl (Fitzgerald)"
        case (I_SWF_MGCL2)
          write(carma%LUNOPRT,*) "    irhswcomp     :    MgCl2 (Fitzgerald)"
        case (I_SWF_LICL)
          write(carma%LUNOPRT,*) "    irhswcomp     :    LiCl (Fitzgerald)"

        case (I_SWG_NH42SO4)
          write(carma%LUNOPRT,*) "    irhswcomp     :    (NH4)2SO4 (Gerber)"
        case (I_SWG_RURAL)
          write(carma%LUNOPRT,*) "    irhswcomp     :    Rural (Gerber)"
        case (I_SWG_SEA_SALT)
          write(carma%LUNOPRT,*) "    irhswcomp     :    Sea Salt (Gerber)"
        case (I_SWG_URBAN)
          write(carma%LUNOPRT,*) "    irhswcomp     :    Urban (Gerber)"

        case default
          write(carma%LUNOPRT,*) "    irhswell      :    unknown, ", irhswcomp
      end select
      
      select case(cnsttype)
        case (0)
          write(carma%LUNOPRT,*) "    cnsttype      :    none"
        case (I_CNSTTYPE_PROGNOSTIC)
          write(carma%LUNOPRT,*) "    cnsttype      :    prognostic"
         case (I_CNSTTYPE_DIAGNOSTIC)
          write(carma%LUNOPRT,*) "    cnsttype      :    diagnostic"
        case default
          write(carma%LUNOPRT,*) "    cnsttype      :    unknown, ", cnsttype
      end select

      select case(ifallrtn)
        case (I_FALLRTN_STD)
          write(carma%LUNOPRT,*) "    cnsttype      :    standard"
        case (I_FALLRTN_STD_SHAPE)
          write(carma%LUNOPRT,*) "    cnsttype      :    standard (shape)"
        case default
          write(carma%LUNOPRT,*) "    cnsttype      :    unknown, ", cnsttype
      end select

  
      write(carma%LUNOPRT,*)   
      write(carma%LUNOPRT,"('    ', a4, 5a12)") "bin",  "r",  "dr",  "rmass",  "dm",  "vol"
      write(carma%LUNOPRT,"('    ', a4, 5a12)") "",  "(cm)",  "(cm)",  "(g)",  "(g)",  "(cm3)"
     
      do i = 1, carma%NBIN
        write(carma%LUNOPRT, "('    ', i4,  5g12.3)") i, r(i), dr(i), rmass(i), dm(i), vol(i)
      end do
    end if
    
    return
  end subroutine CARMAGROUP_Print
end module
