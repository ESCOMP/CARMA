!! The CARMAELEMENT module contains configuration information about a particle
!! element used by CARMA.
!!
!!  @version March-2010
!!  @author  Chuck Bardeen 
module CARMAELEMENT_mod

  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod

  ! CARMA explicitly declares all variables. 
  implicit none

  ! All CARMA variables and procedures are private except those explicitly declared to be public.
  private

  ! Declare the public methods.
  public CARMAELEMENT_Create
  public CARMAELEMENT_Destroy
  public CARMAELEMENT_Get
  public CARMAELEMENT_Print

contains

  !! Defines a gas used by CARMA for nucleation and growth of cloud and 
  !! aerosol particles.
  !!
  !! NOTE: The element density can be specifeid per bin using rhobin; however,
  !! if only the bulk density is provided (rho) then the same value will be used
  !! for all bins. The bulk density allows for backward compatability and ease of
  !! configuration. If rhobin is provided, then rho is ignored.
  !!
  !! @author  Chuck Bardeen
  !! @version March-2010
  !!
  !! @see CARMA_AddGas
  !! @see CARMAELEMENT_Destroy
 subroutine CARMAELEMENT_Create(carma, ielement, igroup, name, rho, itype, icomposition, rc, shortname, isolute, rhobin)
    type(carma_type), intent(inout)       :: carma               !! the carma object
    integer, intent(in)                   :: ielement            !! the element index
    integer, intent(in)                   :: igroup              !! Group to which the element belongs
    character(*), intent(in)              :: name                !! the element name, maximum of 255 characters
    real(kind=f), intent(in)              :: rho                 !! bulk mass density of particle element [g/cm^3]
    integer, intent(in)                   :: itype               !! Particle type specification
    integer, intent(in)                   :: icomposition        !! Particle compound specification
    integer, intent(out)                  :: rc                  !! return code, negative indicates failure
    character(*), optional, intent(in)    :: shortname           !! the element shortname, maximum of 6 characters
    integer, optional, intent(in)         :: isolute             !! Index of solute for the particle element
    real(kind=f), optional, intent(in)    :: rhobin(carma%NBIN)  !! mass density per bin of particle element [g/cm^3]

    ! Local variables
    integer                               :: ier
    
    ! Assume success.
    rc = RC_OK
    
    ! Make sure there are enough elements allocated.
    if (ielement > carma%NELEM) then
      if (carma%do_print) write(carma%LUNOPRT, *) "CARMAELEMENT_Create:: ERROR - The specifed element (", &
        ielement, ") is larger than the number of elements (", carma%NELEM, ")."
      rc = RC_ERROR
      return
    end if
    
    ! Make sure there are enough groups allocated.
    if (igroup > carma%NGROUP) then
      if (carma%do_print) write(carma%LUNOPRT, *) "CARMAELEMENT_Create:: ERROR - The specifed group (", &
        igroup, ") is larger than the number of groups (", carma%NGROUP, ")."
      rc = RC_ERROR
      return
    end if
    
    allocate( &
      carma%element(ielement)%rho(carma%NBIN), &
      stat=ier) 
    if(ier /= 0) then
        if (carma%do_print) write(carma%LUNOPRT, *) "CARMAELEMENT_Add: ERROR allocating, status=", ier
      rc = RC_ERROR
      return
    end if

    ! Save off the settings.
    carma%element(ielement)%igroup       = igroup
    carma%element(ielement)%name         = name
    carma%element(ielement)%rho(:)       = rho
    carma%element(ielement)%itype        = itype
    carma%element(ielement)%icomposition = icomposition
    
    
    ! Defaults for optional parameters
    carma%element(ielement)%shortname   = ""
    carma%element(ielement)%isolute     = 0
    
    ! Set optional parameters.
    if (present(shortname))  carma%element(ielement)%shortname    = shortname
    if (present(isolute)) then
    
      ! Make sure there are enough solutes allocated.
      if (isolute > carma%NSOLUTE) then
        if (carma%do_print) write(carma%LUNOPRT, *) "CARMAELEMENT_Create:: ERROR - The specifed solute (", &
          isolute, ") is larger than the number of solutes (", carma%NSOLUTE, ")."
        rc = RC_ERROR
        return
      end if

      carma%element(ielement)%isolute      = isolute
    end if
    if (present(rhobin)) carma%element(ielement)%rho(:) = rhobin(:)
    
    ! Keep track of the fact that another element has been added to the group.
    carma%group(igroup)%nelem = carma%group(igroup)%nelem + 1
    
    return
  end subroutine CARMAELEMENT_Create
    

  !! Deallocates the memory associated with a CARMAELEMENT object.
  !!
  !! @author  Chuck Bardeen
  !! @version March-2010
  !!
  !! @see CARMAELEMENT_Create
  subroutine CARMAELEMENT_Destroy(carma, ielement, rc)
    type(carma_type), intent(inout)        :: carma         !! the carma object
    integer, intent(in)                    :: ielement      !! the element index
    integer, intent(out)                   :: rc            !! return code, negative indicates failure

    ! Local variables
    integer                               :: ier

    ! Assume success.
    rc = RC_OK
    
    ! Make sure there are enough elements allocated.
    if (ielement > carma%NELEM) then
      if (carma%do_print) write(carma%LUNOPRT, *) "CARMAELEMENT_Destroy:: ERROR - The specifed element (", &
        ielement, ") is larger than the number of elements (", carma%NELEM, ")."
      rc = RC_ERROR
      return
    end if

    if (allocated(carma%element(ielement)%rho)) then
      deallocate( &
        carma%element(ielement)%rho, &
        stat=ier) 
      if(ier /= 0) then
        if (carma%do_print) write(carma%LUNOPRT, *) "CARMAELEMENT_Destroy: ERROR deallocating, status=", ier
        rc = RC_ERROR
        return
      endif
    endif

    return
  end subroutine CARMAELEMENT_Destroy


  !! Gets information about a particle element.
  !!
  !! The group name and other properties are available after a call to
  !! CARMAELEMENT_Create().
  !!
  !! @author  Chuck Bardeen
  !! @version March-2010
  !!
  !! @see CARMAELEMENT_Create
  !! @see CARMA_GetElement
  subroutine CARMAELEMENT_Get(carma, ielement, rc, igroup, name, shortname, rho, itype, icomposition, isolute)
    type(carma_type), intent(in)                :: carma           !! the carma object
    integer, intent(in)                         :: ielement        !! the element index
    integer, intent(out)                        :: rc              !! return code, negative indicates failure
    integer, optional, intent(out)              :: igroup          !! Group to which the element belongs
    character(len=*), optional, intent(out)     :: name            !! the element name
    character(len=*), optional, intent(out)     :: shortname       !! the element short name
    real(kind=f), optional, intent(out)         :: rho(CARMA%NBIN) !! Mass density of particle element [g/cm^3]
    integer, optional, intent(out)              :: itype           !! Particle type specification
    integer, optional, intent(out)              :: icomposition    !! Particle compound specification
    integer, optional, intent(out)              :: isolute         !! Index of solute for the particle element
    
    ! Assume success.
    rc = RC_OK

    ! Make sure there are enough elements allocated.
    if (ielement > carma%NELEM) then
      if (carma%do_print) write(carma%LUNOPRT, *) "CARMAELEMENT_Get:: ERROR - The specifed element (", &
        ielement, ") is larger than the number of elements (", carma%NELEM, ")."
      rc = RC_ERROR
      return
    end if

    ! Return any requested properties of the group.
    if (present(igroup))       igroup       = carma%element(ielement)%igroup
    if (present(name))         name         = carma%element(ielement)%name
    if (present(shortname))    shortname    = carma%element(ielement)%shortname
    if (present(rho))          rho(:)       = carma%element(ielement)%rho(:)
    if (present(itype))        itype        = carma%element(ielement)%itype
    if (present(icomposition)) icomposition = carma%element(ielement)%icomposition
    if (present(isolute))      isolute      = carma%element(ielement)%isolute
        
    return
  end subroutine CARMAELEMENT_Get
  
  
  !! Prints information about an element.
  !!
  !! @author  Chuck Bardeen
  !! @version March-2010
  !!
  !! @see CARMAELEMENT_Get
  subroutine CARMAELEMENT_Print(carma, ielement, rc)
    type(carma_type), intent(in)              :: carma         !! the carma object
    integer, intent(in)                       :: ielement      !! the element index
    integer, intent(out)                      :: rc            !! return code, negative indicates failure
    
    ! Local variables
    character(len=CARMA_NAME_LEN)             :: name             ! name
    character(len=CARMA_SHORT_NAME_LEN)       :: shortname        ! shortname
    real(kind=f)                              :: rho(carma%NBIN)  ! density (g/cm3)
    integer                                   :: igroup           ! Group to which the element belongs
    integer                                   :: itype            ! Particle type specification
    integer                                   :: icomposition     ! Particle compound specification
    integer                                   :: isolute          ! Index of solute for the particle element

    ! Assume success.
    rc = RC_OK

    ! Test out the Get method.
    if (carma%do_print) then
      call CARMAELEMENT_Get(carma, ielement, rc, name=name, shortname=shortname, igroup=igroup, &
                            itype=itype, icomposition=icomposition, rho=rho, isolute=isolute) 
      if (rc < 0) return

    
      write(carma%LUNOPRT,*) "    name          : ", trim(name)
      write(carma%LUNOPRT,*) "    igroup        : ", igroup
      write(carma%LUNOPRT,*) "    shortname     : ", trim(shortname)
      write(carma%LUNOPRT,*) "    rho           : ", rho, " (g/cm3)"

      select case(itype)
        case (I_INVOLATILE)
          write(carma%LUNOPRT,*) "    itype         :    involatile"
        case (I_VOLATILE)
          write(carma%LUNOPRT,*) "    itype         :    volatile"
        case (I_COREMASS)
          write(carma%LUNOPRT,*) "    itype         :    core mass"
        case (I_VOLCORE)
          write(carma%LUNOPRT,*) "    itype         :    volatile core"
        case (I_CORE2MOM)
          write(carma%LUNOPRT,*) "    itype         :    core mass - second moment"
        case default
          write(carma%LUNOPRT,*) "    itype         :    unknown, ", itype
      end select

      write(carma%LUNOPRT,*) "    icomposition  : ", icomposition
      write(carma%LUNOPRT,*) "    isolute       : ", isolute
    end if
    
    return
  end subroutine CARMAELEMENT_Print
end module
