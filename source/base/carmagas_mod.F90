!! The CARMAGAS module contains configuration information about a gas used by CARMA.
!!
!!  @version May-2009 
!!  @author  Chuck Bardeen 
module carmagas_mod

  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod

  ! CARMA explicitly declares all variables. 
  implicit none

  ! All CARMA variables and procedures are private except those explicitly declared to be public.
  private

  ! Declare the public methods.
  public CARMAGAS_Create
  public CARMAGAS_Destroy
  public CARMAGAS_Get
  public CARMAGAS_Print

contains

  !! Defines a gas used by CARMA for nucleation and growth of cloud and 
  !! aerosol particles.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  !!
  !! @see CARMA_AddGas
  !! @see CARMAGAS_Destroy
  subroutine CARMAGAS_Create(carma, igas, name, wtmol, ivaprtn, icomposition, rc, shortname)
    type(carma_type), intent(inout)       :: carma           !! the carma object
    integer, intent(in)                   :: igas            !! the gas index
    character(*), intent(in)              :: name            !! the gas name, maximum of 255 characters
    real(kind=f), intent(in)              :: wtmol           !! the gas molecular weight [g/mol]
    integer, intent(in)                   :: ivaprtn         !! vapor pressure routine for this gas
		integer, intent(in)                   :: icomposition    !! gas compound specification
    integer, intent(out)                  :: rc              !! return code, negative indicates failure
    character(*), optional, intent(in)    :: shortname       !! the gas shortname, maximum of 6 characters

    ! Assume success.
    rc = RC_OK
    
    ! Make sure there are enough gases allocated.
    if (igas > carma%NGAS) then
      if (carma%do_print) write(carma%LUNOPRT, *) "CARMAGAS_GetCreate:: ERROR - The specifed gas (", &
        igas, ") is larger than the number of gases (", carma%NGAS, ")."
      rc = RC_ERROR
      return
    end if

    ! Save off the settings.
    carma%gas(igas)%name         = name
    carma%gas(igas)%wtmol        = wtmol
    carma%gas(igas)%ivaprtn      = ivaprtn
    carma%gas(igas)%icomposition = icomposition
    
    
    ! Defaults for optional parameters
    carma%gas(igas)%shortname    = ""
    
    ! Set optional parameters.
    if (present(shortname))  carma%gas(igas)%shortname    = shortname

    return
  end subroutine CARMAGAS_Create
    

  !! Deallocates the memory associated with a CARMAGAS object.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  !!
  !! @see CARMAGAS_Create
  subroutine CARMAGAS_Destroy(carma, igas, rc)
    type(carma_type), intent(inout)    :: carma         !! the carma object
    integer, intent(in)                :: igas          !! the gas index
    integer, intent(out)               :: rc            !! return code, negative indicates failure

    ! Assume success.
    rc = RC_OK
    
    ! Make sure there are enough gases allocated.
    if (igas > carma%NGAS) then
      if (carma%do_print) write(carma%LUNOPRT, *) "CARMAGAS_Destroy:: ERROR - The specifed gas (", &
        igas, ") is larger than the number of gases (", carma%NGAS, ")."
      rc = RC_ERROR
      return
    end if

    return
  end subroutine CARMAGAS_Destroy


  !! Gets information about a gas.
  !!
  !! The group name and other properties are available after a call to
  !! CARMAGAS_Create().
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  !!
  !! @see CARMAGAS_Create
  !! @see CARMA_GetGas
  subroutine CARMAGAS_Get(carma, igas, rc, name, shortname, wtmol, ivaprtn, icomposition)
    type(carma_type), intent(in)                :: carma        !! the carma object
    integer, intent(in)                         :: igas         !! the gas index
    integer, intent(out)                        :: rc           !! return code, negative indicates failure
    character(len=*), optional, intent(out)     :: name         !! the gas name
    character(len=*), optional, intent(out)     :: shortname    !! the gas short name
    real(kind=f), optional, intent(out)         :: wtmol        !! the gas molecular weight [g/mol]
    integer, optional, intent(out)              :: ivaprtn      !! vapor pressure routine for this gas
		integer, optional, intent(out)              :: icomposition !! gas compound specification

    ! Assume success.
    rc = RC_OK

    ! Make sure there are enough gases allocated.
    if (igas > carma%NGAS) then
      if (carma%do_print) write(carma%LUNOPRT, *) "CARMAGAS_Get:: ERROR - The specifed gas (", &
        igas, ") is larger than the number of gases (", carma%NGAS, ")."
      rc = RC_ERROR
      return
    end if

    ! Return any requested properties of the group.
    if (present(name))         name         = carma%gas(igas)%name
    if (present(shortname))    shortname    = carma%gas(igas)%shortname
    if (present(wtmol))        wtmol        = carma%gas(igas)%wtmol
    if (present(ivaprtn))      ivaprtn      = carma%gas(igas)%ivaprtn
    if (present(icomposition)) icomposition = carma%gas(igas)%icomposition
        
    return
  end subroutine CARMAGAS_Get
  
  
  !! Prints information about a gas.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  !!
  !! @see CARMAGAS_Get
  subroutine CARMAGAS_Print(carma, igas, rc)
    type(carma_type), intent(in)              :: carma         !! the carma object
    integer, intent(in)                       :: igas          !! the gas index
    integer, intent(out)                      :: rc            !! return code, negative indicates failure
    
    ! Local variables
    character(len=CARMA_NAME_LEN)             :: name          !! name
    character(len=CARMA_SHORT_NAME_LEN)       :: shortname     !! shortname
    real(kind=f)                              :: wtmol         !! molecular weight (g/mol)
    integer                                   :: ivaprtn       !! vapor pressure routine for this gas
    integer                                   :: icomposition  !! gas compound specification

    ! Assume success.
    rc = RC_OK

    ! Test out the Get method.
    if (carma%do_print) then
      call CARMAGAS_Get(carma, igas, rc, name=name, shortname=shortname, wtmol=wtmol, &
                        ivaprtn=ivaprtn, icomposition=icomposition)
      if (rc < RC_OK) return

    
      write(carma%LUNOPRT,*) "    name          : ", trim(name)
      write(carma%LUNOPRT,*) "    shortname     : ", trim(shortname)
      write(carma%LUNOPRT,*) "    wtmol         : ", wtmol, " (g/mol)"

      select case(ivaprtn)
        case (I_VAPRTN_H2O_BUCK1981)
          write(carma%LUNOPRT,*) "    ivaprtn       :    Buck [1981]"
        case (I_VAPRTN_H2O_MURPHY2005)
          write(carma%LUNOPRT,*) "    ivaprtn       :    Murphy & Koop [2005]"
        case default
          write(carma%LUNOPRT,*) "    ivaprtn       :    unknown, ", ivaprtn
      end select

      select case(icomposition)
        case (I_GCOMP_H2O)
          write(carma%LUNOPRT,*) "    icomposition  :    H2O"
        case default
          write(carma%LUNOPRT,*) "    icomposition  :    unknown, ", icomposition
      end select
    end if
    
    return
  end subroutine CARMAGAS_Print
end module
