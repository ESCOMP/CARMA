!! Some common utilities used in the CARMA test code.
!!
!! @author Chuck Bardeen
!! @version Jun-2010


subroutine dumpElement(carma, rc)
  use carma_precision_mod 
  use carma_enums_mod 
  use carma_types_mod 
  use carmaelement_mod
  use carma_mod

  implicit none

  type(carma_type), intent(in)     :: carma              !! the carma object
  integer, intent(inout)           :: rc                 !! return code, negative indicates failure
  
  ! Local Variables
  integer                          :: i
  integer                          :: NELEM
	
  write(*,*)  ""
  write(*,*)  "Element Information"

  call CARMA_Get(carma, rc, NELEM=NELEM)
  if (rc /=0) write(*, *) "    *** FAILED ***, rc=", rc
  
  do i = 1, NELEM
    call CARMAELEMENT_Print(carma, i, rc)
    if (rc /=0) write(*, *) "    *** FAILED ***, rc=", rc
    
    write(*,*) ""  
  end do
 
  write(*,*) ""
  return
end subroutine


subroutine dumpGas(carma, rc)
  use carma_precision_mod 
  use carma_enums_mod 
  use carma_types_mod 
  use carmagas_mod
  use carma_mod

  implicit none

  type(carma_type), intent(in)   :: carma              !! the carma object
  integer, intent(inout)         :: rc                 !! return code, negative indicates failure
  
  ! Local Variables
  integer                        :: i
  integer                        :: NGAS
	
  call CARMA_Get(carma, rc, NGAS=NGAS)
  if (rc /=0) write(*, *) "    *** FAILED ***, rc=", rc

  write(*,*)  ""
  write(*,*)  "Gas Information"

  do i = 1, NGAS
   call CARMAGAS_Print(carma, i, rc)
   if (rc /=0) write(*, *) "    *** FAILED ***, rc=", rc
   
   write(*,*) ""  
 end do
 
 write(*,*) ""  
end subroutine


subroutine dumpGroup(carma, rc)
  use carma_precision_mod 
  use carma_enums_mod 
  use carma_types_mod 
  use carmagroup_mod
  use carma_mod

  implicit none

  type(carma_type), intent(in)     :: carma              !! the carma object
  integer, intent(inout)           :: rc                 !! return code, negative indicates failure
  
  ! Local Variables
  integer                          :: i
  integer                          :: NGROUP
	
  write(*,*)  ""
  write(*,*)  "Group Information"
  
  call CARMA_Get(carma, rc, NGROUP=NGROUP)
  if (rc /=0) write(*, *) "    *** FAILED ***, rc=", rc
  
  do i = 1, NGROUP
    call CARMAGROUP_Print(carma, i, rc)
    if (rc /=0) write(*, *) "    *** FAILED ***, rc=", rc
    
    write(*,*) ""  
  end do
 
  write(*,*) ""
  return
end subroutine
