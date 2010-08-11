!! This code is to test the error handling in the CARMA configuration and
!! initialization interface.
!!
!! @author Chuck Bardeen
!! @version Mar-2009

program carma_inittest
  implicit none

  write(*,*) "CARMA Initializtion Test"

  call test_initialization()  
end program


subroutine test_initialization()
  use carma_precision_mod 
  use carma_constants_mod 
  use carma_enums_mod 
  use carma_types_mod 
  use carmaelement_mod
  use carmagroup_mod
  use carmastate_mod
  use carma_mod
  use atmosphere_mod

  implicit none

  integer, parameter    :: NX           = 1
  integer, parameter    :: NY           = 1
  integer, parameter    :: NZ           = 80
  integer, parameter    :: NZP1         = NZ+1
  integer, parameter    :: NELEM        = 3
  integer, parameter    :: NBIN         = 20
  integer, parameter    :: NGROUP       = 2
  integer, parameter    :: NSOLUTE      = 0
  integer, parameter    :: NGAS         = 0
  integer, parameter    :: NWAVE        = 0
  integer, parameter    :: nstep        = 72

  real(kind=f), parameter   :: dtime  = 600._f
  real(kind=f), parameter   :: deltax = 100._f
  real(kind=f), parameter   :: deltay = 100._f
  real(kind=f), parameter   :: deltaz = 100._f
  real(kind=f), parameter   :: zmin   = 0._f
  
  integer, parameter        :: I_BLACKCARBON       = 1
  integer, parameter        :: I_ORGANICCARBON     = 2

  type(carma_type), target  :: carma
  type(carma_type), pointer :: carma_ptr
  type(carmastate_type), allocatable     :: cstate(:)
  integer                   :: rc = 0
  
  real(kind=f), allocatable   :: xc(:,:,:)
  real(kind=f), allocatable   :: dx(:,:,:)
  real(kind=f), allocatable   :: yc(:,:,:)
  real(kind=f), allocatable   :: dy(:,:,:)
  real(kind=f), allocatable   :: zc(:,:,:)
  real(kind=f), allocatable   :: zl(:,:,:)
  real(kind=f), allocatable   :: p(:,:,:)
  real(kind=f), allocatable   :: pl(:,:,:)
  real(kind=f), allocatable   :: t(:,:,:)
  real(kind=f), allocatable   :: rhoa(:,:,:)
  
  real(kind=f), allocatable, target  :: mmr(:,:,:,:,:)
  
  real(kind=f), allocatable          :: lat(:,:)
  real(kind=f), allocatable          :: lon(:,:)
  
  integer               :: i, j
  integer               :: ix
  integer               :: iy
  integer               :: ixy
  integer               :: istep
  integer               :: ielem
  integer               :: ibin
  integer               :: ithread
  integer, parameter    :: lun = 42

  real(kind=f)          :: time
  real(kind=f)          :: rmin, rmrat, ck0
  real(kind=f)          :: r(NBIN)
  real(kind=f)          :: dr(NBIN)
  real(kind=f)          :: rmass(NBIN)
    

  write(*,*) ""

  ! Allocate the arrays that we need for the model
  allocate(xc(NZ,NY,NX), dx(NZ,NY,NX), yc(NZ,NY,NX), dy(NZ,NY,NX), &
           zc(NZ,NY,NX), zl(NZP1,NY,NX), p(NZ,NY,NX), pl(NZP1,NY,NX), &
           t(NZ,NY,NX), rhoa(NZ,NY,NX))
  allocate(mmr(NZ,NY,NX,NELEM,NBIN))
  allocate(lat(NY,NX), lon(NY,NX))  

  ! Define the particle-grid extent of the CARMA test
  call CARMA_Create(carma, NBIN, NELEM, NGROUP, NSOLUTE, NGAS, NWAVE, rc, LUNOPRT=6)
  if (rc /=0) stop "    *** FAILED ***"
	carma_ptr => carma

  ! Define the groups
  ! -----------------
  write(*,*) "  Add Group(s) ..."
  rmrat = 2._f
  rmin = 3.e-7_f
  call CARMAGROUP_Create(carma, 1, 'mixed carbon', rmin, rmrat, I_SPHERE, 1._f, .FALSE., rc)
  if (rc /=0) stop "    *** FAILED ***"

  call CARMAGROUP_Create(carma, 2, 'organic carbon', rmin, rmrat, I_SPHERE, 1._f, .FALSE., rc)
  if (rc /=0) stop "    *** FAILED ***"

  
  ! Define the elements
  ! -------------------
  write(*,*) "  Add Element(s) ..."

  ! Organic Carbon
  call CARMAELEMENT_Create(carma, 1, 1, "bc in mixed", 1._f, I_INVOLATILE, I_BLACKCARBON, rc)
  if (rc /=0) write(*, *) "    *** FAILED ***, rc=", rc

  call CARMAELEMENT_Create(carma, 2, 1, "oc in mixed", 1._f, I_COREMASS, I_ORGANICCARBON, rc)
  if (rc /=0) write(*, *) "    *** FAILED ***, rc=", rc

  call CARMAELEMENT_Create(carma, 3, 2, "organic carbon", 1._f, I_INVOLATILE, I_ORGANICCARBON, rc)
  if (rc /=0) write(*, *) "    *** FAILED ***, rc=", rc


  ! AddCoagulation Tests ...  
  write(*,*) ""
  write(*,*) "  Check for group too large in AddCoagulation ... Should Fail"
  write(*,*) ""
  call CARMA_AddCoagulation(carma, 4, 1, 1, I_COLLEC_DATA, rc) 
  if (rc /=0) then
    rc = 0
  else
    stop "    *** FAILED ***"
  endif
  call CARMA_AddCoagulation(carma, 1, 4, 1, I_COLLEC_DATA, rc) 
  if (rc /=0) then
    rc = 0
  else
    stop "    *** FAILED ***"
  endif
  call CARMA_AddCoagulation(carma, 1, 1, 4, I_COLLEC_DATA, rc) 
  if (rc /=0) then
    rc = 0
  else
    stop "    *** FAILED ***"
  endif
  
  ! Initialization Tests
  write(*,*) ""
  write(*,*) "  Check for order of element list ... Should Fail"
  write(*,*) ""
  call CARMA_AddCoagulation(carma, 1, 1, 1, I_COLLEC_DATA, rc) 
  call CARMA_AddCoagulation(carma, 2, 2, 2, I_COLLEC_DATA, rc) 
  call CARMA_AddCoagulation(carma, 2, 1, 1, I_COLLEC_DATA, rc) 
  if (rc /=0) stop "    *** FAILED ***"

  call CARMA_Initialize(carma, rc, do_coag=.TRUE.)
  if (rc /=0) then
    rc = 0
  else
    stop "    *** FAILED ***"
  endif
  

  write(*,*)  ""
  write(*,*) "  CARMA_Destroy() ..."
  call CARMA_Destroy(carma, rc)
  if (rc /=0) write(*, *) "    *** FAILED ***, rc=", rc  
end subroutine


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
	
  write(*,*)  ""
  write(*,*)  "Element Information"
  
  do i = 1, carma%NELEM
    call CARMAELEMENT_Print(carma, i, rc)
    if (rc /=0) write(carma%LUNOPRT, *) "    *** FAILED ***, rc=", rc
    
    write(carma%LUNOPRT,*) ""  
  end do
 
  write(carma%LUNOPRT,*) ""
  return
end subroutine


subroutine dumpGas(carma, rc)
  use carma_precision_mod 
  use carma_enums_mod 
  use carma_types_mod 
  use carmagas_mod
  use carma_mod

  implicit none

  type(carma_type), pointer, intent(inout)           :: carma              !! the carma object
  integer, intent(inout)                    :: rc                 !! return code, negative indicates failure
  
  ! Local Variables
  integer                       :: i
  character(len=255)            :: gasname
	real(kind=f)                  :: gwtmol
	
  write(*,*)  ""
  write(*,*)  "Gas Information"
  
  do i = 1, carma%NGAS
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
	
  write(*,*)  ""
  write(*,*)  "Group Information"
  
  do i = 1, carma%NGROUP
    call CARMAGROUP_Print(carma, i, rc)
    if (rc /=0) write(carma%LUNOPRT, *) "    *** FAILED ***, rc=", rc
    
    write(carma%LUNOPRT,*) ""  
  end do
 
  write(carma%LUNOPRT,*) ""
  return
end subroutine
