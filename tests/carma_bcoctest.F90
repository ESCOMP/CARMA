!! This code is to demonstrate the CARMA coagulation routines
!! in a test case utilizing three groups, one group that has two
!! components (BC, OC, OC/BC).
!!
!! Upon execution, a text file (carma_bcoctest.txt) is generated.
!! The text file can be read with the IDL procedure read_bcoctest.pro.
!!
!! @author Peter Colarco (based on Chuck Bardeen's code)
!! @version Feb-2009

program carma_bcoctest
  implicit none

  write(*,*) "Coagulation Test (3 groups)"

  call test_coagulation_bcoc()  
  
  write(*,*) "Done"
end program


subroutine test_coagulation_bcoc()
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
  integer, parameter    :: NELEM        = 4
  integer, parameter    :: NBIN         = 20
  integer, parameter    :: NGROUP       = 3
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
  type(carmastate_type)     :: cstate
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
  real(kind=f)          :: rmin, rmrat, rhobc, rhooc, ck0
  real(kind=f)          :: r(NBIN)
  real(kind=f)          :: dr(NBIN)
  real(kind=f)          :: rmass(NBIN)
  

!  write(*,*) ""
!  write(*,*) "Coagulation of Particles"

  ! Open the output text file
  open(unit=lun,file="carma_bcoctest.txt",status="unknown")
  
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
!  write(*,*) "  Add Group(s) ..."
  rmrat = 2._f
  rmin = 3.e-7_f
  
  rhobc = 2._f
  rhooc = 1._f
  
  ! Black Carbon
  call CARMAGROUP_Create(carma, 1, 'black carbon', rmin, rmrat, I_SPHERE, 1._f, .FALSE., rc)
  if (rc /=0) stop "    *** FAILED ***"
  
  ! Organic Carbon
  call CARMAGROUP_Create(carma, 2, 'organic carbon', rmin, rmrat, I_SPHERE, 1._f, .FALSE., rc)
  if (rc /=0) stop "    *** FAILED ***"

  ! Mixed BC/OC
  call CARMAGROUP_Create(carma, 3, 'mixed bc/oc', rmin, rmrat, I_SPHERE, 1._f, .FALSE., rc)
  if (rc /=0) stop "    *** FAILED ***"

  
  ! Define the elements
  ! -------------------
!  write(*,*) "  Add Element(s) ..."

  ! Black Carbon
  call CARMAELEMENT_Create(carma, 1, 1, "black carbon", rhobc, I_INVOLATILE, I_BLACKCARBON, rc)
  if (rc /=0) stop "    *** FAILED ***"

  ! Organic Carbon
  call CARMAELEMENT_Create(carma, 2, 2, "organic carbon", rhooc, I_INVOLATILE, I_ORGANICCARBON, rc)
  if (rc /=0) stop "    *** FAILED ***"

  ! Mixed BC/OC
  call CARMAELEMENT_Create(carma, 3, 3, "mixed bc/oc", rhooc, I_INVOLATILE, I_ORGANICCARBON, rc)
  if (rc /=0) stop "    *** FAILED ***"

  ! Mass of BC in Mixed Group
  call CARMAELEMENT_Create(carma, 4, 3, "mass of bc in mixed", rhobc, I_COREMASS, I_BLACKCARBON, rc)
  if (rc /=0) stop "    *** FAILED ***"

  
  ! Setup the coagulation mapping and kernels
  ! -----------------------------------------
!  write(*,*) "  CARMA_AddCoagulation(carma, 1, 1, 1, I_COLLEC_DATA, rc) ..."
  ! From Jacobson:
! ck0 = 8 * kB * T / 3 / dynamic viscosity of air
  ck0 = 8._f * bk * 298._f / 3._f / 1.85e-4_f
  ! Black Carbon
  call CARMA_AddCoagulation(carma, 1, 1, 1, I_COLLEC_DATA, rc, ck0=ck0) 
  if (rc /=0) stop "    *** FAILED ***"
  
  ! Organic Carbon
  call CARMA_AddCoagulation(carma, 2, 2, 2, I_COLLEC_DATA, rc, ck0=ck0) 
  if (rc /=0) stop "    *** FAILED ***"
  
  ! Organic Carbon & Black Carbon
  call CARMA_AddCoagulation(carma, 2, 1, 3, I_COLLEC_DATA, rc, ck0=ck0) 
  if (rc /=0) stop "    *** FAILED ***"
  
  ! Black Carbon & Mixed
  call CARMA_AddCoagulation(carma, 3, 1, 3, I_COLLEC_DATA, rc, ck0=ck0) 
  if (rc /=0) stop "    *** FAILED ***"
  
  ! Organic Carbon & Mixed
  call CARMA_AddCoagulation(carma, 3, 2, 3, I_COLLEC_DATA, rc, ck0=ck0) 
  if (rc /=0) stop "    *** FAILED ***"
  
  ! Mixed to Mixed
  call CARMA_AddCoagulation(carma, 3, 3, 3, I_COLLEC_DATA, rc, ck0=ck0) 
  if (rc /=0) stop "    *** FAILED ***"
  
  ! Setup the CARMA processes to exercise
!  write(*,*) "  CARMA_Initialize(carma, rc, do_vtran=.FALSE., "// &
!               "do_coag=.TRUE.) ..."
               
               
  call CARMA_Initialize(carma, rc, do_coag=.TRUE.)
  if (rc /=0) stop "    *** FAILED ***"

  ! Print the Group Information
!  write(*,*)  ""
!  call dumpGroup(carma, rc)
!  if (rc /=0) stop "    *** FAILED ***"
  
!  write(*,*) ""
  
  ! Print the Element Information
!  write(*,*)  ""
!  call dumpElement(carma, rc)
!  if (rc /=0) stop "    *** FAILED ***"
  
!  write(*,*) ""
  
  ! For simplicity of setup, do a case with Cartesian coordinates,
  ! which are specified in this interface in meters.
  !
  ! NOTE: For Cartesian coordinates, the first level is the bottom 
  ! of the model (e.g. z = 0), while for sigma and hybrid coordinates
  ! the first level is the top of the model.
  lat(:,:) = 40.0_f
  lon(:,:) = -105.0_f
  
  ! Horizonal centers
  do ix = 1, NX
    do iy = 1, NY
      dx(:,iy,ix) = deltax
      xc(:,iy,ix) = ix*dx(:,iy,ix) / 2._f
      dy(:,iy,ix) = deltay
      yc(:,iy,ix) = iy*dy(:,iy,ix) / 2._f
    end do
  end do
  
  ! Vertical center
  do i = 1, NZ
    zc(i,:,:) = zmin + (deltaz * (i - 0.5_f))
  end do
  call GetStandardAtmosphere(zc, p=p, t=t)

!  write(*,'(a6, 3a12)') "level", "zc", "p", "t"
!  write(*,'(a6, 3a12)') "", "(m)", "(Pa)", "(K)"
!  do i = 1, NZ
!    write(*,'(i6,3f12.3)') i, zc(i,NY,NX), p(i,NY,NX), t(i,NY,NX)
!  end do


  ! Vertical edge
  do i = 1, NZP1
    zl(i,:,:) = zmin + ((i - 1) * deltaz)
  end do
  call GetStandardAtmosphere(zl, p=pl)

!  write(*,*) ""
!  write(*,'(a6, 2a12)') "level", "zl", "pl"
!  write(*,'(a6, 2a12)') "", "(m)", "(Pa)"
!  do i = 1, NZP1
!    write(*,'(i6,2f12.3)') i, zl(i,NY,NX), pl(i,NY,NX)
!  end do

  			
  ! Put a monodisperse aerosol in first bin
  ! 10^12 m-3 (= 10^6 cm-3)
  ! Note: p is in MKS here, but rmass is in CGS, so scale
  mmr(:,:,:,:,:) = 0._f
  
  !initial black carbon
  call CARMAGROUP_Get(carma, 1, rc, rmass=rmass)
  if (rc /=0) stop "    *** FAILED ***"
  mmr(1,:,:,1,1) = rmass(1)/1000._f* 1.e12_f &
                 / (p(1,:,:)/287._f/t(1,:,:))
                 
  ! initial organic carbon
  call CARMAGROUP_Get(carma, 2, rc, rmass=rmass)
  if (rc /=0) stop "    *** FAILED ***"
  mmr(1,:,:,2,1) = rmass(1)/1000._f* 1.e12_f &
                 / (p(1,:,:)/287._f/t(1,:,:))

!  write(*,*)  ""
!  write(*, '(a6, 4a12)') "level", "mmr(i,NY,NX,1)", "mmr(i,NY,NX,2)"
!  do i = 1, NZ
!	  write(*, '(i6, 4g12.3)') i, mmr(i,NY,NX,1,1), mmr(i,NY,NX,1,2)
!  end do
  
  ! Write output for the coagtest (output is scaled to CGS)
  write(lun,*) NBIN, NELEM, NGROUP
  
  do j = 1, NGROUP
		call CARMAGROUP_Get(carma, j, rc, r=r, dr=dr)
		if (rc /=0) stop "    *** FAILED ***"
    
    ! Bin structure
    do i = 1, NBIN
      write(lun,'(i3,2(1x,e12.5))') &
       i, &
       r(i), &
       dr(i)
    end do
  end do

  ! Initial particle mass densities.
  write(lun,*) 0
  
  rhoa(:,:,:) = p(:,:,:)/287._f/t(:,:,:)
  
  do j = 1, NELEM
   do i = 1, NBIN
    write(lun,'(i3,1(1x,e12.5))') &
     i, &
     real(mmr(1,NY,NX,j,i)*rhoa(1,NY,NX)*1e-6_f*1e3_f)
   end do
  end do
		
  ! Iterate the model over a few time steps.
!  write(*,*) ""
  do istep = 1, nstep
  
    ! Calculate the model time.
    time = (istep - 1) * dtime

    ! NOTE: This means that there should not be any looping over NX or NY done
    ! in any other CARMA routines. They should only loop over NZ.
    do ixy = 1, NX*NY
       ix = ((ixy-1) / NY) + 1
       iy = ixy - (ix-1)*NY
			
       ! Create a CARMASTATE for this column.
       call CARMASTATE_Create(cstate, carma_ptr, time, dtime, NZ, &
                              I_CART, I_CART, lat(iy,ix), lon(iy,ix), &
                              xc(:,iy,ix), dx(:,iy,ix), &
                              yc(:,iy,ix), dy(:,iy,ix), &
                              zc(:,iy,ix), zl(:,iy,ix), p(:,iy,ix), &
                              pl(:,iy,ix), t(:,iy,ix), rc)
		
       ! Send the bin mmrs to CARMA
       do ielem = 1, NELEM
        do ibin = 1, NBIN
         call CARMASTATE_SetBin(cstate, ielem, ibin, &
                                mmr(:,iy,ix,ielem,ibin), rc)
        end do
       end do
			
       ! Execute the step
       call CARMASTATE_Step(cstate, rc)

       ! Get the updated bin mmr.
       do ielem = 1, NELEM
        do ibin = 1, NBIN
         call CARMASTATE_GetBin(cstate, ielem, ibin, &
                                mmr(:,iy,ix,ielem,ibin), rc)
        end do
       end do

       ! Get the updated temperature.
       call CARMASTATE_GetState(cstate, rc, t=t(:,iy,ix))

    enddo

    ! Write output for the coagtest (output in CGS)
    write(lun,'(f12.0)') istep*dtime
    rhoa(:,:,:) = p(:,:,:)/287._f/t(:,:,:)
    
    do j = 1, NELEM
     do i = 1, NBIN
      write(lun,'(i3,1(1x,e12.5))') &
       i, &
       real(mmr(1,NY,NX,j,i)*rhoa(1,NY,NX)*1e-6_f*1e3_f)
     end do
    end do

  end do   ! time loop

  ! Cleanup the carma state objects
  call CARMASTATE_Destroy(cstate, rc)
  if (rc /=0) stop "    *** FAILED ***"

  ! Close the output file
  close(unit=lun)	
	
!  write(*,*)  ""
!  write(*,*)  ""
!  write(*, '(a8, 8a14)') "level", "mmr(i,NY,NX,1)", "mmr(i,NY,NX,2)"
!  do i = 1, NZ
!   write(*, '(i8, 8g14.3)') i, mmr(i,NY,NX,1,1), mmr(i,NY,NX,1,2)
!  end do
		
!  write(*,*)  ""
!  write(*, '(a8, 2a12)') "level", "t(:,1,1)", "t(:,NY,NX)"		
!  do i = 1, NZ
!   write(*, '(i8, 2f12.3)') i, t(i,1,1), t(i,NY,NX)
!  end do

  if (rc /=0) stop "    *** FAILED ***"

!  write(*,*)  ""
!  write(*,*) "  CARMA_Destroy() ..."
  call CARMA_Destroy(carma, rc)
  if (rc /=0) stop "    *** FAILED ***"
end subroutine
