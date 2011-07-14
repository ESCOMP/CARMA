!! This code is to demonstrate the CARMA mie routines.
!!
!! Upon execution, a text file (carma_mietest.txt) is generated.
!! The text file can be read with the IDL procedure read_mietest.pro.
!!
!! @author  Chuck Bardeen
!! @version May-2009

program carma_mietest
  implicit none

  write(*,*) "Mie Test"

  call test_mie()  
  
  write(*,*) "Done"
end program


subroutine test_mie()
  use carma_precision_mod 
  use carma_constants_mod 
  use carma_enums_mod 
  use carma_types_mod 
  use carmaelement_mod
  use carmagroup_mod
  use carmastate_mod
  use carma_mod

  implicit none

  integer, parameter    :: NELEM        = 1
  integer, parameter    :: NBIN         = 16
  integer, parameter    :: NGROUP       = 1
  integer, parameter    :: NSOLUTE      = 0
  integer, parameter    :: NGAS         = 0
  integer, parameter    :: NWAVE        = 44
  
  integer, parameter    :: I_DUST       = 1

  type(carma_type), target            :: carma
  type(carma_type), pointer           :: carma_ptr
  integer                             :: rc = 0
  

  integer               :: i
  integer               :: ielem
  integer               :: iwave
  integer               :: ibin
  integer               :: igroup
  integer, parameter    :: lun = 42

  real(kind=f)          :: rmin, rmrat, rho
  
  real(kind=f)          :: wave(NWAVE)
  complex(kind=f)       :: refidx(NWAVE)
  real(kind=f)          :: r_refidx(NWAVE)
  real(kind=f)          :: i_refidx(NWAVE)
  real(kind=f)          :: qext(NWAVE, NBIN)
  real(kind=f)          :: ssa(NWAVE, NBIN)
  real(kind=f)          :: asym(NWAVE, NBIN)
  real(kind=f)          :: r(NBIN)
  
  ! Set the wavelengths
  data wave &
  /  0.340_f,  0.380_f,  0.412_f,  0.440_f,  0.443_f,  0.490_f, &
     0.500_f,  0.531_f,  0.532_f,  0.551_f,  0.555_f,  0.667_f, &
		 0.675_f,  0.870_f,  1.020_f,  1.640_f,  1.111_f,  1.333_f, &
		 1.562_f,  1.770_f,  2.051_f,  2.210_f,  2.584_f,  3.284_f, &
     3.809_f,  4.292_f,  4.546_f,  4.878_f,  5.128_f,  5.405_f, &
     5.714_f,  6.061_f,  6.452_f,  6.897_f,  7.407_f,  8.333_f, &
     9.009_f, 10.309_f, 12.500_f, 13.889_f, 16.667_f, 20.000_f, &
    26.316_f, 35.714_f /  
    
  data r_refidx &
  / 1.343_f, 1.341_f, 1.339_f, 1.337_f, 1.337_f, 1.335_f, &
    1.335_f, 1.334_f, 1.334_f, 1.333_f, 1.333_f, 1.331_f, &
    1.331_f, 1.329_f, 1.327_f, 1.317_f, 1.327_f, 1.323_f, &
    1.319_f, 1.313_f, 1.305_f, 1.295_f, 1.252_f, 1.455_f, &
    1.362_f, 1.334_f, 1.326_f, 1.320_f, 1.308_f, 1.283_f, &
    1.278_f, 1.313_f, 1.326_f, 1.310_f, 1.293_f, 1.270_f, &
    1.227_f, 1.164_f, 1.173_f, 1.287_f, 1.415_f, 1.508_f, &
    1.541_f, 1.669_f /    
      
  data i_refidx &
  /  6.5e-9_f,   4.0e-9_f,   1.86e-9_f,  1.02e-9_f,  1.02e-9_f,  1.0e-9_f,  &
     1.0e-9_f,   1.5e-9_f,   1.5e-9_f,   1.96e-9_f,  1.96e-9_f,  3.35e-8_f, &
     3.35e-8_f,  2.93e-7_f,  2.89e-6_f,  8.55e-5_f,  2.05E-06_f, 2.39E-05,  &
     1.20E-04_f, 1.18E-04_f, 6.79E-04_f, 3.51E-04_f, 2.39E-03_f, 0.0442_f,  &
     0.00339_f,  0.00833_f,  0.0139_f,   0.0125_f,   0.011_f,    0.015_f,   &
     0.075_f,    0.086_f,    0.039_f,    0.035_f,    0.035_f,    0.038,     &
     0.051_f,    0.161_f,    0.308_f,    0.39_f,     0.42_f,     0.395_f,   &
     0.373_f,    0.5_f /
  

!  write(*,*) ""
!  write(*,*) "Mie Calculations"

  ! Open the output text file
  open(unit=lun,file="carma_mietest.txt",status="unknown")
  
  ! Convert wavelength um -> cm
  wave = wave * 1e-4_f
  
  ! Construct a complex refractive index.
  refidx = cmplx(r_refidx, i_refidx)
  
  
  ! Define the particle-grid extent of the CARMA test
!  write(*,*) "  CARMA_Create(carma, ", NBIN,    ", ", NELEM, ", ", NGROUP, &
!                                 ", ", NSOLUTE, ", ", NGAS, ", rc, 6) ..."
  call CARMA_Create(carma, NBIN, NELEM, NGROUP, NSOLUTE, NGAS, NWAVE, rc, LUNOPRT=6, wave=wave)
  if (rc < 0) stop "    *** FAILED ***"
	carma_ptr => carma


  ! Define the group
  rho = 2.65_f
  rmrat = 4.32_f
  rmin = 1e-6_f
!  write(*,*) "  Add Group(s) ..."
  call CARMAGROUP_Create(carma, 1, "dust", rmin, rmrat, &
                         I_SPHERE, 1._f, .FALSE., rc, refidx=refidx, do_mie=.true.)
  if (rc < 0) stop "    *** FAILED ***"

  
  ! Define the element
!  write(*,*) "  Add Element(s) ..."
  call CARMAELEMENT_Create(carma, 1, 1, "dust", rho, I_INVOLATILE, I_DUST, rc)
  if (rc < 0) stop "    *** FAILED ***"
  
  
  ! Setup the CARMA processes to exercise
!  write(*,*) "  Initialize ..."
  call CARMA_Initialize(carma, rc)
  if (rc < 0) stop "    *** FAILED ***"

  
  ! Print the Group Information
!  write(*,*)  ""
!  call dumpGroup(carma, rc)
!  if (rc < 0) stop "    *** FAILED ***"
  
  ! Print the Element Information
!  write(*,*)  ""
!  call dumpElement(carma, rc)
!  if (rc < 0) stop "    *** FAILED ***"

!  write(*,*) ""
  
  ! Write output for the falltest
  write(lun,*) NGROUP, NWAVE, NBIN
  
  do iwave = 1, NWAVE
   write(lun,'(i3,e10.3)') iwave, wave(iwave)
  end do
  
  do igroup = 1, NGROUP
    
    call CARMAGROUP_Get(carma, igroup, rc, r=r, qext=qext, ssa=ssa, asym=asym)
    
    do ibin = 1, NBIN
      write(lun,'(i3,e10.3)') ibin, r(ibin)
    end do
  
    do iwave = 1, NWAVE
      write(lun,'(i3,2e10.3)') iwave, refidx(iwave)
    end do

    do iwave = 1, NWAVE
      do ibin = 1, NBIN
        write(lun,'(2i3,3(x,e10.3))') iwave, ibin, qext(iwave,ibin), ssa(iwave,ibin), asym(iwave,ibin)
      end do
    end do
  end do

  ! Close the output file
  close(unit=lun)	
	
!  write(*,*)  ""
!  write(*,*) "  CARMA_Destroy() ..."
  call CARMA_Destroy(carma, rc)
  if (rc /=0) stop "    *** FAILED ***"
end subroutine

