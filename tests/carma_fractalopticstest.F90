!! This code is to demonstrate the CARMA optical treatment of
!! fractal aggregates composed of identical spheres.
!!
!! Upon execution, a text file (carma_fractalopticstest.txt) is generated.
!! The text file can be read with the IDL procedure read_fractalopticstest.pro.
!!
!! @author  Eric Wolf (based on Chuck Bardeen's code)
!! @version March-2013

program carma_fractalopticstest
  implicit none

  write(*,*) "Fractal Optics Test"

  call test_fractaloptics()

  write(*,*) "Done"
end program


subroutine test_fractaloptics()
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
  integer, parameter    :: NBIN         = 1
  integer, parameter    :: NGROUP       = 1
  integer, parameter    :: NSOLUTE      = 0
  integer, parameter    :: NGAS         = 0
  integer, parameter    :: NWAVE        = 10
  integer, parameter    :: NRUN		= 2
  integer, parameter    :: NREFIDX      = 1  !! Number of refractive indices per element

  integer, parameter    :: I_DUST       = 1

  type(carma_type), target            :: carma
  type(carma_type), pointer           :: carma_ptr
  integer                             :: rc = 0


  integer               :: i
  integer               :: ielem
  integer               :: iwave
  integer               :: ibin
  integer               :: igroup
  integer               :: irun
  integer, parameter    :: lun = 42

  real(kind=f)          :: rmin, rmrat, rho, rmon, nmon, falpha

  real(kind=f)          :: wave(NWAVE)
  complex(kind=f)       :: refidx(NWAVE, NREFIDX)
  real(kind=f)          :: r_refidx(NWAVE)
  real(kind=f)          :: i_refidx(NWAVE)
  real(kind=f)          :: qext_data(NRUN, NWAVE, NBIN)
  real(kind=f)          :: ssa_data(NRUN, NWAVE, NBIN)
  real(kind=f)          :: asym_data(NRUN, NWAVE, NBIN)
  real(kind=f)          :: qext(NWAVE, NBIN)
  real(kind=f)          :: ssa(NWAVE, NBIN)
  real(kind=f)          :: asym(NWAVE, NBIN)
  real(kind=f)          :: r(NBIN)
  real(kind=f)		:: df(NBIN)

  data wave     /0.1181_f, 0.1362_f, 0.1968_f, 0.2952_f, 0.4133_f, 0.5635_f, 0.6888_f, 0.8731_f, 1.016_f, 2.019_f/
  data r_refidx /1.75_f,   1.70_f,   1.66_f,   1.66_f,   1.69_f,   1.70_f,   1.68_f,   1.66_f,   1.65_f,  1.63_f/
  data i_refidx /0.40_f,   0.27_f,   0.22_f,   0.15_f,   0.076_f,  0.023_f,  0.0088_f, 0.0024_f, 0.001_f, 0.00072_f/

!  write(*,*) ""
!  write(*,*) "Calculations"

  ! Open the output text file
  open(unit=lun,file="carma_fractalopticstest.txt",status="unknown")

  ! Convert wavelength um -> cm
  wave = wave * 1e-4_f

  ! Construct a complex refractive index.
  refidx(:,1) = cmplx(r_refidx(:), i_refidx(:), kind=f)


  do irun=1,NRUN

    ! Define the particle-grid extent of the CARMA test
    call CARMA_Create(carma, NBIN, NELEM, NGROUP, NSOLUTE, NGAS, NWAVE, rc, LUNOPRT=6, wave=wave, NREFIDX=NREFIDX)
    if (rc < 0) stop "    *** CARMA_Create FAILED ***"
    carma_ptr => carma

     ! Define the groups
     rho = 2.65_f   ! [g cm-3]
     rmrat = 2.5_f
     df(:)  = 2.0_f
     rmon = 50.0e-7_f  ! [cm]
     nmon = 100._f
     rmin = nmon**(1._f/3._f)*rmon
     falpha = 1.1_f

    if (irun .EQ. 1) then  ! do fractal optics calulcation
      call CARMAGROUP_Create(carma, 1, "fractal", rmin, rmrat, &
                             I_SPHERE, 1._f, .FALSE., rc, is_fractal=.TRUE., &
                             do_mie=.true., imiertn=I_MIERTN_BOTET1997, &
                             rmon=rmon, df=df, falpha=falpha)

      if (rc < 0) stop "    *** irun 1 CARMAGROUP_Create FAILED ***"
      ! Define the element
      call CARMAELEMENT_Create(carma, 1, 1, "dust", rho, I_INVOLATILE, I_DUST, rc, refidx=refidx )
      if (rc < 0) stop "    *** irun 1 CARMAELEMENT_Create FAILED ***"
    end if

    if (irun .EQ. 2) then  ! do standard mie calculation for comparison
       call CARMAGROUP_Create(carma, 1, "sphere", rmin, rmrat, &
                             I_SPHERE, 1._f, .FALSE., rc, is_fractal=.FALSE., &
                             do_mie=.true., imiertn=I_MIERTN_TOON1981)
       if (rc < 0) stop "    ***  irun 2 CARMAGROUP_CreateFAILED ***"
       ! Define the element
       call CARMAELEMENT_Create(carma, 1, 1, "dust", rho, I_INVOLATILE, I_DUST, rc, refidx=refidx )
       if (rc < 0) stop "    *** irun 2 CARMAELEMENT_Create FAILED ***"
    endif


    ! Setup the CARMA processes to exercise
    call CARMA_Initialize(carma, rc, do_pheat=.TRUE.)
    if (rc < 0) stop "    *** CARMA_Initialize FAILED ***"
    do igroup = 1, NGROUP
      call CARMAGROUP_Get(carma, igroup, rc, r=r, qext=qext, ssa=ssa, asym=asym)
      do iwave = 1, NWAVE
        do ibin = 1, NBIN
          qext_data(irun,iwave,ibin) = qext(iwave,ibin)
          ssa_data(irun,iwave,ibin) = ssa(iwave,ibin)
          asym_data(irun,iwave,ibin) = asym(iwave,ibin)
        end do
      end do
    end do


    !  write(*,*)  ""
    !  write(*,*) "  CARMA_Destroy() ..."
    call CARMA_Destroy(carma, rc)
    if (rc /=0) stop "    *** CARMA_Destroy FAILED ***"
  end do  ! irun

  ! Write grid parameters
  write(lun,*) NGROUP, NWAVE, NBIN

  do iwave = 1, NWAVE
    write(lun,'(i3,e10.3)') iwave, wave(iwave)
  end do

  do ibin = 1, NBIN
    write(lun,'(i3,e10.3)') ibin, r(ibin)
  end do

  do iwave = 1, NWAVE
    write(lun,'(i3,2e10.3)') iwave, refidx(iwave,1)
  end do

  do iwave = 1, NWAVE
    do ibin = 1, NBIN
      write(lun,'(2i3,6(x,e10.3))') iwave, ibin, qext_data(1,iwave,ibin), ssa_data(1,iwave,ibin), asym_data(1,iwave,ibin), &
                                                 qext_data(2,iwave,ibin), ssa_data(2,iwave,ibin), asym_data(2,iwave,ibin)
    end do
  end do

  ! Close the output file
  close(unit=lun)

end subroutine
