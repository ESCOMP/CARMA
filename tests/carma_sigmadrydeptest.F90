!! This code is to test the dry deposition routine by comparing
!! sedimentation with and without dry deposition. The model is
!! using sigma coordinates.
!!
!! Upon execution, a text file (carma_drydeptest.txt) is generated.
!! The text file can be read with the IDL procedure read_drydeptest.pro.
!!
!! @author  Tianyi Fan
!! @version Apr-2011


program carma_sigmadrydeptest
  implicit none

  write(*,*) "Sedimentation Test (Sigma Coordinates)"

  call test_sigmadrydep()  
  
  write(*,*) "Done"
end program


!! Create 2 particle groups, one for particles with dry deposition 
!! using arbitary values of ram(aerodynamic resistance) and fv (friction velocity)
!! the other for particles without dry deposition, to see if they make a difference.
subroutine test_sigmadrydep()
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

  integer, parameter    :: NZ           = 72
  integer, parameter    :: NZP1         = NZ+1
  integer, parameter    :: NELEM        = 2
  integer, parameter    :: NBIN         = 16
  integer, parameter    :: NGROUP       = 2
  integer, parameter    :: NSOLUTE      = 0
  integer, parameter    :: NGAS         = 0
  integer, parameter    :: NWAVE        = 0
  integer, parameter    :: nstep        = 100*6
  
  
  ! To keep the file processing simpler, only one bin will get written out
  ! to the output file. 
  integer, parameter        :: OUTBIN       = 14
  
  real(kind=f), parameter   :: dtime  = 1000._f
  real(kind=f), parameter   :: deltax = 100._f
  real(kind=f), parameter   :: deltay = 100._f
  
  integer, parameter        :: I_PART       = 1

  type(carma_type), target  :: carma
  type(carma_type), pointer :: carma_ptr
  type(carmastate_type)     :: cstate
  integer                   :: rc = 0
  
  real(kind=f), allocatable   :: xc(:)
  real(kind=f), allocatable   :: dx(:)
  real(kind=f), allocatable   :: yc(:)
  real(kind=f), allocatable   :: dy(:)
  real(kind=f), allocatable   :: zc(:)
  real(kind=f), allocatable   :: zl(:)
  real(kind=f), allocatable   :: p(:)
  real(kind=f), allocatable   :: pl(:)
  real(kind=f), allocatable   :: t(:)
  
  real(kind=f), allocatable, target  :: mmr(:,:,:)
  
  real(kind=f)          :: lat
  real(kind=f)          :: lon
  
  integer               :: i
  integer               :: istep
  integer               :: ielem
  integer               :: ibin
  integer               :: igroup
  integer, parameter    :: lun = 42
  integer, parameter    :: lun1 = 41

  real(kind=f)          :: time
  real(kind=f)          :: rmin, rmrat, rho
  
  real(kind=f)          :: a72(73), b72(73), t72(72), ze(73)
  real(kind=f)          :: dz(NZ), zm(NZ), rhoa(NZ)
  real(kind=f), parameter :: ps = 98139.8  ! GEOS-5, Omaha, NE, 20090101

  real(kind=f)          :: lndfv    = 1.5_f   ! land friction velocity
  real(kind=f)          :: lndram   = 60._f   ! land aerodynamic resistance
  real(kind=f)          :: lndfrac  = 0.0_f   ! land fraction
  
  real(kind=f)          :: ocnfv    = 2.0_f   ! ocean friction velocity
  real(kind=f)          :: ocnram   = 40._f   ! ocean aerodynamic resistance
  real(kind=f)          :: ocnfrac  = 1.0_f   ! ocean fraction
  
  real(kind=f)          :: icefv    = 2.5_f   ! ice friction velocity
  real(kind=f)          :: iceram   = 20._f   ! ice aerodynamic resistance
  real(kind=f)          :: icefrac  = 0.0_f   ! ice fraction

  real(kind=f)          :: vdry(NBIN, NGROUP)
  
  data a72 / &
       1.0000000,       2.0000002,       3.2700005,       4.7585009,       6.6000011, &
       8.9345014,       11.970302,       15.949503,       21.134903,       27.852606, &
       36.504108,       47.580610,       61.677911,       79.513413,       101.94402, &
       130.05102,       165.07903,       208.49704,       262.02105,       327.64307, &
       407.65710,       504.68010,       621.68012,       761.98417,       929.29420, &
       1127.6902,       1364.3402,       1645.7103,       1979.1604,       2373.0405, &
       2836.7806,       3381.0007,       4017.5409,       4764.3911,       5638.7912, &
       6660.3412,       7851.2316,       9236.5722,       10866.302,       12783.703, &
       15039.303,       17693.003,       20119.201,       21686.501,       22436.301, &
       22389.800,       21877.598,       21214.998,       20325.898,       19309.696, &
       18161.897,       16960.896,       15625.996,       14290.995,       12869.594, &
       11895.862,       10918.171,       9936.5219,       8909.9925,       7883.4220, &
       7062.1982,       6436.2637,       5805.3211,       5169.6110,       4533.9010, &
       3898.2009,       3257.0809,       2609.2006,       1961.3106,       1313.4804, &
       659.37527,       4.8048257,       0.0000000 /

  data b72 / &
       0.0000000,       0.0000000,       0.0000000,       0.0000000,       0.0000000, &
       0.0000000,       0.0000000,       0.0000000,       0.0000000,       0.0000000, &
       0.0000000,       0.0000000,       0.0000000,       0.0000000,       0.0000000, &
       0.0000000,       0.0000000,       0.0000000,       0.0000000,       0.0000000, &
       0.0000000,       0.0000000,       0.0000000,       0.0000000,       0.0000000, &
       0.0000000,       0.0000000,       0.0000000,       0.0000000,       0.0000000, &
       0.0000000,       0.0000000,       0.0000000,       0.0000000,       0.0000000, &
       0.0000000,       0.0000000,       0.0000000,       0.0000000,       0.0000000, &
       0.0000000,   8.1754130e-09,    0.0069600246,     0.028010041,     0.063720063, &
      0.11360208,      0.15622409,      0.20035011,      0.24674112,      0.29440312, &
      0.34338113,      0.39289115,      0.44374018,      0.49459020,      0.54630418, &
      0.58104151,      0.61581843,      0.65063492,      0.68589990,      0.72116594, &
      0.74937819,      0.77063753,      0.79194696,      0.81330397,      0.83466097, &
      0.85601798,      0.87742898,      0.89890800,      0.92038701,      0.94186501, &
      0.96340602,      0.98495195,       1.0000000 /
      
  data t72 / &
      212.161,   210.233,   217.671,   225.546,   232.222,   237.921,   241.836,   243.246, &
      250.032,   265.518,   262.335,   255.389,   253.560,   253.848,   252.496,   247.806, &
      243.108,   237.288,   230.839,   226.233,   221.617,   218.474,   218.014,   218.881, &
      220.297,   222.262,   224.564,   224.059,   221.671,   220.732,   220.200,   218.445, &
      217.424,   215.322,   212.882,   211.080,   210.573,   210.942,   212.593,   214.064, &
      213.704,   209.045,   211.286,   218.995,   227.209,   235.050,   241.144,   246.328, &
      250.606,   254.079,   257.222,   260.012,   262.534,   265.385,   267.348,   267.998, &
      267.964,   267.827,   268.075,   268.397,   268.440,   268.371,   268.302,   268.203, &
      267.943,   266.305,   265.331,   265.628,   266.371,   267.219,   267.981,   268.379 /
  
  data ze / &
      78126.3,   73819.6,   70792.7,   68401.3,   66240.5,   64180.9,   62142.8,   60110.2,   58104.9,   56084.0,  53980.6, &
      51944.7,   50003.9,   48117.8,   46270.4,   44469.8,   42739.0,   41076.6,   39488.8,   37977.8,   36530.2,  35144.5, &
      33810.5,   32511.3,   31238.9,   29990.5,   28750.5,   27517.4,   26306.8,   25128.6,   23974.7,   22842.9,  21739.4, &
      20653.8,   19591.3,   18553.2,   17536.4,   16534.3,   15530.5,   14518.7,   13500.1,   12483.0,   11491.9,  10495.9, &
      9466.47,   8427.35,   7712.39,   7048.43,   6429.18,   5849.30,   5304.75,   4791.16,   4305.46,   3844.44,  3404.89, &
      3122.98,   2850.14,   2586.45,   2331.56,   2084.55,   1892.20,   1750.98,   1612.29,   1476.05,   1342.20,  1210.71, &
      1082.17,   956.177,   832.052,   709.535,   588.531,   469.023,   350.157 /


  ! Open the output text file
  open(unit=lun,file="carma_sigmadrydeptest.txt",status="unknown")
  
  ! Allocate the arrays that we need for the model
  allocate(xc(NZ), dx(NZ), yc(NZ), dy(NZ), &
           zc(NZ), zl(NZP1), p(NZ), pl(NZP1), &
           t(NZ))
  allocate(mmr(NZ,NELEM,NBIN))
  
  
  ! Define the particle-grid extent of the CARMA test
  call CARMA_Create(carma, NBIN, NELEM, NGROUP, NSOLUTE, NGAS, NWAVE, rc, LUNOPRT=6)
  if (rc /=0) write(*, *) "    *** CARMA_Create FAILED ***, rc=", rc
  carma_ptr => carma

  ! Define the group
  rho = 2.65_f
  rmrat = 4.32_f
  rmin = 1e-6_f

  call CARMAGROUP_Create(carma, 1, "DryDep", rmin, rmrat, I_SPHERE, 1._f, .FALSE., rc, &
       do_mie = .FALSE., do_wetdep=.FALSE., do_drydep=.TRUE., do_vtran=.TRUE., shortname="DD") 
  if (rc /=0) stop "    *** CARMAGROUP_Create FAILED ***"
  
  call CARMAGROUP_Create(carma, 2, "NoDryDep", rmin, rmrat, I_SPHERE, 1._f, .FALSE., rc, &
       do_mie = .FALSE., do_wetdep=.FALSE., do_drydep=.FALSE., do_vtran=.TRUE., shortname="NDD")
  if (rc /=0) stop "    *** CARMAGROUP_Create FAILED ***"
  
  ! Define the element
  call CARMAELEMENT_Create(carma, 1, 1, "DryDep", rho, I_INVOLATILE, I_PART, rc)
  if (rc /=0) stop "    *** CARMAELEMENT_Create FAILED ***"
  
  call CARMAELEMENT_Create(carma, 2, 2, "NoDryDep", rho, I_INVOLATILE, I_PART, rc)
  if (rc /=0) stop "    *** CARMAELEMENT_Create FAILED ***"
  
!  write(*,*) "  CARMA_AddCoagulation(carma, 1, 1, 1, I_COLLEC_DATA, rc) ..."
!  call CARMA_AddCoagulation(carma, 1, 1, 1, I_COLLEC_DATA, rc) 
!  if (rc /=0) write(*, *) "    *** FAILED ***, rc=", rc
   
  ! Setup the CARMA processes to exercise
  call CARMA_Initialize(carma, rc, do_vtran=.TRUE., do_drydep=.TRUE.)
  if (rc /=0) write(*, *) "    *** CARMA_Initialize FAILED ***, rc=", rc
  

  ! For simplicity of setup, do a case with Cartesian coordinates,
  ! which are specified in this interface in meters.
  !
  ! NOTE: For Cartesian coordinates, the first level is the bottom 
  ! of the model (e.g. z = 0), while for sigma and hybrid coordinates
  ! the first level is the top of the model.
  lat = 40.0_f
  lon = -105.0_f
  
  ! Horizonal centers
  dx(:) = deltax
  xc(:) = dx(:) / 2._f
  dy(:) = deltay
  yc(:) = dy(:) / 2._f

  ! Layer edges
  pl(:) = a72(:) + b72(:)*ps
  zl(:) = a72(:) / 1e5_f + b72(:)
  t(:)  = t72(:)


  ! Calculate based upon the known fields (edges to middle, ...)  
  p(:)    = (pl(1:NZ) + pl(2:NZP1)) / 2._f
  zc(:)   = (zl(1:NZ) + zl(2:NZP1)) / 2._f
  
  rhoa(:) = p(:) / 287._f / t(:)
  dz(:)   = zl(2:NZP1) - zl(1:NZ)
  zm(:) = (ze(2:NZP1) + ze(1:NZ)) / 2._f
  dz(:)   = ze(1:NZ) - ze(2:NZP1)

  			
  ! Put a blob in the model for all elements and bins at 8 km.
  mmr(:,:,:) = 0._f
  do ielem = 1, NELEM
    do ibin = 1, NBIN        
      mmr(:,ielem,ibin) = 1e-10_f * exp(-((zm(:) - 8.e3_f) / 3.e3_f)**2) / rhoa(:)
    end do
  end do		     

   
  ! Write output for the test
  write(lun,*) NZ, NELEM
  do i = 1, NZ
    write(lun,'(i3,2f10.1)') i, zm(i), dz(i)
  end do
  
  write(lun,*) 0
  do ielem = 1, NELEM
    do i = 1, NZ
      write(lun,'(2i4,e10.3,e10.3)') ielem, i, real(mmr(i,ielem,OUTBIN)), real(mmr(i,ielem,OUTBIN)*rhoa(i))
    end do
  end do
  
  
  ! Iterate the model over a few time steps.
  do istep = 1, nstep
  
    ! Calculate the model time.
    time = (istep - 1) * dtime

    ! Create a CARMASTATE for this column.
    call CARMASTATE_Create(cstate, carma_ptr, time, dtime, NZ, &
                              I_HYBRID, I_CART, lat, lon, &
                              xc(:), dx(:), &
                              yc(:), dy(:), &
                              zc(:), zl(:), &
			      p(:), pl(:), &
			      t(:), rc)
    if (rc /=0) stop "    *** CARMASTATE_Create FAILED ***"
		
    ! Send the bin mmrs to CARMA
    do ielem = 1, NELEM
      do ibin = 1, NBIN
        call CARMASTATE_SetBin(cstate, ielem, ibin, mmr(:,ielem,ibin), rc)
        if (rc /=0) stop "    *** CARMASTATE_SetBin FAILED ***"
      end do
    end do
			
    ! Execute the step
    call CARMASTATE_Step(cstate, rc, &
      lndfv=lndfv,     ocnfv=ocnfv,     icefv=icefv, &
      lndram=lndram,   ocnram=ocnram,   iceram=iceram, &
      lndfrac=lndfrac, ocnfrac=ocnfrac, icefrac=icefrac)
    if (rc /=0) stop "    *** CARMASTATE_Step FAILED ***"
	
    ! Get the updated bin mmr.
    do ielem = 1, NELEM
      do ibin = 1, NBIN
        call CARMASTATE_GetBin(cstate, ielem, ibin, mmr(:,ielem,ibin), rc, vd=vdry(ibin,ielem))
        if (rc /=0) stop "    *** CARMASTATE_GetBin FAILED ***"
      end do
    end do

    ! Get the updated temperature.
    call CARMASTATE_GetState(cstate, rc, t=t(:))
    if (rc /=0) stop "    *** CARMASTATE_GetState FAILED ***"
 
    ! Write output for the test
    write(lun,'(f12.0)') istep*dtime
    
    do ielem = 1, NELEM
      do i = 1, NZ
        write(lun,'(2i4,e10.3,e10.3)') ielem, i, real(mmr(i,ielem,OUTBIN)), real(mmr(i,ielem,OUTBIN)*p(i) / 287._f / t(i))
      end do
    end do
  end do   ! time loop
  
  
  ! Cleanup the carma state objects
  call CARMASTATE_Destroy(cstate, rc)
  if (rc /=0) stop "    *** CARMASTATE_Destroy FAILED ***"

  ! Close the output file
  close(unit=lun)	
  
  
  ! write the dry deposition velocity
  open(unit=lun1,file="carma_sigmavdry.txt",status="unknown")
  
  write(lun1, *) NGROUP
  do igroup = 1, NGROUP
    write(lun1,*) igroup, real(vdry(:, igroup))     
  end do

  close(unit=lun1)
  	
  call CARMA_Destroy(carma, rc)
  if (rc /=0) write(*, *) "    *** CARMA_Destroy FAILED ***, rc=", rc  
end subroutine

