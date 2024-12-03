!! This code is to test changes made to GetBin and SetBin to
!! remove the special behavior of the concentration element from
!! the caller to CARMASTATE. To the caller the concentration element
!! will just be the mass mixing ratio of that element NOT the mass
!! mixing ratio of the entire particle. That representation will
!! still be used internally within CARMA. This avoids errors that can
!! be introduced by advection of CARMA particles.
!!
!! Upon execution, a text file (carma_cetest.txt) is generated.
!! The text file can be read with the IDL procedure read_cetest.pro.
!!
!! @author  Chuck Bardeen
!! @version Feb-2024

program carma_cetest
  implicit none

  write(*,*) "Concentration Element Test"

  call test_conc_elem()

  write(*,*) "Done"
end program

!! Just have one grid box. We don't actually need to step the timestep,
!! we just want to do CARMASTATE_SetBin followed by CARMASTATE_GetBin
!! and confirm that we get the same value back.
subroutine test_conc_elem()
  use carma_precision_mod
  use carma_constants_mod
  use carma_enums_mod
  use carma_types_mod
  use carmaelement_mod
  use carmagroup_mod
  use carmagas_mod
  use carmasolute_mod
  use carmastate_mod
  use carma_mod
  use atmosphere_mod

  implicit none

  integer, parameter        :: NX           = 1
  integer, parameter        :: NY           = 1
  integer, parameter        :: NZ           = 1
  integer, parameter        :: NZP1         = NZ+1
  integer, parameter        :: NELEM        = 3
  integer, parameter        :: NBIN         = 16
  integer, parameter        :: NGROUP       = 2
  integer, parameter        :: NSOLUTE      = 0
  integer, parameter        :: NGAS         = 0
  integer, parameter        :: NWAVE        = 0
  integer, parameter        :: LUNOPRT      = 6
  integer, parameter        :: nstep        = 1


  real(kind=f), parameter   :: dtime  = 1._f
  real(kind=f), parameter   :: deltaz = 100._f
  real(kind=f), parameter   :: zmin   = 3000._f

  real(kind=f), parameter   :: n    = 100._f     !! concentration (cm-3)
  real(kind=f), parameter   :: r0   = 2.5e-6_f   !! mean radius (cm)
  real(kind=f), parameter   :: rsig = 1.5_f      !! distribution width

  real(kind=f)              :: rhop
  real(kind=f)              :: rhoa

  integer, parameter        :: I_H2SO4   = 1      !! sulfate aerosol composition
  integer, parameter        :: I_ICE     = 2      !! ice
  integer, parameter        :: I_WATER   = 3      !! water

  integer, parameter        :: I_GRP_IN  = 1      !! sulfate aerosol composition
  integer, parameter        :: I_GRP_ICE = 2      !! ice

  integer, parameter        :: I_ELEM_IN  = 1      !! sulfate aerosol composition
  integer, parameter        :: I_ELEM_ICE = 2      !! ice
  integer, parameter        :: I_ELEM_CM  = 3      !! sulfate core mass

  type(carma_type), target  :: carma
  type(carma_type), pointer :: carma_ptr
  type(carmastate_type)     :: cstate
  integer                   :: rc = 0

  real(kind=f), allocatable   :: zc(:,:,:)
  real(kind=f), allocatable   :: zl(:,:,:)
  real(kind=f), allocatable   :: p(:,:,:)
  real(kind=f), allocatable   :: pl(:,:,:)
  real(kind=f), allocatable   :: t(:,:,:)
  real(kind=f), allocatable   :: relhum(:,:,:)
  real(kind=f), allocatable   :: rho(:,:,:)

  real(kind=f), allocatable   :: mmr(:,:,:,:,:)
  real(kind=f), allocatable   :: mmr2(:,:,:,:,:)

  real(kind=f), allocatable   :: r(:)
  real(kind=f), allocatable   :: dr(:)
  real(kind=f), allocatable   :: rmass(:)

  real(kind=f), allocatable   :: lat(:,:)
  real(kind=f), allocatable   :: lon(:,:)

  integer               :: i
  integer               :: ix
  integer               :: iy
  integer               :: ixy
  integer               :: iz
  integer               :: istep = 1
  integer               :: igas
  integer               :: igroup
  integer               :: ielem
  integer               :: ibin
  integer               :: ithread
  integer, parameter    :: lun = 42

  real(kind=f)          :: time = 0.
  real(kind=f)          :: rmin, rmrat
  real(kind=f)          :: drh

  real(kind=f)          :: RHO_CN = 1.78_f


!  write(*,*) ""
!  write(*,*) "Particle Growth - Simple"

  ! Open the output text file
  open(unit=lun,file="carma_cetest.txt",status="unknown")

  ! Allocate the arrays that we need for the model
  allocate(zc(NZ,NY,NX), zl(NZP1,NY,NX), p(NZ,NY,NX), pl(NZP1,NY,NX), &
           t(NZ,NY,NX), rho(NZ,NY,NX))
  allocate(mmr(NZ,NY,NX,NELEM,NBIN), mmr2(NZ,NY,NX,NELEM,NBIN))
  allocate(r(NBIN))
  allocate(dr(NBIN))
  allocate(rmass(NBIN))
  allocate(lat(NY,NX), lon(NY,NX))


  ! Define the particle-grid extent of the CARMA test
!  write(*,*) "  CARMA_Create ..."
  call CARMA_Create(carma, NBIN, NELEM, NGROUP, NSOLUTE, NGAS, NWAVE, rc, LUNOPRT=LUNOPRT)
  if (rc /=0) stop "    *** FAILED ***"
  carma_ptr => carma


  ! Define the groups
!  write(*,*) "  Add Group(s) ..."
  call CARMAGROUP_Create(carma, I_GRP_IN, "Sulfate IN", 1.e-7_f, 4._f, I_SPHERE, 1._f, .false., &
                         rc, do_wetdep=.false., do_drydep=.false., solfac=0.3_f, &
                         scavcoef=0.1_f, shortname="CRIN", do_mie=.false.)
  if (rc /=0) stop "    *** FAILED ***"

  call CARMAGROUP_Create(carma, I_GRP_ICE, "Ice Crystal", 5.e-5_f, 4.0_f, I_SPHERE, 3._f, .true., &
                         rc, do_wetdep=.false., do_drydep=.false., solfac=0.3_f, &
                         scavcoef=0.1_f, shortname="CRICE", do_mie=.false.)
  if (rc /=0) stop "    *** FAILED ***"


  ! Define the elements
!  write(*,*) "  Add Element(s) ..."
  call CARMAELEMENT_Create(carma, 1, I_GRP_IN, "Sulfate IN", RHO_CN, I_INVOLATILE, I_H2SO4, rc, shortname="CRIN")
  if (rc /=0) stop "    *** FAILED ***"

  call CARMAELEMENT_Create(carma, 2, I_GRP_ICE, "Ice Crystal", RHO_I, I_VOLATILE, I_ICE, rc, shortname="CRICE")
  if (rc /=0) stop "    *** FAILED ***"

  call CARMAELEMENT_Create(carma, 3, I_GRP_ICE, "Core Mass", RHO_CN, I_COREMASS, I_H2SO4, rc, shortname="CRCORE")
  if (rc /=0) stop "    *** FAILED ***"


  ! Define the Solutes


  ! Define the gases


!  write(*,*) "  Initialize ..."
  call CARMA_Initialize(carma, rc)
  if (rc /=0) stop "    *** FAILED ***"


  ! Print the Group Information
!  write(*,*)  ""
!  call dumpGroup(carma, rc)
!  if (rc /=0) stop "    *** FAILED ***"

  ! Print the Element Information
!  write(*,*)  ""
!  call dumpElement(carma, rc)
!  if (rc /=0) stop "    *** FAILED ***"

  ! Print the Gas Information
!  write(*,*)  ""
!  call dumpGas(carma, rc)
!  if (rc /=0) stop "    *** FAILED ***"

!  write(*,*) ""


  ! For simplicity of setup, do a case with Cartesian coordinates,
  ! which are specified in this interface in meters.
  !
  ! NOTE: For Cartesian coordinates, the first level is the bottom
  ! of the model (e.g. z = 0), while for sigma and hybrid coordinates
  ! the first level is the top of the model.
  lat(:,:) = -40.0_f
  lon(:,:) = -105.0_f

  ! Vertical center
  do i = 1, NZ
    zc(i,:,:) = zmin + (deltaz * (i - 0.5_f))
  end do

  call GetStandardAtmosphere(zc, p=p, t=t)

  ! Vertical edge
  do i = 1, NZP1
    zl(i,:,:) = zmin + ((i - 1) * deltaz)
  end do
  call GetStandardAtmosphere(zl, p=pl)


  ! Write output for the test
  write(lun,*) NGROUP, NELEM, NBIN, NGAS

  do igroup = 1, NGROUP
    call CARMAGROUP_Get(carma, igroup, rc, r=r, dr=dr, rmass=rmass)
    if (rc /=0) stop "    *** FAILED ***"

    do ibin = 1, NBIN
      write(lun,'(2i4,2e10.3)') igroup, ibin, r(ibin) * 1e4_f, rmass(ibin)
    end do
  end do


  ! Try TTL Conditions ...
  !
  ! p, T, z, mmrgas, rmin, particle concentration
  ! 90 hPa, 190 K, 17 km, 3.5e-6 g/g, 5 um, 0.1/cm^3.
  p(1,:,:)         = 90._f * 100._f
  zc(1,:,:)        = 17000._f
!  t(1,:,:)         = 190._f
  t(1,:,:)         = 205._f
  zl(1,:,:)        = zc(1,:,:) - deltaz
  zl(2,:,:)        = zc(1,:,:) + deltaz
  rho(1,:,:)       = (p(1,:,:) * 10._f) / (R_AIR * t(1,:,:)) * (1e-3_f * 1e6_f)
  pl(1,:,:)        = p(1,:,:) - (zl(1,:,:) - zc(1,:,:)) * rho(1,:,:) * (GRAV / 100._f)
  pl(2,:,:)        = p(1,:,:) - (zl(2,:,:) - zc(1,:,:)) * rho(1,:,:) * (GRAV / 100._f)

  ! Initialize the mmr to 0
  mmr(:,:,:,:,:) = 0._f


  ! Put in an intial distribution of sulfates in the IN.
  call CARMAGROUP_Get(carma, I_GRP_IN, rc, r=r, dr=dr, rmass=rmass)
  if (rc /=0) stop "    *** FAILED ***"
  do ibin = 1, NBIN
    rhop = (100._f * dr(ibin) / (sqrt(2._f*PI) * r(ibin) * log(rsig))) * exp(-((log(r(ibin)) - log(r0))**2) / (2._f*(log(rsig))**2)) * rmass(ibin)

    ! We don't know rhoa for the initial condition, but assume something typical of
    ! the conditions at 100 mb and 200K. (mb -> dynes, since R_AIR in cgs)
    rhoa = 100._f * 1000._f / (R_AIR) / (200._f)

    mmr(:,:,:,I_ELEM_IN,ibin)   = rhop / rhoa
  end do

  ! Put in an intial distribution of sulfates and wate in the ICE.
  call CARMAGROUP_Get(carma, I_GRP_IN, rc, r=r, dr=dr, rmass=rmass)
  if (rc /=0) stop "    *** FAILED ***"
  do ibin = 1, NBIN
    rhop = (10._f * dr(ibin) / (sqrt(2._f*PI) * r(ibin) * log(rsig))) * exp(-((log(r(ibin)) - log(r0))**2) / (2._f*(log(rsig))**2)) * rmass(ibin)

    ! We don't know rhoa for the initial condition, but assume something typical of
    ! the conditions at 100 mb and 200K. (mb -> dynes, since R_AIR in cgs)
    rhoa = 100._f * 1000._f / (R_AIR) / (200._f)

    mmr(:,:,:,I_ELEM_CM,ibin)   = rhop / rhoa
    mmr(:,:,:,I_ELEM_ICE,ibin)  = 2. * rhop / rhoa
  end do

  write(lun,*) 0

  do ielem = 1, NELEM
    do ibin = 1, NBIN
      write(lun,'(2i4,e10.3)') ielem, ibin, real(mmr(1,NY,NX,ielem,ibin))
    end do
  end do

  do ixy = 1, NX*NY
    ix = ((ixy-1) / NY) + 1
    iy = ixy - (ix-1)*NY

    ! Create a CARMASTATE for this column.
    call CARMASTATE_Create(cstate, carma_ptr, time, dtime, NZ, &
                      I_CART, lat(iy,ix), lon(iy,ix), &
                      zc(:,iy,ix), zl(:,iy,ix), p(:,iy,ix), &
                      pl(:,iy,ix), t(:,iy,ix), rc)
    if (rc /=0) stop "    *** FAILED ***"

    ! Send the bin mmrs to CARMA
    do ielem = 1, NELEM
      do ibin = 1, NBIN
        call CARMASTATE_SetBin(cstate, ielem, ibin, &
                            mmr(:,iy,ix,ielem,ibin), rc)
        if (rc /=0) stop "    *** FAILED ***"
      end do
    end do

    ! Get the mmr back from state (should be the same).
    do ielem = 1, NELEM
      do ibin = 1, NBIN
        call CARMASTATE_GetBin(cstate, ielem, ibin, &
                            mmr2(:,iy,ix,ielem,ibin), rc)
        if (rc /=0) stop "    *** FAILED ***"
      end do
    end do

    ! The value from mmr2 should be the same as mmr.
    do ielem = 1, NELEM
      do ibin = 1, NBIN
        do iz = 1, NZ
          if ((abs((mmr2(iz,iy,ix,ielem,ibin) - mmr(iz,iy,ix,ielem,ibin)) / mmr(iz,iy,ix,ielem,ibin))) .gt. 1e-15) then
            write(LUNOPRT,*) ielem, ibin, mmr(:,iy,ix,ielem,ibin), mmr2(:,iy,ix,ielem,ibin), mmr2(:,iy,ix,ielem,ibin) - mmr(:,iy,ix,ielem,ibin)
!            stop "    *** FAILED GetBin did not return the same valued as SetBin***"
          end if
        end do
      end do
    end do

    ! Write output for the concentration element test
    write(lun,'(f12.0)') istep*dtime
    do ielem = 1, NELEM
      do ibin = 1, NBIN
        write(lun,'(2i4,e12.3)') ielem, ibin, real(mmr(1,NY,NX,ielem,ibin))
      end do
    end do
  end do

  ! Cleanup the carma state objects
  call CARMASTATE_Destroy(cstate, rc)
  if (rc /=0) stop "    *** FAILED ***"

  ! Close the output file
  close(unit=lun)

!  write(*,*)  ""

  if (rc /=0) stop "    *** FAILED ***"

!  write(*,*)  ""
!  write(*,*) "  CARMA_Destroy() ..."
  call CARMA_Destroy(carma, rc)
  if (rc /=0) stop "    *** FAILED ***"
end subroutine
