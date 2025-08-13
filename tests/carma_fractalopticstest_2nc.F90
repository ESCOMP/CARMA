!! This code is to demonstrate the CARMA optical treatment of
!! fractal aggregates composed of identical spheres.
!!
!! Upon execution, a netcdf file (carma_fractalopticstest.nc) is generated.
!! The text file can be plotted in ./tests/python/plots_optics.ipynb
!!
!! @author  Eric Wolf (based on Chuck Bardeen's code)
!! @version March-2013
!! Modified by Ilaria Quaglia June-2026

program carma_fractalopticstest
  implicit none

  write(*,*) "Fractal Optics Test"

  call test_fractaloptics()

  write(*,*) "Done"
end program


subroutine check(status)
  use netcdf

  implicit none
  integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
end subroutine check  

subroutine test_fractaloptics()
  use carma_precision_mod
  use carma_constants_mod
  use carma_enums_mod
  use carma_types_mod
  use carmaelement_mod
  use carmagroup_mod
  use carmastate_mod
  use carma_mod
  use test2nc_mod
  use netcdf

  implicit none

  ! Model dimensions
  integer, parameter    :: NELEM        = 1
  integer, parameter    :: NBIN         = 5
  integer, parameter    :: NGROUP       = 1
  integer, parameter    :: NSOLUTE      = 0
  integer, parameter    :: NGAS         = 0 
  integer, parameter    :: NREFIDX      = 1  ! Number of refractive indices per element
  integer, parameter    :: NRUN         = 2  ! NRUN max value is 2, one for fractal and one for standard
  integer, parameter    :: NWAVE        = 30 
  
  ! Model grid
  integer               :: bins(NBIN)

  ! CARMA object
  type(carma_type), target            :: carma
  type(carma_type), pointer           :: carma_ptr
  integer                             :: rc = 0

  ! Netcdf variables
  type(nc_type)      :: ncflgs
  character(len=80)  :: filename_in = "/glade/u/home/iquaglia/Solid_SAI/rindex/alumina.nc" ! input file for refractive index
  ! refractive index can be define later  without loading any file
  character(len=80)  :: filename_out = "carma_fractalopticstest.nc"
  integer            :: inid, outid
  integer            :: wave_varid, real_varid, imag_varid
  integer            :: date = 00010101

  ! Indexes for iterations
  integer               :: i
  integer               :: ielem
  integer               :: iwave
  integer               :: ibin
  integer               :: igroup
  integer               :: irun

  ! Input for optics 
  real(kind=f)          :: wave(NWAVE)
  real(kind=f)          :: r_refidx(NWAVE)
  real(kind=f)          :: i_refidx(NWAVE)
  
  ! Input for Group/element
  integer, parameter    :: I_ALUMINUM = 1
  
  character(len=80)     :: name_run(NRUN) = (/'fractal', 'standard'/)
  real(kind=f)          :: rmin   = 15.0e-6_f   ! minimum radius [cm]
  real(kind=f)          :: rmrat  = 2.0_f       ! the volume ratio between bins
  real(kind=f)          :: rho    = 3.51_f      ! dry density of particles [g cm-3]
  
  real(kind=f)          :: df(NBIN,NRUN)        ! aerosol fractal dimension
  real(kind=f)          :: rmon   = 15.0e-6_f   ! monomer radius [cm] == rmin
  real(kind=f)          :: falpha = 1._f        ! satellite aerosol fractal packing coefficient

  ! Output for group variables - varible definitions in test2nc_mod.F90 
  real(kind=f)          :: r(NBIN, NRUN), rmass(NBIN, NRUN), vol(NBIN, NRUN), dr(NBIN, NRUN), dm(NBIN, NRUN)
  real(kind=f)          :: rrat(NBIN, NRUN), rprat(NBIN, NRUN), arat(NBIN, NRUN), nmon(NBIN, NRUN) 
  real(kind=f)          :: rg(NBIN, NRUN), rp(NBIN, NRUN)
  complex(kind=f)       :: refidx(NWAVE, NREFIDX)

  real(kind=f)          :: qext(NWAVE, NBIN, NRUN)
  real(kind=f)          :: ssa(NWAVE, NBIN, NRUN)
  real(kind=f)          :: asy(NWAVE, NBIN, NRUN)

  ! Working variables
  real(kind=f)          :: rf(NBIN), vpor(NBIN), upor(NBIN), gamma(NBIN), happel(NBIN)
  real(kind=f)          :: perm(NBIN), brinkman(NBIN), epsil(NBIN), omega(NBIN)

  

  df(:,1) = 1.6_f ! fractal
  df(:,2) = 3._f  ! standard

  ! Define Wavelenght and refractive index 
  ! Open netCDF file  
  CALL check(nf90_open(filename_in, nf90_nowrite, inid))
  
  ! Get the varids of the wavelength, real and imag refractive index
  call check( nf90_inq_varid(inid, "wavelength", wave_varid) )
  call check( nf90_inq_varid(inid, "real", real_varid) )
  call check( nf90_inq_varid(inid, "imag", imag_varid) )

  ! Get the values
  call check(nf90_get_var(inid,wave_varid,wave)) 
  call check(nf90_get_var(inid,real_varid,r_refidx))
  call check(nf90_get_var(inid,imag_varid,i_refidx))
  
  ! Close the file
  call check( nf90_close(inid) )

  ! Define wavelengths and refractive index here insted of opening a netCDF file 
  ! wave     = (/.., ../)  
  ! r_refidx = (/.., ../)
  ! i_refidx = (/.., ../)
  
  ! Construct a complex refractive index.
  refidx(:, 1) = cmplx(r_refidx(:), i_refidx(:), kind=f)


  ! Create bin 
  do ibin = 1,NBIN
        bins(ibin) = ibin
  end do


  ! NetCDF definition
  call ncdef_optics( filename_out, NWAVE, NBIN, NRUN, name_run, outid)

  
  ! Define the particle-grid extent of the CARMA test 
  do irun = 1, NRUN 
        call CARMA_Create(carma, NBIN, NELEM, NGROUP, NSOLUTE, NGAS, NWAVE, rc, wave=wave, NREFIDX=NREFIDX)
        if (rc < 0) stop "    *** CARMA_Create FAILED ***"
        carma_ptr => carma

        ! Test1: do fractal optics calulcation
        if (name_run(irun) .EQ. 'fractal') then
               print*,name_run(irun) 
                call CARMAGROUP_Create(carma, 1, "fractal", rmin, rmrat, &
                        I_SPHERE, 1._f, .FALSE., rc, is_fractal=.TRUE., &
                        do_mie=.true., imiertn=I_MIERTN_BOTET1997, &
                        rmon=rmon, df=df, falpha=falpha)

                if (rc < 0) stop "    *** irun 1 CARMAGROUP_Create FAILED ***"

                ! Define the element
                call CARMAELEMENT_Create(carma, 1, 1, "Aerosol", rho, I_INVOLATILE, I_ALUMINUM, rc, refidx=refidx(:,1) )
                if (rc < 0) stop "    *** irun 1 CARMAELEMENT_Create FAILED ***"
           

        ! Test2: do standard mie calculation for comparison
        else if (name_run(irun) .EQ. 'standard') then  
                call CARMAGROUP_Create(carma, 1, "sphere", rmin, rmrat, &
                         I_SPHERE, 1._f, .FALSE., rc, is_fractal=.FALSE., &
                         do_mie=.true., imiertn=I_MIERTN_TOON1981)
                if (rc < 0) stop "    ***  irun 2 CARMAGROUP_CreateFAILED ***"

                ! Define the element
                call CARMAELEMENT_Create(carma, 1, 1, "Aerosol", rho, I_INVOLATILE, I_ALUMINUM, rc, refidx=refidx(:,1) )
                if (rc < 0) stop "    *** irun 2 CARMAELEMENT_Create FAILED ***"
        end if

      
        ! Setup the CARMA processes to exercise
        call CARMA_Initialize(carma, rc, do_pheat=.TRUE.)
       
        if (rc < 0) stop "    *** CARMA_Initialize FAILED ***"
        do igroup = 1, NGROUP           
                call CARMAGROUP_Get(carma, igroup, rc, r=r(:,irun), dr=dr(:,irun), rmass=rmass(:,irun), dm=dm(:,irun), vol=vol(:,irun), & 
                        qext=qext(:,:,irun), ssa=ssa(:,:,irun), asym=asy(:,:,irun), arat=arat(:,irun), rrat=rrat(:,irun), &
                        rprat=rprat(:,irun), nmon=nmon(:,irun)) 
        end do
        
        ! Calculation of agglomerate radius and mobility radius
        if (name_run(irun) .EQ. 'fractal') then

        ! agglomerate radius
        rg(:,irun) = falpha * rmon * nmon(:,irun)**(1.0_f/df(:,irun))

        ! mobility radius for permeable aggregates from Vainshtein et al 2003
                rf(:) = (1.0_f/falpha)**(1.0_f/df(:,irun)) * r(:,irun)**(3.0_f/df(:,irun)) * rmon**(1.0_f-3.0_f/df(:,irun))
                vpor(:) = 1.0_f - (nmon(:,irun))**(1.0_f - 3.0_f/df(:,irun))           ! Volume average porosity (eq. 3.2)
                upor(:) = 1.0_f - (1.0_f - vpor(:)) * sqrt(df(:,irun)/3.0_f)                ! Uniform poroisty (eq. 3.10)
                gamma(:) = (1.0_f - upor(:))**(1.0_f/3.0_f)
                happel(:) = 2.0_f/(9.0_f*(1.0_f-upor(:)))*   &                            ! Happel permeability model
                         (3.0_f-4.5_f*gamma(:)+4.5_f*gamma(:)**5.0_f-3.0_f*gamma(:)**6.0_f)/  &
                         (3.0_f+2.0_f*gamma(:)**5.0_f)
                perm(:) = happel(:) * rmon**2.0_f                                         ! Permeability (eq. 3.3) 
                brinkman(:) = nmon(:,irun)**(1.0_f/df(:,irun)) * 1.0_f/sqrt(happel(:))    ! Brinkman parameter (eq. 3.9)
                epsil(:) = 1.0_f - brinkman(:)**(-1.0_f)*tanh(brinkman(:))                    !
                omega(:) = 2.0_f/3.0_f*epsil(:)/(2.0_f/3.0_f+epsil(:)/brinkman(:)**2.0_f)     ! drag coefficient (eq. 2.7)
                rp(:,irun) = rf(:) * omega(:)                                                                   ! mobility radius for permeable aggregates
        else
            rg(:,irun) = falpha * rmon * nmon(:,irun)**(1.0_f/df(:,irun))
            rp(:,irun) = r(:,irun)
        end if
        
        call CARMA_Destroy(carma, rc)
        if (rc /=0) stop "    *** CARMA_Destroy FAILED ***"
   
  end do

  ! Write NETCDF file
  call ncwrt_optics(outid, NWAVE, NBIN, NRUN, wave, bins, name_run,  &
                r_refidx=r_refidx, i_refidx=i_refidx, r=r, rmass=rmass, vol=vol, dr=dr, dm=dm, & 
                rrat=rrat, rprat=rprat, arat=arat, df=df,rp=rp, rg=rg, qext=qext, asy=asy, ssa=ssa, nmon=nmon)

        
  ! Close the output NETCDF file
  call ncclose(outid)
  
  if (rc /=0) write(*, *) "    *** FAILED ***, rc=", rc
          
end subroutine
