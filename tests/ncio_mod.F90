! NETCDF I/O
! C. Bardeen
! June 2018
!
! This file contains the following subroutines, related to NETCDF input and
! output of model results:
!
!   ncabort
!   ncclose
!   nccreate
!   ncdefint
!   ncdefreal
!   ncrdint
!   ncrdint2d
!   ncrdint3d
!   ncrdndim
!   ncrdreal
!   ncrdreal2d
!   ncrdreal3d
!   ncopen
!   ncwrtdims
!   wcwrtint
!   ncwrttime
!   wcwrtreal
!
! These routines are wrappers on the F77 version of the NETCF libraries from
! UNIDATA, which need to built and installed outside of this project. The paths
! to the include files and libraries need to be added to the Makefile. This
! source file is freeform, and either needs to be modified to follow F77 standard
! or the compiler switches for freeform need to be added to the makefile.
!
!
! Pseudocode examples for reading & writing NETCDF files ...
!
!
! Reading input ...
!
!     ! Open an existing file
!     ncopen(filename, NF_NOWRITE)
!
!     ! Get the size of the dimensions
!     foreach dimension
!       ndim = ncrdndim(fileid, name)
!       dim(ndim) = ncrdint(fileid, name, ndim)
!       dim(ndim) = ncrdreal(fileid, name, ndim)
!     endfor
!
!     ! Input once per time loop
!     foreach time
!
!       ! Read in the data
!       foreach variable to be read
!
!         ! Use one of these calls (ncrdint, ncrdint2d, ncrdint3d
!         2dvar(nlat, nlon) = ncrdint2d(fileid, name, itime, nlat, nlon)
!         2dvar(nlat, nlon) = ncrdreal2d(fileid, name, itime, nlat, nlon)
!         3dvar(nlev, lat, nlon) = ncrdint3d(fileid, name, itime, nlev, nlat, nlon)
!         3dvar(nlev, lat, nlon) = ncrdreal3d(fileid, name, itime, nlev, nlat, nlon)
!       endfor
!     endfor
!
!     ! Close the file
!     ncclose(fileid)
!
!
! Writing output ...
!
!     ! Create a new file and define the dimensions
!     nccreate(filename, time_name, time_lname, time_untis, nlevs, nlats, nlons)
!
!     ! Define each variable    
!     foreach variable to be output
!       ncdefint(fileid, nlevs, nlats, nlons, name, long_name, units)
!       ncdefreal(fileid, nlevs, nlats, nlons, name, long_name, units)
!     endfor
!     ncdefend(fileid)
!
!     ! Output the dimensions
!     ncwrtdims(fileid, nlevs, levs, nlats, lats, nlons, lons)
!
!     ! Output once per time loop
!     foreach time
!
!       ! Add a time to the time dimension
!       ncwrttime(fileid, itime, time, date)
!
!       ! Output the data
!       foreach variable to be output
!
!         ! Use one of these calls (ncwrtint, ncwrtreal)
!         ncwrtxxx(fileid, name, itime, nlev, nlat, nlon, ilev, ilat, ilon, value)
!       endfor
!     endfor
!
!     ! Close the file
!     ncclose(fileid)
!     
!     
! Note: In all of these routines, time is considered to be an unlimited
! dimension.

module ncio_mod

  ! types
  use carma_precision_mod

  implicit none
  
  private

  interface ncrdint
    module procedure ncrdint1d
    module procedure ncrdint2d
    module procedure ncrdint3d
  end interface

  public ncrdint

  interface ncrdreal
    module procedure ncrdreal1d
    module procedure ncrdreal2d
    module procedure ncrdreal3d
    module procedure ncrdreal3dw
  end interface

  public ncrdreal
  
  public ncabort
  public ncabortname
  public ncclose
  public ncdefend
  public ncdefint
  public ncdefreal
  public ncdefbinreal
  public ncdefreal_optics
  public ncdefreal_height
  public ncdefbin
  public ncsync
  public ncwrtdims
  public ncwrtdims_optics
  public ncwrtdims_height
  public ncwrttime
  public ncwrtint
  public ncwrtreal
  public ncwrtbinreal
  public ncwrtreal_optics
  public ncwrtbin
  
  public nccreate
  public nccreate_optics
  public nccreate_height
  public ncopen
  public ncrdndim

  contains


      
  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Reports an error that occurred while trying to manipulate a NETCDF file. =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  routine  - character[*], string indicating where the problem was found   =*                                                           =*
  !=  status   - integer, error indicator                                      =*                                                           =*
  !-----------------------------------------------------------------------------*
  subroutine ncabort(routine, status)
    include 'netcdf.inc'

    ! input:
    integer status
    character(LEN=*) routine

    ! local:
      
    write(*,'(A,A)') 'NETCDF error detected in ', routine
    if (status .ne. NF_NOERR) then
      write(*,'(A)') trim(NF_STRERROR(status))
    end if
    stop

  end subroutine ncabort


  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Reports an error that occurred while trying to manipulate a NETCDF file. =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  routine  - character[*], string indicating where the problem was found   =*
  !=  name     - character[*], string indicating field that had problems       =*
  !=  status   - integer, error indicator                                      =*
  !-----------------------------------------------------------------------------*

  subroutine ncabortname(routine, name, status)
   include 'netcdf.inc'

    ! input:
    integer status
    character(LEN=*) routine
    character(LEN=*) name

    ! local:
  
    write(*,'(A,A)') 'NETCDF error detected in ', routine
    write(*,'(A,A)') 'for name = ', name
    if (status .ne. NF_NOERR) then
      write(*,'(A)') trim(NF_STRERROR(status))
    end if
    stop

  end subroutine ncabortname


  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Creates a new NETCDF file and defines the dimensions that will be used   =*
  !=  with variables written to this file. An unlimited time dimension is      =*
  !=  always assumed to exist. The levels, lats, and lons are optional and can =*
  !=  be set to zero if not used; however, if lons are defined then lats must  =*
  !=  exist too.                                                               =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  filename   - character[*], name of the file                              =*
  !=  time_name  - character[*], name for the time                             =*
  !=  time_lname - character[*], long_name for the time                        =*
  !=  time_units - character[*], units for the time                            =*
  !=  nlevs      - integer, number of vertical levels (midpoints)              =*
  !=  lev_name   - character[*], name for the levels                           =*
  !=  lev_lname  - character[*], long_name for the levels                      =*
  !=  lev_units  - character[*], units for the levels                          =*
  !=  nlats      - integer, number of latitudes                                =*
  !=  nlons      - integer, number of longitudes                               =*
  !-----------------------------------------------------------------------------*
  !=  RETURNS:                                                                 =*
  !=  fileid  - integer, handle for open netcdf file                           =*
  !-----------------------------------------------------------------------------*

  function nccreate(filename, time_name, time_lname, time_units, nlevs, nlats, nlons, nbins)
   include 'netcdf.inc'

    ! input:
    character(LEN=*) filename
    character(LEN=*) time_name
    character(LEN=*) time_lname
    character(LEN=*) time_units
    integer nlevs
    integer nlats
    integer nlons
    integer nbins

    ! function value:
    integer nccreate

    ! local:
    integer status
    integer fileid
    integer timedid, levdid, ilevdid, latdid, londid, bindid
    integer timeid, levid, ilevid, latid, lonid, binid
    integer dims(1)
 
    ! Create a new NETCDF file.
    status = NF_CREATE(filename, NF_CLOBBER, fileid)
    if (status .ne. NF_NOERR) call ncabortname('nccreate', filename, status)
 
    ! Always have an unlimited time dimension
    status = NF_DEF_DIM(fileid, time_name, NF_UNLIMITED, timedid)
    if (status .ne. NF_NOERR) call ncabortname('nccreate', 'time', status)
  
    dims(1) = timedid
#ifdef SINGLE
    status = NF_DEF_VAR(fileid, time_name, NF_FLOAT, 1, dims(1), timeid)
#else
    status = NF_DEF_VAR(fileid, time_name, NF_DOUBLE, 1, dims(1), timeid)
#endif
    if (status .ne. NF_NOERR) call ncabortname('nccreate', 'time', status)
    status = NF_PUT_ATT_TEXT(fileid, timeid, 'long_name', LEN(time_lname), time_lname)
    if (status .ne. NF_NOERR) call ncabortname('nccreate', 'time', status)
    status = NF_PUT_ATT_TEXT(fileid, timeid, 'units', LEN(time_units), time_units)
    if (status .ne. NF_NOERR) call ncabortname('nccreate', 'time', status)

    ! And a date variable to go with it
    dims(1) = timedid
    status = NF_DEF_VAR(fileid, "date", NF_INT, 1, dims(1), timeid)
    if (status .ne. NF_NOERR) call ncabortname('nccreate', 'date', status)
    status = NF_PUT_ATT_TEXT(fileid, timeid, 'long_name', LEN('Current Date'), 'Current Date')
    if (status .ne. NF_NOERR) call ncabortname('nccreate', 'date', status)
    status = NF_PUT_ATT_TEXT(fileid, timeid, 'units', LEN('yyyymmdd'), 'yyyymmdd')
    if (status .ne. NF_NOERR) call ncabortname('nccreate', 'date', status)

    ! Define levs, ilevs, lats, and lons as needed.
    if (nlevs .gt. 1) then
      status = NF_DEF_DIM(fileid, 'lev', nlevs, levdid)
      if (status .ne. NF_NOERR) call ncabortname('nccreate', 'lev', status)

      dims(1) = levdid
#ifdef SINGLE
      status = NF_DEF_VAR(fileid, 'lev', NF_FLOAT, 1, dims(1), levid)
#else
      status = NF_DEF_VAR(fileid, 'lev', NF_DOUBLE, 1, dims(1), levid)
#endif
      if (status .ne. NF_NOERR) call ncabortname('nccreate', 'lev', status)
      status = NF_PUT_ATT_TEXT(fileid, levid, 'long_name', LEN('hybrid level at midpoints'), 'hybrid level at midpoints')
      if (status .ne. NF_NOERR) call ncabortname('nccreate', 'lev', status)
      status = NF_PUT_ATT_TEXT(fileid, levid, 'units', LEN('hPa'), 'hPa')
      if (status .ne. NF_NOERR) call ncabortname('nccreate', 'lev', status)


      status = NF_DEF_DIM(fileid, 'ilev', nlevs, levdid)
      if (status .ne. NF_NOERR) call ncabortname('nccreate', 'ilev', status)

      dims(1) = levdid
#ifdef SINGLE
      status = NF_DEF_VAR(fileid, 'ilev', NF_FLOAT, 1, dims(1), levid)
#else
      status = NF_DEF_VAR(fileid, 'ilev', NF_DOUBLE, 1, dims(1), levid)
#endif
      if (status .ne. NF_NOERR) call ncabortname('nccreate', 'ilev', status)
      status = NF_PUT_ATT_TEXT(fileid, levid, 'long_name', LEN('hybrid layer at interfaces'), 'hybrid layer at interfaces')
      if (status .ne. NF_NOERR) call ncabortname('nccreate', 'ilev', status)
      status = NF_PUT_ATT_TEXT(fileid, levid, 'units', LEN('hPa'), 'hPa')
      if (status .ne. NF_NOERR) call ncabortname('nccreate', 'ilev', status)

    end if

    if (nlats .gt. 1) then
      status = NF_DEF_DIM(fileid, 'lat', nlats, latdid)
      if (status .ne. NF_NOERR) call ncabortname('nccreate', 'lat', status)

      dims(1) = latdid
#ifdef SINGLE
      status = NF_DEF_VAR(fileid, 'lat', NF_FLOAT, 1, dims(1), latid)
#else
      status = NF_DEF_VAR(fileid, 'lat', NF_DOUBLE, 1, dims(1), latid)
#endif
      if (status .ne. NF_NOERR) call ncabortname('nccreate', 'lat', status)
      status = NF_PUT_ATT_TEXT(fileid, latid, 'long_name', LEN('Latitude'), 'Latitude')
      if (status .ne. NF_NOERR) call ncabortname('nccreate', 'lat', status)
      status = NF_PUT_ATT_TEXT(fileid, latid, 'units', LEN('degN'), 'degN')
      if (status .ne. NF_NOERR) call ncabortname('nccreate', 'lat', status)
    end if

    if (nlons .gt. 1) then
      status = NF_DEF_DIM(fileid, 'lon', nlons, londid)
      if (status .ne. NF_NOERR) call ncabortname('nccreate', 'lon', status)

      dims(1) = londid
#ifdef SINGLE
      status = NF_DEF_VAR(fileid, 'lon', NF_FLOAT, 1, dims(1), lonid)
#else
      status = NF_DEF_VAR(fileid, 'lon', NF_DOUBLE, 1, dims(1), lonid)
#endif
      if (status .ne. NF_NOERR) call ncabortname('nccreate', 'lon', status)
      status = NF_PUT_ATT_TEXT(fileid, lonid, 'long_name', LEN('Longitude'),  'Longitude')
      if (status .ne. NF_NOERR) call ncabortname('nccreate', 'lon', status)
      status = NF_PUT_ATT_TEXT(fileid, lonid, 'units', LEN('degE'), 'degE')
      if (status .ne. NF_NOERR) call ncabortname('nccreate', 'lon', status)
    end if

    if (nbins .gt. 1) then
      status = NF_DEF_DIM(fileid, 'bin', nbins, bindid)
      if (status .ne. NF_NOERR) call ncabortname('nccreate', 'bin', status)

      dims(1) = bindid
#ifdef SINGLE
      status = NF_DEF_VAR(fileid, 'bin', NF_FLOAT, 1, dims(1), binid)
#else
      status = NF_DEF_VAR(fileid, 'bin', NF_DOUBLE, 1, dims(1), binid)
#endif
      if (status .ne. NF_NOERR) call ncabortname('nccreate', 'bin', status)
      status = NF_PUT_ATT_TEXT(fileid, binid, 'long_name', LEN('bin'),  'bin')
      if (status .ne. NF_NOERR) call ncabortname('nccreate', 'bin', status)
      status = NF_PUT_ATT_TEXT(fileid, binid, 'units', LEN('bin'), 'bin')
      if (status .ne. NF_NOERR) call ncabortname('nccreate', 'bin', status)
    end if

    nccreate = fileid

  end function nccreate

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Creates a new NETCDF file and defines the dimensions that will be used   =*
  !=  with variables written to this file. An unlimited time dimension is      =*
  !=  always assumed to exist. The levels, lats, and lons are optional and can =*
  !=  be set to zero if not used; however, if lons are defined then lats must  =*
  !=  exist too.                                                               =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  filename   - character[*], name of the file                              =*
  !=  time_name  - character[*], name for the time                             =*
  !=  time_lname - character[*], long_name for the time                        =*
  !=  time_units - character[*], units for the time                            =*
  !=  nlevs      - integer, number of vertical levels (midpoints)              =*
  !=  lev_name   - character[*], name for the levels                           =*
  !=  lev_lname  - character[*], long_name for the levels                      =*
  !=  lev_units  - character[*], units for the levels                          =*
  !=  nlats      - integer, number of latitudes                                =*
  !=  nlons      - integer, number of longitudes                               =*
  !-----------------------------------------------------------------------------*
  !=  RETURNS:                                                                 =*
  !=  fileid  - integer, handle for open netcdf file                           =*
  !-----------------------------------------------------------------------------*

  function nccreate_height(filename, time_name, time_lname, time_units, nlevs, nlats, nlons, nbins)
  include 'netcdf.inc'

  ! input:
  character(LEN=*) filename
  character(LEN=*) time_name
  character(LEN=*) time_lname
  character(LEN=*) time_units
  integer nlevs
  integer nlats
  integer nlons
  integer nbins

  ! function value:
  integer nccreate_height

  ! local:
  integer status
  integer fileid
  integer timedid, levdid, ilevdid, latdid, londid, bindid
  integer timeid, levid, ilevid, latid, lonid, binid
  integer dims(1)

  ! Create a new NETCDF file.
  status = NF_CREATE(filename, NF_CLOBBER, fileid)
  if (status .ne. NF_NOERR) call ncabortname('nccreate_height', filename, status)

  ! Always have an unlimited time dimension
  status = NF_DEF_DIM(fileid, time_name, NF_UNLIMITED, timedid)
  if (status .ne. NF_NOERR) call ncabortname('nccreate_height', 'time', status)

  dims(1) = timedid
#ifdef SINGLE
  status = NF_DEF_VAR(fileid, time_name, NF_FLOAT, 1, dims(1), timeid)
#else
  status = NF_DEF_VAR(fileid, time_name, NF_DOUBLE, 1, dims(1), timeid)
#endif
  if (status .ne. NF_NOERR) call ncabortname('nccreate_height', 'time', status)
  status = NF_PUT_ATT_TEXT(fileid, timeid, 'long_name', LEN(time_lname), time_lname)
  if (status .ne. NF_NOERR) call ncabortname('nccreate_height', 'time', status)
  status = NF_PUT_ATT_TEXT(fileid, timeid, 'units', LEN(time_units), time_units)
  if (status .ne. NF_NOERR) call ncabortname('nccreate_height', 'time', status)

  ! And a date variable to go with it
  dims(1) = timedid
  status = NF_DEF_VAR(fileid, "date", NF_INT, 1, dims(1), timeid)
  if (status .ne. NF_NOERR) call ncabortname('nccreate_height', 'date', status)
  status = NF_PUT_ATT_TEXT(fileid, timeid, 'long_name', LEN('Current Date'), 'Current Date')
  if (status .ne. NF_NOERR) call ncabortname('nccreate_height', 'date', status)
  status = NF_PUT_ATT_TEXT(fileid, timeid, 'units', LEN('yyyymmdd'), 'yyyymmdd')
  if (status .ne. NF_NOERR) call ncabortname('nccreate_height', 'date', status)

  ! Define levs, ilevs, lats, and lons as needed.
  if (nlevs .gt. 1) then
    status = NF_DEF_DIM(fileid, 'Height', nlevs, levdid)
    if (status .ne. NF_NOERR) call ncabortname('nccreate_height', 'Height', status)

    dims(1) = levdid
#ifdef SINGLE
    status = NF_DEF_VAR(fileid, 'Height', NF_FLOAT, 1, dims(1), levid)
#else
    status = NF_DEF_VAR(fileid, 'Height', NF_DOUBLE, 1, dims(1), levid)
#endif
    if (status .ne. NF_NOERR) call ncabortname('nccreate_height', 'Height', status)
    status = NF_PUT_ATT_TEXT(fileid, levid, 'long_name', LEN('Geographic Height'), 'Geographic Height')
    if (status .ne. NF_NOERR) call ncabortname('nccreate_height', 'Height', status)
    status = NF_PUT_ATT_TEXT(fileid, levid, 'units', LEN('km'), 'km')
    if (status .ne. NF_NOERR) call ncabortname('nccreate_height', 'Height', status)

  end if

  if (nlats .gt. 1) then
    status = NF_DEF_DIM(fileid, 'lat', nlats, latdid)
    if (status .ne. NF_NOERR) call ncabortname('nccreate_height', 'lat', status)

    dims(1) = latdid
#ifdef SINGLE
    status = NF_DEF_VAR(fileid, 'lat', NF_FLOAT, 1, dims(1), latid)
#else
    status = NF_DEF_VAR(fileid, 'lat', NF_DOUBLE, 1, dims(1), latid)
#endif
    if (status .ne. NF_NOERR) call ncabortname('nccreate_height', 'lat', status)
    status = NF_PUT_ATT_TEXT(fileid, latid, 'long_name', LEN('Latitude'), 'Latitude')
    if (status .ne. NF_NOERR) call ncabortname('nccreate_height', 'lat', status)
    status = NF_PUT_ATT_TEXT(fileid, latid, 'units', LEN('degN'), 'degN')
    if (status .ne. NF_NOERR) call ncabortname('nccreate_height', 'lat', status)
  end if

  if (nlons .gt. 1) then
    status = NF_DEF_DIM(fileid, 'lon', nlons, londid)
    if (status .ne. NF_NOERR) call ncabortname('nccreate_height', 'lon', status)

    dims(1) = londid
#ifdef SINGLE
    status = NF_DEF_VAR(fileid, 'lon', NF_FLOAT, 1, dims(1), lonid)
#else
    status = NF_DEF_VAR(fileid, 'lon', NF_DOUBLE, 1, dims(1), lonid)
#endif
    if (status .ne. NF_NOERR) call ncabortname('nccreate_height', 'lon', status)
    status = NF_PUT_ATT_TEXT(fileid, lonid, 'long_name', LEN('Longitude'),  'Longitude')
    if (status .ne. NF_NOERR) call ncabortname('nccreate_height', 'lon', status)
    status = NF_PUT_ATT_TEXT(fileid, lonid, 'units', LEN('degE'), 'degE')
    if (status .ne. NF_NOERR) call ncabortname('nccreate_height', 'lon', status)
  end if

  if (nbins .gt. 1) then
    status = NF_DEF_DIM(fileid, 'bin', nbins, bindid)
    if (status .ne. NF_NOERR) call ncabortname('nccreate_height', 'bin', status)

    dims(1) = bindid
#ifdef SINGLE
    status = NF_DEF_VAR(fileid, 'bin', NF_FLOAT, 1, dims(1), binid)
#else
    status = NF_DEF_VAR(fileid, 'bin', NF_DOUBLE, 1, dims(1), binid)
#endif
    if (status .ne. NF_NOERR) call ncabortname('nccreate_height', 'bin', status)
    status = NF_PUT_ATT_TEXT(fileid, binid, 'long_name', LEN('bin'),  'bin')
    if (status .ne. NF_NOERR) call ncabortname('nccreate_height', 'bin', status)
    status = NF_PUT_ATT_TEXT(fileid, binid, 'units', LEN('bin'), 'bin')
    if (status .ne. NF_NOERR) call ncabortname('nccreate_height', 'bin', status)
  end if

  nccreate_height = fileid

end function nccreate_height

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Creates a new NETCDF file and defines the dimensions that will be used   =*
  !=  with variables written to this file. An unlimited time dimension is      =*
  !=  always assumed to exist. The levels, lats, and lons are optional and can =*
  !=  be set to zero if not used; however, if lons are defined then lats must  =*
  !=  exist too.                                                               =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  filename   - character[*], name of the file                              =*
  !=  time_name  - character[*], name for the time                             =*
  !=  time_lname - character[*], long_name for the time                        =*
  !=  time_units - character[*], units for the time                            =*
  !=  nruns      - integer, number of runs                                     =*
  !=  nwaves     - integer, number of latitudes                                =*
  !=  nbins      - integer, number of longitudes                               =*
  !-----------------------------------------------------------------------------*
  !=  RETURNS:                                                                 =*
  !=  fileid  - integer, handle for open netcdf file                           =*
  !-----------------------------------------------------------------------------*

  function nccreate_optics(filename, nwaves, nbins)
    include 'netcdf.inc'

    ! input:
    character(LEN=*) filename
    integer nwaves
    integer nbins

    ! function value:
    integer nccreate_optics

    ! local:
    integer status
    integer fileid
    integer wavedid, bindid
    integer waveid, binid
    integer dims(1)

    ! Create a new NETCDF file.
    status = NF_CREATE(filename, NF_CLOBBER, fileid)
    if (status .ne. NF_NOERR) call ncabortname('nccreate_optics', filename, status)

    ! Define levs, ilevs, lats, and lons as needed.
    if (nwaves .gt. 1) then
      status = NF_DEF_DIM(fileid, 'wave', nwaves, wavedid)
      if (status .ne. NF_NOERR) call ncabortname('nccreate_optics', 'wave', status)

      dims(1) = wavedid
#ifdef SINGLE
        status = NF_DEF_VAR(fileid, 'wave', NF_FLOAT, 1, dims(1), waveid)
#else
        status = NF_DEF_VAR(fileid, 'wave', NF_DOUBLE, 1, dims(1), waveid)
#endif
      if (status .ne. NF_NOERR) call ncabortname('nccreate_optics', 'wave', status)
      status = NF_PUT_ATT_TEXT(fileid, waveid, 'long_name', LEN('wavelength'), 'wavelength')
      if (status .ne. NF_NOERR) call ncabortname('nccreate_optics', 'wave', status)
      status = NF_PUT_ATT_TEXT(fileid, waveid, 'units', LEN('cm'), 'cm')
      if (status .ne. NF_NOERR) call ncabortname('nccreate_optics', 'wave', status)
    end if

    if (nbins .gt. 1) then
      status = NF_DEF_DIM(fileid, 'bin', nbins, bindid)
      if (status .ne. NF_NOERR) call ncabortname('nccreate_optics', 'bin', status)

      dims(1) = bindid
#ifdef SINGLE
        status = NF_DEF_VAR(fileid, 'bin', NF_FLOAT, 1, dims(1), binid)
#else
        status = NF_DEF_VAR(fileid, 'bin', NF_DOUBLE, 1, dims(1), binid)
#endif
      if (status .ne. NF_NOERR) call ncabortname('nccreate_optics', 'bin', status)
      status = NF_PUT_ATT_TEXT(fileid, binid, 'long_name', LEN('bin'),  'bin')
      if (status .ne. NF_NOERR) call ncabortname('nccreate_optics', 'bin', status)
      status = NF_PUT_ATT_TEXT(fileid, binid, 'units', LEN('bin'), 'bin')
      if (status .ne. NF_NOERR) call ncabortname('nccreate_optics', 'bin', status)
    end if

    nccreate_optics= fileid

  end function nccreate_optics

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Ends the definition phase of creating the output file.                   =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  fileid    - integer, handle for open netcdf file                         =*
  !-----------------------------------------------------------------------------*
  subroutine ncdefend(fileid)
    include 'netcdf.inc'

    ! input:
    integer fileid

    ! local:
    integer  status

    ! End defintion mode.
    status = NF_endDEF(fileid)
    if (status .ne. NF_NOERR) call ncabort('ncdefend', status)
  
  end subroutine ncdefend

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Writes the bin dimension for a previously created NETCDF file.           =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  fileid    - integer, handle for open netcdf file                         =*
  !=  igroup    - integer, group index                                         =*
  !=  ibin      - integer, bin index                                           =*
  !-----------------------------------------------------------------------------*

  subroutine ncdefbin(fileid, nbins, name, long_name, units)
    include 'netcdf.inc'

    ! input:
    integer fileid
    integer nbins
    character(LEN=*) name
    character(LEN=*) long_name
    character(LEN=*) units

    ! local:
    integer status,varid,binid
    integer dims(1)

      status = NF_INQ_DIMID(fileid, 'bin', binid)
      if (status .ne. NF_NOERR) call ncabortname('ncdefbin', name, status)
        dims(1) = binid
#ifdef SINGLE
        status = NF_DEF_VAR(fileid, name, NF_FLOAT, 1, dims(1), varid)
#else
        status = NF_DEF_VAR(fileid, name, NF_DOUBLE, 1, dims(1), varid)
#endif
        if (status .ne. NF_NOERR) call ncabortname('ncdefbin', name, status)

    ! add attributes
      status = NF_PUT_ATT_TEXT(fileid, varid, 'long_name', LEN(long_name), long_name)
      if (status .ne. NF_NOERR) call ncabortname('nccreate', name, status)
      status = NF_PUT_ATT_TEXT(fileid, varid, 'units', LEN(units), units)
      if (status .ne. NF_NOERR) call ncabortname('nccreate', name, status)

  end subroutine ncdefbin

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Defines an integer output variable using previously defined dimensions.  =*
  !=  An unlimited time dimension is always assumed to exist. The name,        =*
  !=  long_name, and units strings describe the variable.                      =*
  !=  NOTE: If defined, nlevs, nlats, and nlons should be the same as the      =*
  !=  values used in nccreate.                                                 =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  fileid    - integer, handle for open netcdf file                         =*
  !=  nlevs     - integer, number of vertical levels                           =*
  !=  nlats     - integer, number of latitudes                                 =*
  !=  nlons     - integer, number of longitudes                                =*
  !=  name      - character[*], name of field                                  =*
  !=  long_name - character[*], description of field                           =*
  !=  units     - character[*], units of data in the field                     =*
  !-----------------------------------------------------------------------------*
  subroutine ncdefint(fileid, nlevs, nlats, nlons, name, long_name, units)
    include 'netcdf.inc'

    ! input:
    integer fileid
    integer nlevs
    integer nlats
    integer nlons
    character(LEN=*) name
    character(LEN=*) long_name
    character(LEN=*) units
  
    ! output: (on converted grid)

    ! local:
    integer status, varid, timeid, levid, latid, lonid
    integer dims(4)
    integer length
 
    ! Always have an unlimited time dimension
    status = NF_INQ_DIMID(fileid, 'time', timeid)
    if (status .ne. NF_NOERR) call ncabortname('ncdefint', name, status)

    ! Find the appropriate dimensions and link then to the definition
    ! on this field.
    !
    ! NOTE: To have lons, you must also have lats.
    if (nlevs .gt. 1) then
      status = NF_INQ_DIMID(fileid, 'lev', levid)
      if (status .ne. NF_NOERR) call ncabortname('ncdefint', name, status)
   
      ! Check the number of levels to see if the dimension that should be used is
      ! lev (midpoints) or ilev (interfaces).
      status = NF_INQ_DIMLEN(fileid, levid, length)
      if (status .ne. NF_NOERR) call ncabortname('ncdefint', name, status)
   
      if (nlevs .ne. length) then
        if (nlevs .eq. (length + 1)) then
          status = NF_INQ_DIMID(fileid, 'ilev', levid)
          if (status .ne. NF_NOERR) call ncabortname('ncdefint', name, status)
        else
          write(*,*) 'ncdefint: ERROR - Number of levels does not match lev or ilev.'
          call ncabortname('ncdefint', name, NF_NOERR)
        end if
      end if

      if (nlats .gt. 1) then
        status = NF_INQ_DIMID(fileid, 'lat', latid)
        if (status .ne. NF_NOERR) call ncabortname('ncdefint', name, status)
   
        if (nlons .gt. 1) then
          status = NF_INQ_DIMID(fileid, 'lon', lonid)
          if (status .ne. NF_NOERR) call ncabortname('ncdefint', name, status)

          dims(1) = lonid
          dims(2) = latid
          dims(3) = levid
          dims(4) = timeid

          status = NF_DEF_VAR(fileid, name, NF_INT, 4, dims(1:4), varid)
          if (status .ne. NF_NOERR) call ncabortname('ncdefint', name, status)
        else
          dims(1) = latid
          dims(2) = levid
          dims(3) = timeid

          status = NF_DEF_VAR(fileid, name, NF_INT, 3, dims(1:3), varid)
          if (status .ne. NF_NOERR) call ncabortname('ncdefint', name, status)
        end if
      else
        dims(1) = levid
        dims(2) = timeid

        status = NF_DEF_VAR(fileid, name, NF_INT, 2, dims(1:2), varid)
        if (status .ne. NF_NOERR) call ncabortname('ncdefint', name, status)
      end if
    else
      if (nlats .gt. 1) then
        status = NF_INQ_DIMID(fileid, 'lat', latid)
        if (status .ne. NF_NOERR) call ncabortname('ncdefint', name, status)
   
        if (nlons .gt. 1) then
          status = NF_INQ_DIMID(fileid, 'lon', lonid)
          if (status .ne. NF_NOERR) call ncabortname('ncdefint', name, status)

          dims(1) = lonid
          dims(2) = latid
          dims(3) = timeid

          status = NF_DEF_VAR(fileid, name, NF_INT, 3, dims(1:3), varid)
          if (status .ne. NF_NOERR) call ncabortname('ncdefint', name, status)
        else
          dims(1) = latid
          dims(2) = timeid

          status = NF_DEF_VAR(fileid, name, NF_INT, 2, dims(1:2), varid)
          if (status .ne. NF_NOERR) call ncabortname('ncdefint', name, status)
        end if
      else
        dims(1) = timeid

        status = NF_DEF_VAR(fileid, name, NF_INT, 1, dims(1), varid)
        if (status .ne. NF_NOERR) call ncabortname('ncdefint', name, status)
      end if
    end if

  end subroutine ncdefint


  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Defines a real output variable using previously defined dimensions. An   =*
  !=  unlimited time dimension is always assumed to exist. The name,           =*
  !=  long_name, and units strings describe the variable.                      =*
  !=  NOTE: If defined, nlevs, nlats, and nlons should be the same as the      =*
  !=  values used in nccreate.                                                 =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  fileid    - integer, handle for open netcdf file                         =*
  !=  nlevs     - integer, number of vertical levels (midpoint or interface)   =*
  !=  nlats     - integer, number of latitudes                                 =*
  !=  nlons     - integer, number of longitudes                                =*
  !=  name      - character[*], name of field                                  =*
  !=  long_name - character[*], description of field                           =*
  !=  units     - character[*], units of data in the field                     =*
  !-----------------------------------------------------------------------------*
  subroutine ncdefreal(fileid, nlevs, nlats, nlons, name, long_name, units, fill_value)
    include 'netcdf.inc'

    ! input:
    integer fileid
    integer nlevs
    integer nlats
    integer nlons
    character(LEN=*) name
    character(LEN=*) long_name
    character(LEN=*) units
    real(kind=f)     fill_value
    ! output: (on converted grid)

    ! local:
    integer status, varid, timeid, levid, latid, lonid
    integer dims(4)
    integer length
 
    ! Always have an unlimited time dimension
    status = NF_INQ_DIMID(fileid, 'time', timeid)
    if (status .ne. NF_NOERR) call ncabortname('ncdefreal', name, status)

    ! Find the appropriate dimensions and link then to the definition
    ! on this field.
    !
    ! NOTE: To have lons, you must also have lats.
    if (nlevs .gt. 1) then
      status = NF_INQ_DIMID(fileid, 'lev', levid)
      if (status .ne. NF_NOERR) call ncabortname('ncdefreal', name, status)
    
      ! Check the number of levels to see if the dimension that should be used is
      ! lev (midpoints) or ilev (interfaces).
      status = NF_INQ_DIMLEN(fileid, levid, length)
      if (status .ne. NF_NOERR) call ncabortname('ncdefreal', name, status)
   
      if (nlevs .ne. length) then
        if (nlevs .eq. (length + 1)) then
          status = NF_INQ_DIMID(fileid, 'ilev', levid)
          if (status .ne. NF_NOERR) call ncabortname('ncdefreal', name, status)
        else
          write(*,*) 'ncdefreal: ERROR - Number of levels does not match lev or ilev.'
          call ncabortname('ncdefreal', name, NF_NOERR)
        end if
      end if
    
      if (nlats .gt. 1) then
        status = NF_INQ_DIMID(fileid, 'lat', latid)
        if (status .ne. NF_NOERR) call ncabortname('ncdefreal', 'name', status)
   
        if (nlons .gt. 1) then
          status = NF_INQ_DIMID(fileid, 'lon', lonid)
          if (status .ne. NF_NOERR) call ncabortname('ncdefreal', name, status)

          dims(1) = lonid
          dims(2) = latid
          dims(3) = levid
          dims(4) = timeid

#ifdef SINGLE
          status = NF_DEF_VAR(fileid, name, NF_FLOAT, 4, dims(1:4), varid)
#else
          status = NF_DEF_VAR(fileid, name, NF_DOUBLE, 4, dims(1:4), varid)
#endif
          if (status .ne. NF_NOERR) call ncabortname('ncdefreal', name, status)
        else
          dims(1) = latid
          dims(2) = levid
          dims(3) = timeid

#ifdef SINGLE
          status = NF_DEF_VAR(fileid, name, NF_FLOAT, 3, dims(1:3), varid)
#else
          status = NF_DEF_VAR(fileid, name, NF_DOUBLE, 3, dims(1:3), varid)
#endif
          if (status .ne. NF_NOERR) call ncabortname('ncdefreal', name, status)
        end if
      else
        dims(1) = levid
        dims(2) = timeid

#ifdef SINGLE
        status = NF_DEF_VAR(fileid, name, NF_FLOAT, 2, dims(1:2), varid)
#else
        status = NF_DEF_VAR(fileid, name, NF_DOUBLE, 2, dims(1:2), varid)
#endif
        if (status .ne. NF_NOERR) call ncabortname('ncdefreal', name, status)
      end if
    else
      if (nlats .gt. 1) then
        status = NF_INQ_DIMID(fileid, 'lat', latid)
        if (status .ne. NF_NOERR) call ncabortname('ncdefreal', name, status)
   
        if (nlons .gt. 1) then
          status = NF_INQ_DIMID(fileid, 'lon', lonid)
          if (status .ne. NF_NOERR) call ncabortname('ncdefreal', name, status)

          dims(1) = lonid
          dims(2) = latid
          dims(3) = timeid

#ifdef SINGLE
          status = NF_DEF_VAR(fileid, name, NF_FLOAT, 3, dims(1:3), varid)
#else
          status = NF_DEF_VAR(fileid, name, NF_DOUBLE, 3, dims(1:3), varid)
#endif
          if (status .ne. NF_NOERR) call ncabortname('ncdefreal', name, status)
        else
          dims(1) = latid
          dims(2) = timeid

#ifdef SINGLE
          status = NF_DEF_VAR(fileid, name, NF_FLOAT, 2, dims(1:2), varid)
#else
          status = NF_DEF_VAR(fileid, name, NF_DOUBLE, 2, dims(1:2), varid)
#endif
          if (status .ne. NF_NOERR) call ncabortname('ncdefreal', name, status)
        end if
      else
        dims(1) = timeid

#ifdef SINGLE
        status = NF_DEF_VAR(fileid, name, NF_FLOAT, 1, dims(1), varid)
#else
        status = NF_DEF_VAR(fileid, name, NF_DOUBLE, 1, dims(1), varid)
#endif
        if (status .ne. NF_NOERR) call ncabortname('ncdefreal', name, status)
      end if
    end if  
 
    status = NF_DEF_VAR_FILL(fileid, varid, 0, fill_value)
    if (status .ne. NF_NOERR) call ncabortname('ncdefreal', name, status)

    ! Add the attributes      
    status = NF_PUT_ATT_TEXT(fileid, varid, 'long_name', LEN(long_name), long_name)
    if (status .ne. NF_NOERR) call ncabortname('nccreate', name, status)
    status = NF_PUT_ATT_TEXT(fileid, varid, 'units', LEN(units), units)
    if (status .ne. NF_NOERR) call ncabortname('nccreate', name, status)
    
  end subroutine ncdefreal


  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Defines a real output variable using previously defined dimensions. An   =*
  !=  unlimited time dimension is always assumed to exist. The name,           =*
  !=  long_name, and units strings describe the variable.                      =*
  !=  NOTE: If defined, nlevs, nlats, and nlons should be the same as the      =*
  !=  values used in nccreate.                                                 =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  fileid    - integer, handle for open netcdf file                         =*
  !=  nlevs     - integer, number of vertical levels (midpoint or interface)   =*
  !=  nlats     - integer, number of latitudes                                 =*
  !=  nlons     - integer, number of longitudes                                =*
  !=  name      - character[*], name of field                                  =*
  !=  long_name - character[*], description of field                           =*
  !=  units     - character[*], units of data in the field                     =*
  !-----------------------------------------------------------------------------*
  subroutine ncdefreal_height(fileid, nlevs, nlats, nlons, name, long_name, units)
  include 'netcdf.inc'

  ! input:
  integer fileid
  integer nlevs
  integer nlats
  integer nlons
  character(LEN=*) name
  character(LEN=*) long_name
  character(LEN=*) units

  ! output: (on converted grid)

  ! local:
  integer status, varid, timeid, levid, latid, lonid
  integer dims(4)
  integer length

  ! Always have an unlimited time dimension
  status = NF_INQ_DIMID(fileid, 'time', timeid)
  if (status .ne. NF_NOERR) call ncabortname('ncdefreal_height', name, status)

  ! Find the appropriate dimensions and link then to the definition
  ! on this field.
  !
  ! NOTE: To have lons, you must also have lats.
  if (nlevs .gt. 1) then
    status = NF_INQ_DIMID(fileid, 'Height', levid)
    if (status .ne. NF_NOERR) call ncabortname('ncdefreal_height', name, status)

    ! Check the number of levels to see if the dimension that should be used is
    ! lev (midpoints) or ilev (interfaces).
    status = NF_INQ_DIMLEN(fileid, levid, length)
    if (status .ne. NF_NOERR) call ncabortname('ncdefreal_height', name, status)

    if (nlats .gt. 1) then
      status = NF_INQ_DIMID(fileid, 'lat', latid)
      if (status .ne. NF_NOERR) call ncabortname('ncdefreal_height', name, status)

      if (nlons .gt. 1) then
        status = NF_INQ_DIMID(fileid, 'lon', lonid)
        if (status .ne. NF_NOERR) call ncabortname('ncdefreal_height', name, status)

        dims(1) = lonid
        dims(2) = latid
        dims(3) = levid
        dims(4) = timeid

#ifdef SINGLE
        status = NF_DEF_VAR(fileid, name, NF_FLOAT, 4, dims(1:4), varid)
#else
        status = NF_DEF_VAR(fileid, name, NF_DOUBLE, 4, dims(1:4), varid)
#endif
        if (status .ne. NF_NOERR) call ncabortname('ncdefreal_height', name, status)
      else
        dims(1) = latid
        dims(2) = levid
        dims(3) = timeid

#ifdef SINGLE
        status = NF_DEF_VAR(fileid, name, NF_FLOAT, 3, dims(1:3), varid)
#else
        status = NF_DEF_VAR(fileid, name, NF_DOUBLE, 3, dims(1:3), varid)
#endif
        if (status .ne. NF_NOERR) call ncabortname('ncdefreal_height', name, status)
      end if
    else
      dims(1) = levid
      dims(2) = timeid

#ifdef SINGLE
      status = NF_DEF_VAR(fileid, name, NF_FLOAT, 2, dims(1:2), varid)
#else
      status = NF_DEF_VAR(fileid, name, NF_DOUBLE, 2, dims(1:2), varid)
#endif
      if (status .ne. NF_NOERR) call ncabortname('ncdefreal_height', name, status)
    end if
  else
    if (nlats .gt. 1) then
      status = NF_INQ_DIMID(fileid, 'lat', latid)
      if (status .ne. NF_NOERR) call ncabortname('ncdefreal_height', name, status)

      if (nlons .gt. 1) then
        status = NF_INQ_DIMID(fileid, 'lon', lonid)
        if (status .ne. NF_NOERR) call ncabortname('ncdefreal_height', name, status)

        dims(1) = lonid
        dims(2) = latid
        dims(3) = timeid

#ifdef SINGLE
        status = NF_DEF_VAR(fileid, name, NF_FLOAT, 3, dims(1:3), varid)
#else
        status = NF_DEF_VAR(fileid, name, NF_DOUBLE, 3, dims(1:3), varid)
#endif
        if (status .ne. NF_NOERR) call ncabortname('ncdefreal_height', name, status)
      else
        dims(1) = latid
        dims(2) = timeid

#ifdef SINGLE
        status = NF_DEF_VAR(fileid, name, NF_FLOAT, 2, dims(1:2), varid)
#else
        status = NF_DEF_VAR(fileid, name, NF_DOUBLE, 2, dims(1:2), varid)
#endif
        if (status .ne. NF_NOERR) call ncabortname('ncdefreal_height', name, status)
      end if
    else
      dims(1) = timeid

#ifdef SINGLE
      status = NF_DEF_VAR(fileid, name, NF_FLOAT, 1, dims(1), varid)
#else
      status = NF_DEF_VAR(fileid, name, NF_DOUBLE, 1, dims(1), varid)
#endif
      if (status .ne. NF_NOERR) call ncabortname('ncdefreal_height', name, status)
    end if
  end if

  ! Add the attributes
  status = NF_PUT_ATT_TEXT(fileid, varid, 'long_name', LEN(long_name), long_name)
  if (status .ne. NF_NOERR) call ncabortname('ncdefreal_height', name, status)
  status = NF_PUT_ATT_TEXT(fileid, varid, 'units', LEN(units), units)
  if (status .ne. NF_NOERR) call ncabortname('ncdefreal_height', name, status)

end subroutine ncdefreal_height


  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Defines a real output variable using previously defined dimensions. An   =*
  !=  unlimited time dimension is always assumed to exist. The name,           =*
  !=  long_name, and units strings describe the variable.                      =*
  !=  NOTE: If defined, nbins, nlats, and nlons should be the same as the      =*
  !=  values used in nccreate.                                                 =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  fileid    - integer, handle for open netcdf file                         =*
  !=  nbins     - integer, number of vertical levels (midpoint or interface)   =*
  !=  nlats     - integer, number of latitudes                                 =*
  !=  nlons     - integer, number of longitudes                                =*
  !=  name      - character[*], name of field                                  =*
  !=  long_name - character[*], description of field                           =*
  !=  units     - character[*], units of data in the field                     =*
  !-----------------------------------------------------------------------------*
  subroutine ncdefbinreal(fileid, nbins, nlats, nlons, name, long_name, units)
    include 'netcdf.inc'

    ! input:
    integer fileid
    integer nbins
    integer nlats
    integer nlons
    character(LEN=*) name
    character(LEN=*) long_name
    character(LEN=*) units

    ! output: (on converted grid)

    ! local:
    integer status, varid, timeid, binid, latid, lonid
    integer dims(4)
    integer length

    ! Always have an unlimited time dimension
    status = NF_INQ_DIMID(fileid, 'time', timeid)
    if (status .ne. NF_NOERR) call ncabortname('ncdefreal', name, status)

    ! Find the appropriate dimensions and link then to the definition
    ! on this field.
    !
    ! NOTE: To have lons, you must also have lats.
    if (nbins .gt. 1) then
      status = NF_INQ_DIMID(fileid, 'bin', binid)
      if (status .ne. NF_NOERR) call ncabortname('ncdefbinreal', name, status)

      if (nlats .gt. 1) then
        status = NF_INQ_DIMID(fileid, 'lat', latid)
        if (status .ne. NF_NOERR) call ncabortname('ncdefbinreal', name, status)

        if (nlons .gt. 1) then
          status = NF_INQ_DIMID(fileid, 'lon', lonid)
          if (status .ne. NF_NOERR) call ncabortname('ncdefbinreal', name, status)

          dims(1) = lonid
          dims(2) = latid
          dims(3) = binid
          dims(4) = timeid

#ifdef SINGLE
          status = NF_DEF_VAR(fileid, name, NF_FLOAT, 4, dims(1:4), varid)
#else
          status = NF_DEF_VAR(fileid, name, NF_DOUBLE, 4, dims(1:4), varid)
#endif
          if (status .ne. NF_NOERR) call ncabortname('ncdefbinreal', name, status)
        else
          dims(1) = latid
          dims(2) = binid
          dims(3) = timeid

#ifdef SINGLE
          status = NF_DEF_VAR(fileid, name, NF_FLOAT, 3, dims(1:3), varid)
#else
          status = NF_DEF_VAR(fileid, name, NF_DOUBLE, 3, dims(1:3), varid)
#endif
          if (status .ne. NF_NOERR) call ncabortname('ncdefbinreal', name, status)
        end if
      else
        dims(1) = binid
        dims(2) = timeid

#ifdef SINGLE
        status = NF_DEF_VAR(fileid, name, NF_FLOAT, 2, dims(1:2), varid)
#else
        status = NF_DEF_VAR(fileid, name, NF_DOUBLE, 2, dims(1:2), varid)
#endif
        if (status .ne. NF_NOERR) call ncabortname('ncdefbinreal', name, status)
      end if
    else
      if (nlats .gt. 1) then
        status = NF_INQ_DIMID(fileid, 'lat', latid)
        if (status .ne. NF_NOERR) call ncabortname('ncdefbinreal', name, status)

        if (nlons .gt. 1) then
          status = NF_INQ_DIMID(fileid, 'lon', lonid)
          if (status .ne. NF_NOERR) call ncabortname('ncdefbinreal', name, status)

          dims(1) = lonid
          dims(2) = latid
          dims(3) = timeid

#ifdef SINGLE
          status = NF_DEF_VAR(fileid, name, NF_FLOAT, 3, dims(1:3), varid)
#else
          status = NF_DEF_VAR(fileid, name, NF_DOUBLE, 3, dims(1:3), varid)
#endif
          if (status .ne. NF_NOERR) call ncabortname('ncdefbinreal', name, status)
        else
          dims(1) = latid
          dims(2) = timeid

#ifdef SINGLE
          status = NF_DEF_VAR(fileid, name, NF_FLOAT, 2, dims(1:2), varid)
#else
          status = NF_DEF_VAR(fileid, name, NF_DOUBLE, 2, dims(1:2), varid)
#endif
          if (status .ne. NF_NOERR) call ncabortname('ncdefbinreal', name, status)
        end if
      else
        dims(1) = timeid

#ifdef SINGLE
        status = NF_DEF_VAR(fileid, name, NF_FLOAT, 1, dims(1), varid)
#else
        status = NF_DEF_VAR(fileid, name, NF_DOUBLE, 1, dims(1), varid)
#endif
        if (status .ne. NF_NOERR) call ncabortname('ncdefbinreal', name, status)
      end if
    end if

    ! Add the attributes
    status = NF_PUT_ATT_TEXT(fileid, varid, 'long_name', LEN(long_name), long_name)
    if (status .ne. NF_NOERR) call ncabortname('nccreate', name, status)
    status = NF_PUT_ATT_TEXT(fileid, varid, 'units', LEN(units), units)
    if (status .ne. NF_NOERR) call ncabortname('nccreate', name, status)

  end subroutine ncdefbinreal


  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Defines a real output variable using previously defined dimensions. An   =*
  !=  unlimited time dimension is always assumed to exist. The name,           =*
  !=  long_name, and units strings describe the variable.                      =*
  !=  NOTE: If defined, nlevs, nlats, and nlons should be the same as the      =*
  !=  values used in nccreate.                                                 =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  fileid    - integer, handle for open netcdf file                         =*
  !=  nruns     - integer, number of runs                                      =*
  !=  nwaves    - integer, number of wavelengths                               =*
  !=  nbins     - integer, number of bins                                      =*
  !=  name      - character[*], name of field                                  =*
  !=  long_name - character[*], description of field                           =*
  !=  units     - character[*], units of data in the field                     =*
  !-----------------------------------------------------------------------------*
  subroutine ncdefreal_optics(fileid, nwaves, nbins, name, long_name, units)
    include 'netcdf.inc'

    ! input:
    integer fileid
    integer nruns
    integer nwaves
    integer nbins
    character(LEN=*) name
    character(LEN=*) long_name
    character(LEN=*) units

    ! output: (on converted grid)

    ! local:
    integer status, varid, runid, waveid, binid
    integer dims(3)
    integer length

    ! Find the appropriate dimensions and link then to the definition
    ! on this field.
    !
    ! NOTE: To have lons, you must also have lats.
    if (nwaves .gt. 1) then
      status = NF_INQ_DIMID(fileid, 'wave', waveid)
      if (status .ne. NF_NOERR) call ncabortname('ncdefreal_optics', name, status)

      if (nbins .gt. 1) then
        status = NF_INQ_DIMID(fileid, 'bin', binid)
        if (status .ne. NF_NOERR) call ncabortname('ncdefreal_optics', name, status)

        dims(1) = waveid !binid
        dims(2) = binid !waveid

#ifdef SINGLE
        status = NF_DEF_VAR(fileid, name, NF_FLOAT, 2, dims(1:2), varid)
#else
        status = NF_DEF_VAR(fileid, name, NF_DOUBLE, 2, dims(1:2), varid)
#endif
        if (status .ne. NF_NOERR) call ncabortname('ncdefreal_optics', name, status)
      else
        dims(1) = waveid

#ifdef SINGLE
        status = NF_DEF_VAR(fileid, name, NF_FLOAT, 1, dims(1), varid)
#else
        status = NF_DEF_VAR(fileid, name, NF_DOUBLE, 1, dims(1), varid)
#endif
      end if
    else
      if (nbins .gt. 1) then
        status = NF_INQ_DIMID(fileid, 'bin', binid)
        if (status .ne. NF_NOERR) call ncabortname('ncdefreal_optics', name, status)

        dims(1) = binid

#ifdef SINGLE
        status = NF_DEF_VAR(fileid, name, NF_FLOAT, 1, dims(1), varid)
#else
        status = NF_DEF_VAR(fileid, name, NF_DOUBLE, 1, dims(1), varid)
#endif
      endif
    endif

    ! Add the attributes
    status = NF_PUT_ATT_TEXT(fileid, varid, 'long_name', LEN(long_name), long_name)
    if (status .ne. NF_NOERR) call ncabortname('ncdefreal_optics', name, status)
    status = NF_PUT_ATT_TEXT(fileid, varid, 'units', LEN(units), units)
    if (status .ne. NF_NOERR) call ncabortname('ncdefreal_optics', name, status)

  end subroutine ncdefreal_optics

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Writes the dimensions for a previously created NETCDF file.              =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  fileid    - integer, handle for open netcdf file                         =*
  !=  levs      - real(kind=f)[*], levels, midpoints                                   =*
  !=  ilevs     - real(kind=f)[*], levels, interfaces                                  =*
  !=  lats      - real(kind=f)[*], latitudes                                           =*
  !=  lons      - real(kind=f)[*], longitudes                                          =*
  !-----------------------------------------------------------------------------*
  subroutine ncwrtdims(fileid, nlevs, levs, nlats, lats, nlons, lons, nbins, bins)
    include 'netcdf.inc'

    ! input:
    integer fileid, nlevs, nlats, nlons, nbins
    real(kind=f) levs(nlevs)
    real(kind=f) ilevs(nlevs+1)
    real(kind=f) lats(nlats)
    real(kind=f) lons(nlons)
    integer      bins(nbins)

    ! local:
    integer  status, timeid, levid, ilevid, latid, lonid, binid

    ! Define levs, lats, and lons as needed.
    status = NF_INQ_VARID(fileid, 'lev', levid)
    If (status .eq. NF_NOERR) then
#ifdef SINGLE
      status = NF_PUT_VAR_REAL(fileid, levid, levs)
#else
      status = NF_PUT_VAR_DOUBLE(fileid, levid, levs)
#endif
      if (status .ne. NF_NOERR) call ncabortname('ncwrtdims', 'lev', status)
    end if

    status = NF_INQ_VARID(fileid, 'ilev', ilevid)
    If (status .eq. NF_NOERR) then
#ifdef SINGLE
      status = NF_PUT_VAR_REAL(fileid, ilevid, ilevs)
#else
      status = NF_PUT_VAR_DOUBLE(fileid, ilevid, ilevs)
#endif
      if (status .ne. NF_NOERR) call ncabortname('ncwrtdims', 'ilev', status)
    end if

    status = NF_INQ_VARID(fileid, 'lat', latid)
    If (status .eq. NF_NOERR) then
      status = NF_PUT_VAR_REAL(fileid, latid, lats)
      if (status .ne. NF_NOERR) call ncabortname('ncwrtdims', 'lat', status)
    end if
 
    status = NF_INQ_VARID(fileid, 'lon', lonid)
    If (status .eq. NF_NOERR) then
#ifdef SINGLE
      status = NF_PUT_VAR_REAL(fileid, lonid, lons)
#else
      status = NF_PUT_VAR_DOUBLE(fileid, lonid, lons)
#endif
      if (status .ne. NF_NOERR) call ncabortname('ncwrtdims', 'lon', status)
    end if

    status = NF_INQ_VARID(fileid, 'bin', binid)
    If (status .eq. NF_NOERR) then
      status = NF_PUT_VAR_INT(fileid, binid, bins)
      if (status .ne. NF_NOERR) call ncabortname('ncwrtdims', 'bin', status)
    end if

  end subroutine ncwrtdims

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Writes the dimensions for a previously created NETCDF file.              =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  fileid    - integer, handle for open netcdf file                         =*
  !=  levs      - real(kind=f)[*], levels, midpoints                                   =*
  !=  ilevs     - real(kind=f)[*], levels, interfaces                                  =*
  !=  lats      - real(kind=f)[*], latitudes                                           =*
  !=  lons      - real(kind=f)[*], longitudes                                          =*
  !-----------------------------------------------------------------------------*
  subroutine ncwrtdims_height(fileid, nlevs, levs, nlats, lats, nlons, lons, nbins, bins)
  include 'netcdf.inc'

  ! input:
  integer fileid, nlevs, nlats, nlons, nbins
  real(kind=f) levs(nlevs)
  real(kind=f) lats(nlats)
  real(kind=f) lons(nlons)
  integer      bins(nbins)

  ! local:
  integer  status, timeid, levid, latid, lonid, binid

  ! Define levs, lats, and lons as needed.
  status = NF_INQ_VARID(fileid, 'Height', levid)
  If (status .eq. NF_NOERR) then
#ifdef SINGLE
    status = NF_PUT_VAR_REAL(fileid, levid, levs)
#else
    status = NF_PUT_VAR_DOUBLE(fileid, levid, levs)
#endif
    if (status .ne. NF_NOERR) call ncabortname('ncwrtdims_height', 'Height', status)
  end if

  status = NF_INQ_VARID(fileid, 'lat', latid)
  If (status .eq. NF_NOERR) then
    status = NF_PUT_VAR_REAL(fileid, latid, lats)
    if (status .ne. NF_NOERR) call ncabortname('ncwrtdims_height', 'lat', status)
  end if

  status = NF_INQ_VARID(fileid, 'lon', lonid)
  If (status .eq. NF_NOERR) then
#ifdef SINGLE
    status = NF_PUT_VAR_REAL(fileid, lonid, lons)
#else
    status = NF_PUT_VAR_DOUBLE(fileid, lonid, lons)
#endif
    if (status .ne. NF_NOERR) call ncabortname('ncwrtdims_height', 'lon', status)
  end if

  status = NF_INQ_VARID(fileid, 'bin', binid)
  If (status .eq. NF_NOERR) then
    status = NF_PUT_VAR_INT(fileid, binid, bins)
    if (status .ne. NF_NOERR) call ncabortname('ncwrtdims_height', 'bin', status)
  end if

  end subroutine ncwrtdims_height

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Writes the dimensions for a previously created NETCDF file.              =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  fileid    - integer, handle for open netcdf file                         =*
  !=  waves      - real(kind=f)[*], waves                                      =*
  !=  bins      - real(kind=f)[*], bins                                        =*
  !-----------------------------------------------------------------------------*
  subroutine ncwrtdims_optics(fileid, nwaves, waves, nbins, bins)
  include 'netcdf.inc'

  ! input:
  integer fileid, nwaves, nbins
  real(kind=f) waves(nwaves)
  integer      bins(nbins)

  ! local:
  integer  status, waveid, binid

  ! Define waves, and bins as needed.

  status = NF_INQ_VARID(fileid, 'wave', waveid)
  If (status .eq. NF_NOERR) then
#ifdef SINGLE
    status = NF_PUT_VAR_REAL(fileid, waveid, waves)
#else
    status = NF_PUT_VAR_DOUBLE(fileid, waveid, waves)
#endif
    if (status .ne. NF_NOERR) call ncabortname('ncwrtdims_optics', 'wave', status)
  end if
  !write(*,*) "1: wave",waves

  status = NF_INQ_VARID(fileid, 'bin', binid)
  If (status .eq. NF_NOERR) then
    status = NF_PUT_VAR_INT(fileid, binid, bins)
    if (status .ne. NF_NOERR) call ncabortname('ncwrtdims', 'bin', status)
  end if

  end subroutine ncwrtdims_optics

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Writes the time dimension for a previously created NETCDF file. Since    =*
  !=  is an unlimited dimension, it is expanded by one each iteration of the   =*
  !=  time loop.                                                               =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  fileid    - integer, handle for open netcdf file                         =*
  !=  itime     - integer, time index                                          =*
  !=  time      - real(kind=f), time value                                     =*
  !=  date      - real(kind=f), date value                                     =*
  !-----------------------------------------------------------------------------*
  subroutine ncwrttime(fileid, itime, time, date)
    include 'netcdf.inc'

    ! input:
    integer fileid
    integer itime
    real(kind=f)    time
    integer date

    ! local:
    integer status, timeid
    integer start(1)
    integer count(1)
    real(kind=f)    times(1)
    integer dates(1)

    ! Write out a time value.
    status = NF_INQ_VARID(fileid, 'time', timeid)
    if (status .ne. NF_NOERR) call ncabortname('ncwrttime', 'time', status)

    start(1) = itime
    count(1) = 1
    times(1) = time
#ifdef SINGLE
    status = NF_PUT_VARA_REAL(fileid, timeid, start, count, times)
#else
    status = NF_PUT_VARA_DOUBLE(fileid, timeid, start, count, times)
#endif
    if (status .ne. NF_NOERR) call ncabortname('ncwrttime', 'time', status)

    ! Write out a date value.
    status = NF_INQ_VARID(fileid, 'date', timeid)
    if (status .ne. NF_NOERR) call ncabortname('ncwrttime', 'date', status)

    start(1) = itime
    count(1) = 1
    dates(1) = date
 
    status = NF_PUT_VARA_INT(fileid, timeid, start, count, dates)
    if (status .ne. NF_NOERR) call ncabortname('ncwrttime', 'date', status)
  
  end subroutine ncwrttime

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Writes the next value for an real field defined for a previously         =*
!=  created NETCDF file. Indexed to dimensions that aren't used should be    =*
!=  set to -1.                                                               =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  fileid   - integer, handle for open netcdf file                          =*
!=  name     - character[*], name of field                                   =*
!=  nbin     - integer, number of bins (midpoint or interface)               =*
!-----------------------------------------------------------------------------*
  subroutine ncwrtbin(fileid, name, nbins, value)
    include 'netcdf.inc'

    ! input:
    integer fileid
    character(LEN=*) name
    integer nbins
    real(kind=f) value(nbins)

    ! local:
    integer status, varid

    status = NF_INQ_VARID(fileid, name, varid)
    if (status .eq. NF_NOERR) then
#ifdef SINGLE
      status = NF_PUT_VAR_REAL(fileid, varid, value)
#else
      status = NF_PUT_VAR_DOUBLE(fileid, varid, value)
#endif
      if (status .ne. NF_NOERR) call ncabortname('ncwrtdims', 'lev', status)
    end if

  end subroutine ncwrtbin

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Writes the next value for an integer field defined for a previously      =*
  !=  created NETCDF file. Indexed to dimensions that aren't used should be    =*
  !=  set to -1.                                                               =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  fileid   - integer, handle for open netcdf file                          =*
  !=  name     - character[*], name of field                                   =*
  !=  itime    - integer, time index                                           =*
  !=  nlev     - integer, number of levels (midpoint or interface)             =*
  !=  nlat     - integer, number of latitudes                                  =*
  !=  nlon     - integer, number of longitudes                                 =*
  !=  ilev     - integer, level index                                          =*
  !=  ilat     - integer, latitude index                                       =*
  !=  ilon     - integer, longitude index                                      =*
  !=  value    - integer, value                                                =*
  !-----------------------------------------------------------------------------*
  subroutine ncwrtint(fileid, name, itime, nlev, nlat, nlon, ilat, ilon, value)
    include 'netcdf.inc'

    ! input:
    integer fileid
    character(LEN=*) name
    integer itime
    integer nlev
    integer nlat
    integer nlon
    integer ilev
    integer ilat
    integer ilon
    integer value(nlev)

    ! local:
    integer status, varid
    integer start(4)
    integer count(4)
    integer values(1,1,nlev,1)
    integer dimids(NF_MAX_VAR_DIMS)
    character*(NF_MAX_NAME) dimname
    integer i

    ! Define levs, lats, and lons as needed.
    status = NF_INQ_VARID(fileid, name, varid)
    if (status .ne. NF_NOERR) call ncabortname('ncwrtint', name, status)

    values(1,1,:,1) = value(:)

    if (nlev .gt. 1) then
      start(1) = ilon
      start(2) = ilat
      start(3) = 1
      start(4) = itime

      count(1) = 1
      count(2) = 1
      count(3) = nlev
      count(4) = 1

      if (nlat .gt. 1) then
        if (nlon .gt. 1) then
          status = NF_PUT_VARA_INT(fileid, varid, start, count, values)
          if (status .ne. NF_NOERR) call ncabortname('ncwrtint', name, status)
        else
          status = NF_PUT_VARA_INT(fileid, varid, start(2:4), count(2:4), values(1,:,:,:))
          if (status .ne. NF_NOERR) call ncabortname('ncwrtint', name, status)
        end if
      else
        status = NF_PUT_VARA_INT(fileid, varid, start(3:4), count(3:4), values(1,1,:,:))
        if (status .ne. NF_NOERR) call ncabortname('ncwrtint', name, status)
      end if
    else
      start(2) = ilon
      start(3) = ilat
      start(4) = itime

      count(2) = 1
      count(3) = 1
      count(4) = 1

      if (nlat .gt. 1) then
        if (nlon .gt. 1) then
          status = NF_PUT_VARA_INT(fileid, varid, start(2:4), count(2:4), values(:,:,1,:))
          if (status .ne. NF_NOERR) call ncabortname('ncwrtint', name, status)
        else
          status = NF_PUT_VARA_INT(fileid, varid, start(3:4), count(3:4), values(1,:,1,:))
          if (status .ne. NF_NOERR) call ncabortname('ncwrtint', name, status)
        end if
      else
        status = NF_PUT_VARA_INT(fileid, varid, start(4), count(4), values(1,1,1,:))
        if (status .ne. NF_NOERR) call ncabortname('ncwrtint', name, status)
      end if
    end if
  
  end subroutine ncwrtint


  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Writes the next value for a real field defined for a previously          =*
  !=  created NETCDF file. Indexed to dimensions that aren't used should be    =*
  !=  set to -1.                                                               =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  fileid   - integer, handle for open netcdf file                          =*
  !=  name     - character[*], name of field                                   =*
  !=  itime    - integer, time index                                           =*
  !=  nlev     - integer, number of levels (midpoint or interface)             =*
  !=  nlat     - integer, number of latitudes                                  =*
  !=  nlon     - integer, number of longitudes                                 =*
  !=  ilat     - integer, latitude index                                       =*
  !=  ilon     - integer, longitude index                                      =*
  !=  value    - real(kind=f), value(nlev)                                             =*
  !-----------------------------------------------------------------------------*
  subroutine ncwrtreal(fileid, name, itime, nlev, nlat, nlon, ilat, ilon, value)
    include 'netcdf.inc'

    ! input:
    integer fileid
    character(LEN=*) name
    integer itime
    integer nlev
    integer nlat
    integer nlon
    integer ilat
    integer ilon
    real(kind=f) value(nlev)

    ! local:
    integer status, varid
    integer start(4)
    integer count(4)
    real(kind=f)    values(1,1,nlev,1)
    integer dimids(NF_MAX_VAR_DIMS)
    character*(NF_MAX_NAME) dimname
    integer i
  
    ! Define levs, lats, and lons as needed.
    status = NF_INQ_VARID(fileid, name, varid)
    if (status .ne. NF_NOERR) call ncabortname('ncwrtreal', name, status)

    values(1,1,:,1) = value(:)

    if (nlev .gt. 1) then
      start(1) = ilon
      start(2) = ilat
      start(3) = 1
      start(4) = itime

      count(1) = 1
      count(2) = 1
      count(3) = nlev
      count(4) = 1

      if (nlat .gt. 1) then
        if (nlon .gt. 1) then
#ifdef SINGLE
          status = NF_PUT_VARA_REAL(fileid, varid, start, count, values)
#else
          status = NF_PUT_VARA_DOUBLE(fileid, varid, start, count, values)
#endif
          if (status .ne. NF_NOERR) call ncabortname('ncwrtreal', name, status)
        else
#ifdef SINGLE
          status = NF_PUT_VARA_REAL(fileid, varid, start(2:4), count(2:4), values(1,:,:,:))
#else
          status = NF_PUT_VARA_DOUBLE(fileid, varid, start(2:4), count(2:4), values(1,:,:,:))
#endif
          if (status .ne. NF_NOERR) call ncabortname('ncwrtreal', name, status)
        end if
      else
#ifdef SINGLE
        status = NF_PUT_VARA_REAL(fileid, varid, start(3:4), count(3:4), values(1,1,:,:))
#else
        status = NF_PUT_VARA_DOUBLE(fileid, varid, start(3:4), count(3:4), values(1,1,:,:))
#endif
        if (status .ne. NF_NOERR) call ncabortname('ncwrtreal', name, status)
      end if
    else
      start(2) = ilon
      start(3) = ilat
      start(4) = itime

      count(2) = 1
      count(3) = 1
      count(4) = 1

      if (nlat .gt. 1) then
        if (nlon .gt. 1) then
#ifdef SINGLE
          status = NF_PUT_VARA_REAL(fileid, varid, start(2:4), count(2:4), values(:,:,1,:))
#else
          status = NF_PUT_VARA_DOUBLE(fileid, varid, start(2:4), count(2:4), values(:,:,1,:))
#endif
          if (status .ne. NF_NOERR) call ncabortname('ncwrtreal', name, status)
        else
#ifdef SINGLE
          status = NF_PUT_VARA_REAL(fileid, varid, start(3:4), count(3:4), values(1,:,1,:))
#else
          status = NF_PUT_VARA_DOUBLE(fileid, varid, start(3:4), count(3:4), values(1,:,1,:))
#endif
          if (status .ne. NF_NOERR) call ncabortname('ncwrtreal', name, status)
        end if
      else
#ifdef SINGLE
        status = NF_PUT_VARA_REAL(fileid, varid, start(4), count(4), values(1,1,1,:))
#else
        status = NF_PUT_VARA_DOUBLE(fileid, varid, start(4), count(4), values(1,1,1,:))
#endif
        if (status .ne. NF_NOERR) call ncabortname('ncwrtreal', name, status)
      end if
    end if
  
  end subroutine ncwrtreal

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Writes the next value for a real field defined for a previously          =*
  !=  created NETCDF file. Indexed to dimensions that aren't used should be    =*
  !=  set to -1.                                                               =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  fileid   - integer, handle for open netcdf file                          =*
  !=  name     - character[*], name of field                                   =*
  !=  itime    - integer, time index                                           =*
  !=  nlev     - integer, number of levels (midpoint or interface)             =*
  !=  nlat     - integer, number of latitudes                                  =*
  !=  nlon     - integer, number of longitudes                                 =*
  !=  ilat     - integer, latitude index                                       =*
  !=  ilon     - integer, longitude index                                      =*
  !=  value    - real(kind=f), value(nlev)                                             =*
  !-----------------------------------------------------------------------------*
  subroutine ncwrtbinreal(fileid, name, itime, nbin, nlat, nlon, ilat, ilon, value)
    include 'netcdf.inc'

    ! input:
    integer fileid
    character(LEN=*) name
    integer itime
    integer nbin
    integer nlat
    integer nlon
    integer ilat
    integer ilon
    real(kind=f) value(nbin)

    ! local:
    integer status, varid
    integer start(4)
    integer count(4)
    real(kind=f)    values(1,1,nbin,1)
    integer dimids(NF_MAX_VAR_DIMS)
    character*(NF_MAX_NAME) dimname
    integer i

    ! Define bins, lats, and lons as needed.
    status = NF_INQ_VARID(fileid, name, varid)
    if (status .ne. NF_NOERR) call ncabortname('ncwrtbinreal', name, status)

    values(1,1,:,1) = value(:)

    if (nbin .gt. 1) then
      start(1) = ilon
      start(2) = ilat
      start(3) = 1
      start(4) = itime

      count(1) = 1
      count(2) = 1
      count(3) = nbin
      count(4) = 1

      if (nlat .gt. 1) then
        if (nlon .gt. 1) then
#ifdef SINGLE
          status = NF_PUT_VARA_REAL(fileid, varid, start, count, values)
#else
          status = NF_PUT_VARA_DOUBLE(fileid, varid, start, count, values)
#endif
          if (status .ne. NF_NOERR) call ncabortname('ncwrtbinreal', name, status)
        else
#ifdef SINGLE
          status = NF_PUT_VARA_REAL(fileid, varid, start(2:4), count(2:4), values(1,:,:,:))
#else
          status = NF_PUT_VARA_DOUBLE(fileid, varid, start(2:4), count(2:4), values(1,:,:,:))
#endif
          if (status .ne. NF_NOERR) call ncabortname('ncwrtbinreal', name, status)
        end if
      else
#ifdef SINGLE
        status = NF_PUT_VARA_REAL(fileid, varid, start(3:4), count(3:4), values(1,1,:,:))
#else
        status = NF_PUT_VARA_DOUBLE(fileid, varid, start(3:4), count(3:4), values(1,1,:,:))
#endif
        if (status .ne. NF_NOERR) call ncabortname('ncwrtbinreal', name, status)
      end if
    else
      start(2) = ilon
      start(3) = ilat
      start(4) = itime

      count(2) = 1
      count(3) = 1
      count(4) = 1

      if (nlat .gt. 1) then
        if (nlon .gt. 1) then
#ifdef SINGLE
          status = NF_PUT_VARA_REAL(fileid, varid, start(2:4), count(2:4), values(:,:,1,:))
#else
          status = NF_PUT_VARA_DOUBLE(fileid, varid, start(2:4), count(2:4), values(:,:,1,:))
#endif
          if (status .ne. NF_NOERR) call ncabortname('ncwrtbinreal', name, status)
        else
#ifdef SINGLE
          status = NF_PUT_VARA_REAL(fileid, varid, start(3:4), count(3:4), values(1,:,1,:))
#else
          status = NF_PUT_VARA_DOUBLE(fileid, varid, start(3:4), count(3:4), values(1,:,1,:))
#endif
          if (status .ne. NF_NOERR) call ncabortname('ncwrtbinreal', name, status)
        end if
      else
#ifdef SINGLE
        status = NF_PUT_VARA_REAL(fileid, varid, start(4), count(4), values(1,1,1,:))
#else
        status = NF_PUT_VARA_DOUBLE(fileid, varid, start(4), count(4), values(1,1,1,:))
#endif
        if (status .ne. NF_NOERR) call ncabortname('ncwrtbinreal', name, status)
      end if
    end if

  end subroutine ncwrtbinreal

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Writes the next value for a real field defined for a previously          =*
  !=  created NETCDF file. Indexed to dimensions that aren't used should be    =*
  !=  set to -1.                                                               =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  fileid   - integer, handle for open netcdf file                          =*
  !=  name     - character[*], name of field                                   =*
  !=  itime    - integer, time index                                           =*
  !=  nrun     - integer, number of levels (midpoint or interface)             =*
  !=  nwave    - integer, number of latitudes                                  =*
  !=  nbin     - integer, number of longitudes                                 =*
  !=  value    - real(kind=f), value(nlev)                                     =*
  !-----------------------------------------------------------------------------*
  subroutine ncwrtreal_optics(fileid, name, nwave, nbin, value)
    include 'netcdf.inc'

    ! input:
    integer fileid
    character(LEN=*) name
    integer nwave
    integer nbin
    real(kind=f) value(nwave,nbin)

    ! local:
    integer status, varid
    integer start(2)
    integer count(2)
    real(kind=f)    values(nwave,nbin)
    integer i

    ! Define levs, lats, and lons as needed.
    status = NF_INQ_VARID(fileid, name, varid)
    if (status .ne. NF_NOERR) call ncabortname('ncwrtreal_optics', name, status)

    values(:,:) = value(:,:)

    if (nwave .gt. 1) then
      if (nbin .gt. 1) then
        start(1) = 1
        start(2) = 1

        count(1) = nwave !nbin
        count(2) = nbin !nwave
        
#ifdef SINGLE
        status = NF_PUT_VARA_REAL(fileid, varid, start, count, values)
#else
        status = NF_PUT_VARA_DOUBLE(fileid, varid, start, count, values)
#endif
        if (status .ne. NF_NOERR) call ncabortname('ncwrtreal_optics', name, status)
      else
        start(2) = 1

        count(2) = nwave

#ifdef SINGLE
        status = NF_PUT_VARA_REAL(fileid, varid, start(2), count(2), values)
#else
        status = NF_PUT_VARA_DOUBLE(fileid, varid, start(2), count(2), values)
#endif
        if (status .ne. NF_NOERR) call ncabortname('ncwrtreal_optics', name, status)
      end if
    else
      if (nbin .gt. 1) then
        start(2) = 1

        count(2) = nbin

#ifdef SINGLE
        status = NF_PUT_VARA_REAL(fileid, varid, start(2), count(2), values)
#else
        status = NF_PUT_VARA_DOUBLE(fileid, varid, start(2), count(2), values)
#endif
      end if
    end if


  end subroutine ncwrtreal_optics

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Closes a NETCDF file that was previously created or opened.              =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  fileid   - integer, handle for open netcdf file                          =*
  !-----------------------------------------------------------------------------*
  subroutine ncclose(fileid)
    include 'netcdf.inc'

    ! input:
    integer fileid

    ! local:
    integer  status

    ! Define levs, lats, and lons as needed.
    status = NF_CLOSE(fileid)
    if (status .ne. NF_NOERR) call ncabort('ncclose', status)

  end subroutine ncclose


  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Opens an existing NETCDF file for read only access.                      =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  filename   - character[*], name of the file                              =*
  !=  omode      - integer, NETCDF open mode (0=READ ONLY, 1=READ/WRITE)       =*
  !-----------------------------------------------------------------------------*
  !=  RETURNS:                                                                 =*
  !=  fileid  - integer, handle for open netcdf file                           =*
  !-----------------------------------------------------------------------------*
  function ncopen(filename, omode)
    include 'netcdf.inc'

    ! input:
    character(LEN=*) filename
    integer omode

    ! function value:
    integer ncopen

    ! local:
    integer status
    integer fileid
 
    ! Open an exisitng NETCDF file.
    status = NF_OPEN(filename, omode, fileid)
    if (status .ne. NF_NOERR) call ncabortname('ncopen', filename, status)
 
    ncopen = fileid
 
  end function ncopen


  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Flushes an existing NETCDF file to disk.                                 =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  fileid  - integer, handle for open netcdf file                           =*
  !-----------------------------------------------------------------------------*
  subroutine ncsync(fileid)
    include 'netcdf.inc'

    ! input:
    integer fileid

    ! local:
    integer status
 
    ! Flush an exisitng NETCDF file.
    status = NF_SYNC(fileid)
    if (status .ne. NF_NOERR) call ncabort('ncsync', status)

  end subroutine ncsync


  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Reads the size of a dimension in a NETCDF file.                          =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  fileid    - integer, handle for open netcdf file                         =*
  !=  name      - character[*], name of the dimension                          =*
  !-----------------------------------------------------------------------------*
  !=  RETURNS:                                                                 =*
  !=  ndim      - integer, the size of the dimension                           =*
  !-----------------------------------------------------------------------------*
  function ncrdndim(fileid, name)
    include 'netcdf.inc'

    ! input:
    integer fileid
    character(LEN=*) name

    ! function value:
    integer ncrdndim

    ! local:
    integer status
    integer dimid
    integer ndim
 
    ! Open an exisitng NETCDF file.
    status = NF_INQ_DIMID(fileid, name, dimid)
    if (status .ne. NF_NOERR) call ncabort('ncrdndim', status)
  
    status = NF_INQ_DIMLEN(fileid, dimid, ndim)
    if (status .ne. NF_NOERR) call ncabort('ncrdndim', status)

    ncrdndim = ndim
 
  end function ncrdndim

      
  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Reads an integer value of a variable from a NETCDF file.                 =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  fileid    - integer, handle for open netcdf file                         =*
  !=  name      - character[*], name of the dimension                          =*
  !=  offset    - integer, offset of first entry to read                       =*
  !=  length    - integer, number of entries in variable                       =*
  !=  value     - integer(length), the values of the field                     =*
  !-----------------------------------------------------------------------------*
  subroutine ncrdint1d(fileid, name, offset, length, value)
    include 'netcdf.inc'

    ! input:
    integer fileid
    character(LEN=*) name
    integer offset
    integer length

    ! function value:
    integer value(length)

    ! local:
    integer status, varid
    integer dimids(NF_MAX_VAR_DIMS)
    integer dimlens(NF_MAX_VAR_DIMS)
    integer ndims
    integer start(1)
    integer count(1)

    ! Lookup variable.
    status = NF_INQ_VARID(fileid, name, varid)
    if (status .ne. NF_NOERR) call ncabortname('ncrdint', name, status)
  
    ! See what dimensions it has.      
    status = NF_INQ_VARNDIMS(fileid, varid, ndims)
    if (status .ne. NF_NOERR) call ncabortname('ncrdint', name, status)
  
    if (ndims .ne. 1) then
      write(*,*) 'ncrdint: ERROR - Fields are required to have 1 dimension.'
      call ncabortname(trim('ncrdint'), name, NF_NOERR)
    end if
  
    ! Get the actual size of each dimension.
    status = NF_INQ_VARDIMID(fileid, varid, dimids)
    if (status .ne. NF_NOERR) call ncabortname('ncrdint', name, status)
  
    status = NF_INQ_DIMLEN(fileid, dimids(1), dimlens(1))
    if (status .ne. NF_NOERR) call ncabortname('ncrdint', name, status)

    ! Check the dimensions in the file to make sure they match the dimensions
    ! defined for the variable.
    if ((offset + length - 1) .gt. dimlens(1)) then
      write(*,*) 'ncrdreal: ERROR - offset + length too big.'
      call ncabortname(trim('ncrdint'), name, NF_NOERR)
    end if
  
    start(1) = offset
    count(1) = length

    status = NF_GET_VARA_INT(fileid, varid, start(1), count(1), value(:))
    if (status .ne. NF_NOERR) call ncabortname('ncrdint', name, status)
      
  end subroutine ncrdint1d
      
      
  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Reads integer values of a variable for a particular time from a NETCDF   =*
  !=  file.                                                                    =*
  !=  NOTE: The minimum value for itime, nlat, and nlon is 1.                  =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  fileid    - integer, handle for open netcdf file                         =*
  !=  name      - character[*], name of the dimension                          =*
  !=  itime     - integer, time index                                          =*
  !=  nlat      - integer, number of latitudes                                 =*
  !=  nlon      - integer, number of longitudes                                =*
  !=  value     - integer(nlat, nlon), the values of the field                 =*
  !-----------------------------------------------------------------------------*
  subroutine ncrdint2d(fileid, name, itime, nlat, nlon, value)
    include 'netcdf.inc'

    ! input:
    integer fileid
    character(LEN=*) name
    integer itime
    integer nlat
    integer nlon

    ! output:
    integer value(nlat, nlon)

    ! local:
    integer status, varid, ndims, idim
    integer dimids(NF_MAX_VAR_DIMS)
    integer dimlens(NF_MAX_VAR_DIMS)
    integer start(3)
    integer count(3)
    integer values(nlon, nlat, 1)

    ! Define levs, lats, and lons as needed.
    status = NF_INQ_VARID(fileid, name, varid)
    if (status .ne. NF_NOERR) call ncabort('ncrdint2d', status)
  
    ! See what dimensions it has. It is required to have a time and a lev
    ! dimension, but it may not have lat and lon.      
    status = NF_INQ_VARNDIMS(fileid, varid, ndims)
    if (status .ne. NF_NOERR) call ncabort('ncrdint2d', status)
  
    if (ndims .lt. 1) then
      write(*,*) 'ncrdint2d: ERROR - 2D Fields are required to have at least 1 dimension.'
      call ncabortname('ncrdint2d', name, NF_NOERR)
    end if
  
    ! Get the actual size of each dimension.
    status = NF_INQ_VARDIMID(fileid, varid, dimids)
    if (status .ne. NF_NOERR) call ncabortname('ncrdint2d', name, status)
  
    do idim = 1, ndims
      status = NF_INQ_DIMLEN(fileid, dimids(idim), dimlens(idim))
      if (status .ne. NF_NOERR) call ncabortname('ncrdint2d', name, status)
    end do

    ! Check the dimensions in the file to make sure they match the dimensions
    ! defined for the output variable.
    if (itime .gt. dimlens(ndims)) then
      write(*,*) 'ncrdint2d: ERROR - Time index bigger than dimension.'
      call ncabortname('ncrdint2d', name, NF_NOERR)
    end if
    start(3) = itime
    count(3) = 1

    if (ndims .ge. 1) then        
      if (ndims .eq. 1) then
        status = NF_GET_VARA_INT(fileid, varid, start(3), count(3), values(1,1,:))
        if (status .ne. NF_NOERR) call ncabortname('ncrdint2d', name, status)
      else if (ndims .gt. 1) then
        start(2) = 1
        count(2) = nlat

        if (ndims .eq. 2) then
          status = NF_GET_VARA_INT(fileid, varid, start(2:3), count(2:3), values(1,:,:))
          if (status .ne. NF_NOERR) call ncabortname('ncrdint2d', name, status)
        else if (ndims .eq. 3) then
          start(1) = 1
          count(1) = nlon

          status = NF_GET_VARA_INT(fileid, varid, start, count, values(:,:,:))
          if (status .ne. NF_NOERR) call ncabortname('ncrdint2d', name, status)
        else
          write(*,*) 'ncrdint2d: ERROR - Variable has too many dimensions.'
          call ncabortname('ncrdint2d', name, NF_NOERR)
        end if
      end if
    end if
  
    value(:,:) = values(:,:,1)
      
  end subroutine ncrdint2d
            

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Reads integer values of a variable for a particular time from a NETCDF   =*
  !=  file.                                                                    =*
  !=  NOTE: The minimum value for itime, nlev, nlat, and nlon is 1.            =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  fileid    - integer, handle for open netcdf file                         =*
  !=  name      - character[*], name of the dimension                          =*
  !=  itime     - integer, time index                                          =*
  !=  nlev      - integer, number of levels                                    =*
  !=  nlat      - integer, number of latitudes                                 =*
  !=  nlon      - integer, number of longitudes                                =*
  !=  value     - integer(nlon, nlat, nlev), the values of the field           =*
  !-----------------------------------------------------------------------------*
  subroutine ncrdint3d(fileid, name, itime, nlev, nlat, nlon, value)
    include 'netcdf.inc'

    ! input:
    integer fileid
    character(LEN=*) name
    integer itime
    integer nlev
    integer nlat
    integer nlon

    ! function value:
    integer value(nlon, nlat, nlev)

   ! local:
    integer status, varid, ndims, idim
    integer dimids(NF_MAX_VAR_DIMS)
    integer dimlens(NF_MAX_VAR_DIMS)
    integer start(4)
    integer count(4)
    integer values(nlon, nlat, nlev, 1)

    ! Define levs, lats, and lons as needed.
    status = NF_INQ_VARID(fileid, name, varid)
    if (status .ne. NF_NOERR) call ncabortname('ncrdint3d', name, status)
  
    ! See what dimensions it has. It is required to have a time and a lev
    ! dimension, but it may not have lat and lon.      
    status = NF_INQ_VARNDIMS(fileid, varid, ndims)
    if (status .ne. NF_NOERR) call ncabortname('ncrdint3d', name, status)
  
    if (ndims .lt. 2) then
      write(*,*) 'ncrdint3d: ERROR - Fields are required to have at least 2 dimensions.'
      call ncabort(trim('ncrdint3d'), NF_NOERR)
    end if
  
    ! Get the actual size of each dimension.
    status = NF_INQ_VARDIMID(fileid, varid, dimids)
    if (status .ne. NF_NOERR) call ncabortname('ncrdint3d', name, status)
  
    do idim = 1, ndims
      status = NF_INQ_DIMLEN(fileid, dimids(idim), dimlens(idim))
      if (status .ne. NF_NOERR) call ncabortname('ncrdint3d', name, status)
    end do

    ! Check the dimensions in the file to make sure they match the dimensions
    ! defined for the output variable.
    if (itime .gt. dimlens(ndims)) then
      write(*,*) 'ncrdint3d: ERROR - Time index bigger than dimension.'
      call ncabortname('ncrdint3d', name, NF_NOERR)
    end if
    start(4) = itime
    count(4) = 1

    if (dimlens(ndims-1) .ne. nlev) then
      write(*,*) 'ncrdint3d: ERROR - Level dimensions do not match.'
      call ncabortname('ncrdint3d', name, NF_NOERR)
    end if
  
    if (ndims .ge. 2) then
      start(3) = 1
      count(3) = nlev
    
      if (ndims .eq. 2) then
        status = NF_GET_VARA_INT(fileid, varid, start(3:4), count(3:4), values(1,1,:,:))
        if (status .ne. NF_NOERR) call ncabortname('ncrdint3d', name, status)
      else if (ndims .gt. 2) then
        start(2) = 1
        count(2) = nlat

        if (ndims .eq. 3) then
          status = NF_GET_VARA_INT(fileid, varid, start(2:4), count(2:4), values(1,:,:,:))
          if (status .ne. NF_NOERR) call ncabortname('ncrdint3d', name, status)
        else if (ndims .eq. 4) then
          start(1) = 1
          count(1) = nlon

          status = NF_GET_VARA_INT(fileid, varid, start, count, values(:,:,:,:))
          if (status .ne. NF_NOERR) call ncabortname('ncrdint3d', name, status)
        else
          write(*,*) 'ncrdint3d: ERROR - Variable has too many dimensions.'
          call ncabortname('ncrdint3d', name, NF_NOERR)
        end if
      end if
    end if
  
    value(:,:,:) = values(:,:,:,1)
      
  end subroutine ncrdint3d
      
      
  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Reads an integer value of a variable from a NETCDF file.                 =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  fileid    - integer, handle for open netcdf file                         =*
  !=  name      - character[*], name of the dimension                          =*
  !=  offset    - integer, offset of first entry to read                       =*
  !=  length    - integer, number entries to read                              =*
  !=  value     - real(kind=f)(length), the values of the field                        =*
  !-----------------------------------------------------------------------------*
  subroutine ncrdreal1d(fileid, name, offset, length, value)
    include 'netcdf.inc'

    ! input:
    integer fileid
    character(LEN=*) name
    integer offset
    integer length

    ! output:
    real(kind=f) value(length)

    ! local:
    integer status, varid, idim
    integer dimids(NF_MAX_VAR_DIMS)
    integer dimlens(NF_MAX_VAR_DIMS)
    integer ndims
    integer start(1)
    integer count(1)

    ! Lookup variable.
    status = NF_INQ_VARID(fileid, name, varid)
    if (status .ne. NF_NOERR) call ncabortname('ncrdreal', name, status)
  
    ! See what dimensions it has.      
    status = NF_INQ_VARNDIMS(fileid, varid, ndims)
    if (status .ne. NF_NOERR) call ncabortname('ncrdreal', name, status)
  
    if (ndims .ne. 1) then
      write(*,*) 'ncrdreal: ERROR - Fields are required to have 1 dimension.'
      call ncabortname('ncrdreal', name, NF_NOERR)
    end if
  
    ! Get the actual size of each dimension.
    status = NF_INQ_VARDIMID(fileid, varid, dimids)
    if (status .ne. NF_NOERR) call ncabortname('ncrdreal', name, status)
  
    status = NF_INQ_DIMLEN(fileid, dimids(1), dimlens(1))
    if (status .ne. NF_NOERR) call ncabortname('ncrdreal', name, status)

    ! Check the dimensions in the file to make sure they match the dimensions
    ! defined for the variable.
    if ((offset + length - 1) .gt. dimlens(1)) then
      write(*,*) 'ncrdreal: ERROR - offset + length too big : ', offset, length
      call ncabortname('ncrdreal', name, NF_NOERR)
    end if
  
    start(1) = offset
    count(1) = length

#ifdef SINGLE
    status = NF_GET_VARA_REAL(fileid, varid, start(1), count(1), value(:))
#else
    status = NF_GET_VARA_DOUBLE(fileid, varid, start(1), count(1), value(:))
#endif
    if (status .ne. NF_NOERR) call ncabortname('ncrdreal', name, status)
      
  end subroutine ncrdreal1d
      
      
  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Reads real values of a variable for a particular time from a NETCDF      =*
  !=  file.                                                                    =*
  !=  NOTE: The minimum value for itime, nlat, and nlon is 1.                  =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  fileid    - integer, handle for open netcdf file                         =*
  !=  name      - character[*], name of the dimension                          =*
  !=  itime     - integer, time index                                          =*
  !=  nlat      - integer, number of latitudes                                 =*
  !=  nlon      - integer, number of longitudes                                =*
  !=  value     - real(kind=f)(nlon, nlat), the values of the field                    =*
  !-----------------------------------------------------------------------------*
  subroutine ncrdreal2d(fileid, name, itime, nlat, nlon, value)
    include 'netcdf.inc'

    ! input:
    integer fileid
    character(LEN=*) name
    integer itime
    integer nlat
    integer nlon

    ! output:
    real(kind=f) value(nlon, nlat)

    ! local:
    integer status, varid, ndims, idim
    integer dimids(NF_MAX_VAR_DIMS)
    integer dimlens(NF_MAX_VAR_DIMS)
    integer start(3)
    integer count(3)
    real(kind=f) values(nlon, nlat, 1)

    ! Define levs, lats, and lons as needed.
    status = NF_INQ_VARID(fileid, name, varid)
    if (status .ne. NF_NOERR) call ncabortname('ncrdreal2d', name, status)
  
    ! See what dimensions it has. It is required to have a time and a lev
    ! dimension, but it may not have lat and lon.      
    status = NF_INQ_VARNDIMS(fileid, varid, ndims)
    if (status .ne. NF_NOERR) call ncabortname('ncrdreal2d', name, status)
  
    if (ndims .lt. 1) then
      write(*,*) 'ncrdreal2d: ERROR - 2D Fields are required to have at least 1 dimension.'
      call ncabortname('ncrdreal2d', name, NF_NOERR)
    end if
  
    ! Get the actual size of each dimension.
    status = NF_INQ_VARDIMID(fileid, varid, dimids)
    if (status .ne. NF_NOERR) call ncabortname('ncrdreal2d', name, status)
  
    do idim = 1, ndims
      status = NF_INQ_DIMLEN(fileid, dimids(idim), dimlens(idim))
      if (status .ne. NF_NOERR) call ncabortname('ncrdreal2d', name, status)
    end do

    ! Check the dimensions in the file to make sure they match the dimensions
    ! defined for the output variable.
    if (itime .gt. dimlens(ndims)) then
      write(*,*) 'ncrdreal2d: ERROR - Time index bigger than dimension.'
      call ncabortname('ncrdreal2d', name, NF_NOERR)
    end if
    start(3) = itime
    count(3) = 1

    if (ndims .ge. 1) then        
      if (ndims .eq. 1) then
#ifdef SINGLE
        status = NF_GET_VARA_REAL(fileid, varid, start(3), count(3), values(1,1,:))
#else
        status = NF_GET_VARA_DOUBLE(fileid, varid, start(3), count(3), values(1,1,:))
#endif
        if (status .ne. NF_NOERR) call ncabortname('ncrdreal2d', name, status)
      else if (ndims .gt. 1) then
        start(2) = 1
        count(2) = nlat

        if (ndims .eq. 2) then
#ifdef SINGLE
          status = NF_GET_VARA_REAL(fileid, varid, start(2:3), count(2:3), values(1,:,:))
#else
          status = NF_GET_VARA_DOUBLE(fileid, varid, start(2:3), count(2:3), values(1,:,:))
#endif
          if (status .ne. NF_NOERR) call ncabortname('ncrdreal2d', name, status)
        else if (ndims .eq. 3) then
          start(1) = 1
          count(1) = nlon

#ifdef SINGLE
          status = NF_GET_VARA_REAL(fileid, varid, start, count, values(:,:,:))
#else
          status = NF_GET_VARA_DOUBLE(fileid, varid, start, count, values(:,:,:))
#endif
          if (status .ne. NF_NOERR) call ncabortname('ncrdreal2d', name, status)
        else
          write(*,*) 'ncrdreal2d: ERROR - Variable has too many dimensions.'
          call ncabortname('ncrdreal2d', name, NF_NOERR)
        end if
      end if
    end if
  
    value(:,:) = values(:,:,1)
      
  end subroutine ncrdreal2d
            

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Reads real values of a variable for a particular time from a NETCDF      =*
  !=  file.                                                                    =*
  !=  NOTE: The minimum value for itime, nlev, nlat, and nlon is 1.            =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  fileid    - integer, handle for open netcdf file                         =*
  !=  name      - character[*], name of the dimension                          =*
  !=  itime     - integer, time index                                          =*
  !=  nlev      - integer, number of levels                                    =*
  !=  nlat      - integer, number of latitudes                                 =*
  !=  nlon      - integer, number of longitudes                                =*
  !=  value     - real(kind=f)(nlon, nlat, nlev), the values of the field              =*
  !-----------------------------------------------------------------------------*
  subroutine ncrdreal3d(fileid, name, itime, nlev, nlat, nlon, value)
    include 'netcdf.inc'

    ! input:
    integer fileid
    character(LEN=*) name
    integer itime
    integer nlev
    integer nlat
    integer nlon

    ! output:
    real(kind=f) value(nlon, nlat, nlev)

    ! local:
    integer status, varid, ndims, idim
    integer dimids(NF_MAX_VAR_DIMS)
    integer dimlens(NF_MAX_VAR_DIMS)
    integer start(4)
    integer count(4)
    real(kind=f) values(nlon, nlat, nlev, 1)

    ! Define levs, lats, and lons as needed.
    status = NF_INQ_VARID(fileid, name, varid)
    if (status .ne. NF_NOERR) call ncabortname('ncrdreal3d', name, status)
  
    ! See what dimensions it has. It is required to have a time and a lev
    ! dimension, but it may not have lat and lon.      
    status = NF_INQ_VARNDIMS(fileid, varid, ndims)
    if (status .ne. NF_NOERR) call ncabortname('ncrdreal3d', name, status)
  
    if (ndims .lt. 2) then
      write(*,*) 'ncrdreal3d: ERROR - Fields are required to have at least 2 dimensions.'
      call ncabortname('ncrdreal3d', name, NF_NOERR)
    end if
  
    ! Get the actual size of each dimension.
    status = NF_INQ_VARDIMID(fileid, varid, dimids)
    if (status .ne. NF_NOERR) call ncabortname('ncrdreal3d', name, status)
    do idim = 1, ndims
      status = NF_INQ_DIMLEN(fileid, dimids(idim), dimlens(idim))
      if (status .ne. NF_NOERR) call ncabortname('ncrdreal3d', name, status)
    end do

    ! Check the dimensions in the file to make sure they match the dimensions
    ! defined for the output variable.
    if (itime .gt. dimlens(ndims)) then
      write(*,*) 'ncrdreal3d: ERROR - Time index bigger than dimension.'
      call ncabortname('ncrdreal3d', name, NF_NOERR)
    end if
    start(4) = itime
    count(4) = 1

    if (dimlens(ndims-1) .ne. nlev) then
      write(*,*) 'ncrdreal3d: ERROR - Level dimensions do not match.'
      call ncabortname('ncrdreal3d', name, NF_NOERR)
    end if
  
    if (ndims .ge. 2) then
      start(3) = 1
      count(3) = nlev
    
      if (ndims .eq. 2) then
#ifdef SINGLE
        status = NF_GET_VARA_REAL(fileid, varid, start(3:4), count(3:4), values(1,1,:,:))
#else
        status = NF_GET_VARA_DOUBLE(fileid, varid, start(3:4), count(3:4), values(1,1,:,:))
#endif
        if (status .ne. NF_NOERR) call ncabortname('ncrdreal3d', name, status)
      else if (ndims .gt. 2) then
        start(2) = 1
        count(2) = nlat

        if (ndims .eq. 3) then
#ifdef SINGLE
          status = NF_GET_VARA_REAL(fileid, varid, start(2:4), count(2:4), values(1,:,:,:))
#else
          status = NF_GET_VARA_DOUBLE(fileid, varid, start(2:4), count(2:4), values(1,:,:,:))
#endif
          if (status .ne. NF_NOERR) call ncabortname('ncrdreal3d', name, status)
        else if (ndims .eq. 4) then
          start(1) = 1
          count(1) = nlon

#ifdef SINGLE
          status = NF_GET_VARA_REAL(fileid, varid, start, count, values(:,:,:,:))
#else
          status = NF_GET_VARA_DOUBLE(fileid, varid, start, count, values(:,:,:,:))
#endif
          if (status .ne. NF_NOERR) call ncabortname('ncrdreal3d', name, status)
        else
          write(*,*) 'ncrdreal3d: ERROR - Variable has too many dimensions.'
          call ncabortname('ncrdreal3d', name, NF_NOERR)
        end if
      end if
    end if
  
    value(:,:,:) = values(:,:,:,1)
      
  end subroutine ncrdreal3d


  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Reads real values of a variable for a particular time from a NETCDF      =*
  !=  file.                                                                    =*
  !=  NOTE: The minimum value for itime, nlev, nlat, and nlon is 1.            =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  fileid    - integer, handle for open netcdf file                         =*
  !=  name      - character[*], name of the dimension                          =*
  !=  itime     - integer, time index                                          =*
  !=  nlev      - integer, number of levels                                    =*
  !=  nlat      - integer, number of latitudes                                 =*
  !=  nlon      - integer, number of longitudes                                =*
  !=  nwave     - integer, number of wavelengths                               =*
  !=  value     - real(kind=f)(nlon, nlat, nlev), the values of the field              =*
  !-----------------------------------------------------------------------------*

  subroutine ncrdreal3dw(fileid, name, itime, nlev, nlat, nlon, nwave, value)
    include 'netcdf.inc'

    ! input:
    integer fileid
    character(LEN=*) name
    integer itime
    integer nlev
    integer nlat
    integer nlon
    integer nwave

    ! output:
    real(kind=f) value(nwave, nlon, nlat, nlev)

    ! local:
    integer status, varid, ndims, idim
    integer dimids(NF_MAX_VAR_DIMS)
    integer dimlens(NF_MAX_VAR_DIMS)
    integer start(5)
    integer count(5)
    real(kind=f) values(nwave, nlon, nlat, nlev, 1)

    ! Define levs, lats, and lons as needed.
    status = NF_INQ_VARID(fileid, name, varid)
    if (status .ne. NF_NOERR) call ncabortname('ncrdreal3dw', name, status)
  
    ! See what dimensions it has. It is required to have a time and a lev
    ! dimension, but it may not have lat and lon.      
    status = NF_INQ_VARNDIMS(fileid, varid, ndims)
    if (status .ne. NF_NOERR) call ncabortname('ncrdreal3dw', name, status)
  
    if (ndims .lt. 3) then
      write(*,*) 'ncrdreal3dw: ERROR - Fields are required to have at least 3 dimensions.'
      call ncabortname('ncrdreal3dw', name, NF_NOERR)
    end if
  
    ! Get the actual size of each dimension.
    status = NF_INQ_VARDIMID(fileid, varid, dimids)
    if (status .ne. NF_NOERR) call ncabortname('ncrdreal3dw', name, status)
    do idim = 1, ndims
      status = NF_INQ_DIMLEN(fileid, dimids(idim), dimlens(idim))
      if (status .ne. NF_NOERR) call ncabortname('ncrdreal3dw', name, status)
    end do

    ! Check the dimensions in the file to make sure they match the dimensions
    ! defined for the output variable.
    if (itime .gt. dimlens(ndims)) then
      write(*,*) 'ncrdreal3dw: ERROR - Time index bigger than dimension.'
      call ncabortname('ncrdreal3dw', name, NF_NOERR)
    end if
    start(5) = itime
    count(5) = 1

    if (dimlens(ndims-1) .ne. nlev) then
      write(*,*) 'ncrdreal3dw: ERROR - Level dimensions do not match.'
      call ncabortname('ncrdreal3dw', name, NF_NOERR)
    end if
  
    if (ndims .ge. 3) then
      start(4) = 1
      count(4) = nlev
    
      if (ndims .eq. 3) then
        start(3) = 1
        count(3) = nwave
      
#ifdef SINGLE
        status = NF_GET_VARA_REAL(fileid, varid, start(3:5), count(3:5), values(:,1,1,:,:))
#else
        status = NF_GET_VARA_DOUBLE(fileid, varid, start(3:5), count(3:5), values(:,1,1,:,:))
#endif
        if (status .ne. NF_NOERR) call ncabortname('ncrdreal3dw', name, status)
      else if (ndims .gt. 3) then
        start(3) = 1
        count(3) = nlat

        if (ndims .eq. 4) then
          start(2) = 1
          count(2) = nwave

#ifdef SINGLE
          status = NF_GET_VARA_REAL(fileid, varid, start(2:5), count(2:5), values(:,1,:,:,:))
#else
          status = NF_GET_VARA_DOUBLE(fileid, varid, start(2:5), count(2:5), values(:,1,:,:,:))
#endif
          if (status .ne. NF_NOERR) call ncabortname('ncrdreal3dw', name, status)
        else if (ndims .eq. 5) then
          start(2) = 1
          count(2) = nlon
          start(1) = 1
          count(1) = nwave

#ifdef SINGLE
          status = NF_GET_VARA_REAL(fileid, varid, start, count, values(:,:,:,:,:))
#else
          status = NF_GET_VARA_DOUBLE(fileid, varid, start, count, values(:,:,:,:,:))
#endif
          if (status .ne. NF_NOERR) call ncabortname('ncrdreal3dw', name, status)
        else
          write(*,*) 'ncrdreal3dw: ERROR - Variable has too many dimensions.'
          call ncabortname('ncrdreal3dw', name, NF_NOERR)
        end if
      end if
    end if
  
    value(:,:,:,:) = values(:,:,:,:,1)
      
  end subroutine ncrdreal3dw

end module
