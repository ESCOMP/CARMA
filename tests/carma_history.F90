module history

  contains

		subroutine write_history(filename, xc, yc, zc, p, t, pc, gc, rc)
			use netcdf
			implicit none
		
			! This is the name of the data file we will create.
			character (len = *), parameter :: FILE_NAME = "pres_temp_4D.nc"
			integer :: ncid
		
			! We are writing 4D data, a 2 x 6 x 12 lvl-lat-lon grid, with 2
			! timesteps of data.
			integer, parameter :: NDIMS = 4, NRECS = 2
			integer, parameter :: NLVLS = 2, NLATS = 6, NLONS = 12
			integer :: lvl_dimid, lon_dimid, lat_dimid, rec_dimid
		
			! The start and count arrays will tell the netCDF library where to
			! write our data.
			integer :: start(NDIMS), count(NDIMS)
		
			! These program variables hold the latitudes and longitudes.
			real    :: lats(NLATS), lons(NLONS)
			integer :: lon_varid, lat_varid
		
			! We will create two netCDF variables, one each for temperature and
			! pressure fields.
			character (len = *), parameter :: PRES_NAME="pressure"
			character (len = *), parameter :: TEMP_NAME="temperature"
			integer :: pres_varid, temp_varid
			integer :: dimids(NDIMS)
		
			! We recommend that each variable carry a "units" attribute.
		
			! Program variables to hold the data we will write out. We will only
			! need enough space to hold one timestep of data; one record.
			real :: pres_out(NLONS, NLATS, NLVLS)
			real :: temp_out(NLONS, NLATS, NLVLS)
			real, parameter :: SAMPLE_PRESSURE = 900.0
			real, parameter :: SAMPLE_TEMP = 9.0
		
			! Use these to construct some latitude and longitude data for this
			! example.
			real, parameter :: START_LAT = 25.0, START_LON = -125.0
		
			! Loop indices
			integer :: lvl, lat, lon, rec, i
		
			! Create pretend data. If this wasn't an example program, we would
			! have some real data to write, for example, model output.
			do lat = 1, NLATS
				 lats(lat) = START_LAT + (lat - 1) * 5.0
			end do
			do lon = 1, NLONS
				 lons(lon) = START_LON + (lon - 1) * 5.0
			end do
			i = 0
			do lvl = 1, NLVLS
				 do lat = 1, NLATS
						do lon = 1, NLONS
							 pres_out(lon, lat, lvl) = SAMPLE_PRESSURE + i
							 temp_out(lon, lat, lvl) = SAMPLE_TEMP + i
							 i = i + 1
						end do
				 end do
			end do
			
			return
		end subroutine	

		HISTORY_Create(history, filename, carma)
			character (len = *), parameter :: LVL_NAME = "level"
			character (len = *), parameter :: LAT_NAME = "latitude"
			character (len = *), parameter :: LON_NAME = "longitude"
			character (len = *), parameter :: REC_NAME = "time"
		
			character (len = *), parameter :: UNITS = "units"
			character (len = *), parameter :: PRES_UNITS = "Pa"
			character (len = *), parameter :: TEMP_UNITS = "K"
			character (len = *), parameter :: LAT_UNITS = "degrees_north"
			character (len = *), parameter :: LON_UNITS = "degrees_east"
	
      integer ::  NGAS
      integer ::  NELEM
      integer ::  NBIN
      
      call CARMA_GetDims(NGROUP, NELEM, NBIN, NSOLUTE, NGAS, NX, NY, NZ)

			! Create the file. 
			call check(nf90_create(filename, nf90_clobber, history%ncid))
		
			! Define the dimensions. The record dimension is defined to have
			! unlimited length - it can grow as needed. In this example it is
			! the time dimension.
			call check(nf90_def_dim(history%ncid, LVL_NAME, NLVLS, history%lvl_dimid))
			call check(nf90_def_dim(history%ncid, LAT_NAME, NLATS, history%lat_dimid))
			call check(nf90_def_dim(history%ncid, LON_NAME, NLONS, history%lon_dimid))
			call check(nf90_def_dim(history%ncid, REC_NAME, NF90_UNLIMITED, history%rec_dimid))
	
			! Define the coordinate variables. We will only define coordinate
			! variables for lat and lon.  Ordinarily we would need to provide
			! an array of dimension IDs for each variable's dimensions, but
			! since coordinate variables only have one dimension, we can
			! simply provide the address of that dimension ID (lat_dimid) and
			! similarly for (lon_dimid).
			call check(nf90_def_var(history%ncid, LAT_NAME, NF90_REAL, history%lat_dimid, history%lat_varid))
			call check(nf90_def_var(history%ncid, LON_NAME, NF90_REAL, history%lon_dimid, history%lon_varid))
	
			! Assign units attributes to coordinate variables.
			call check(nf90_put_att(history%ncid, lat_varid, UNITS, LAT_UNITS))
			call check(nf90_put_att(history%ncid, lon_varid, UNITS, LON_UNITS))
	
			! The dimids array is used to pass the dimids of the dimensions of
			! the netCDF variables. Both of the netCDF variables we are creating
			! share the same four dimensions. In Fortran, the unlimited
			! dimension must come last on the list of dimids.
			dimids = (/ lon_dimid, lat_dimid, lvl_dimid, rec_dimid /)
	
			! Define the netCDF variables for the pressure and temperature data.
			call check(nf90_def_var(history%ncid, PRES_NAME, NF90_REAL, dimids, history%pres_varid))
			call check(nf90_def_var(history%ncid, TEMP_NAME, NF90_REAL, dimids, history%temp_varid))
	
			! Assign units attributes to the netCDF variables.
			call check(nf90_put_att(history%ncid, history%pres_varid, UNITS, PRES_UNITS))
			call check(nf90_put_att(history%ncid, history%temp_varid, UNITS, TEMP_UNITS))
		
			! End define mode.
			call check(nf90_enddef(history%ncid))
		
			! Write the coordinate variable data. This will put the latitudes
			! and longitudes of our data grid into the netCDF file.
			call check(nf90_put_var(history%ncid, history%lat_varid, lats))
			call check(nf90_put_var(history%ncid, history%lon_varid, lons))
		end
  
    subroutine History_WriteStep(history, carma)
      real(kind=f)            :: mmr(history%NX, history%NY, history%NZ)
      real(kind=f), pointer   :: mmrptr(:,:,:)
      
      ! These settings tell netcdf to write one timestep of data. (The
      ! setting of start(4) inside the loop below tells netCDF which
      ! timestep to write.)
      count = (/ history%NX, history%NY, history%NZ, 1 /)
      countl = (/ history%NX, history%NY, history%NZ+1, 1 /)
      start = (/ 1, 1, 1, 1 /)

      ! Write the data.
      start(4) = history%time
      
      ! Write out the atmospheric state.
      call check(nf90_put_var(history%ncid, history%p_varid, p, start=start, count=count))
      call check(nf90_put_var(history%ncid, history%t_varid, t, start=start, count=count))
      call check(nf90_put_var(history%ncid, history%z_varid, z, start=start, count=count))
      call check(nf90_put_var(history%ncid, history%pl_varid, pl, start=start, count=countl))
      call check(nf90_put_var(history%ncid, history%zl_varid, zl, start=start, count=countl))
                              
      ! Write out the bins
      mmrptr => mmr
      do ielem = 1, history%NELEM
        do ibin = 1, history%NBIN
          call CARMA_GetBin(carma, ielem, ibin, mmrptr, rc)
          call check(nf90_put_var(history%ncid, history%bin_varid(ielem,ibin), mmr, start=start, count=count))
        end do
      end do
      
      ! Write out the gases
      do igas = 1, history%NGAS
        call CARMA_GetGas(carma, igas, mmrptr, rc)
        call check(nf90_put_var(history%ncid, history%gas_varid(igas), mmr, start=start, count=count))
      end do
      
      ! Write out statistics?
      
      history%time = history%time+1
    end
    
		!! Close the history file. This causes netCDF to flush all buffers and make
		!! sure your data are really written to disk.
    subroutine History_Close(history))
      call check(nf90_close(history%ncid))
      return
    end
  
    subroutine check(status)
      integer, intent ( in) :: status
    
      if(status /= nf90_noerr) then 
        print *, trim(nf90_strerror(status))
        stop "Stopped"
      end if
    end  
end
