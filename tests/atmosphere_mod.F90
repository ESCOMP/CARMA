!! Public domain code to calculate the US Standard Atmosphere pressure and temperature.
!! This underlying routine, Atmosphere was downloaded from the internet at:
!!
!!  http://www.pdas.com/programs/atmos.f90
!!
!! This module wraps the routine in an interface similar to the one used by the CARMA
!! module, which allows columns, vectors of columns, and arrays of columns via the
!! method GetStandardAtmosphere().
module atmosphere_mod
  ! types
  use carma_precision_mod

  implicit none
  
  private

  ! NOTE: From US Standard Atmosphere 1976, the standard pressure and
  ! temperature at sea-level are:
	!   P0 = 1.013250e5 (Pa)
	!   T0 = 288.15 (K)
  real(kind=f), public, parameter       :: P0 = 1.013250e5_f        !! Standard sea-level pressure
  real(kind=f), public, parameter       :: T0 = 288.15_f            !! Standard sea-level temperature

	interface GetStandardAtmosphere
		module procedure GetStandardAtmosphere_1D
		module procedure GetStandardAtmosphere_2D
		module procedure GetStandardAtmosphere_3D
	end interface
	
	public GetStandardAtmosphere

  contains

		subroutine GetStandardAtmosphere_1D(z, p, t)
			real(kind=f), intent(in)              :: z(:)   !! Geometric Altitude (m)
			real(kind=f), optional, intent(out)   :: p(:)   !! pressure (Pa)
			real(kind=f), optional, intent(out)   :: t(:)   !! temperature (K)
			
			! Local variables
			real    :: sigma
			real    :: delta
			real    :: theta
			integer :: i
			integer :: NZ
			
			NZ = size(z, 1)
	
			do i = 1, NZ
				
				! Get the scaling of the pressure and temperature at the altitude.
				call Atmosphere(real(z(i) / 1000._f), sigma, delta, theta)     !! Convert from m -> km
				
				if (present(p)) p(i) = p0 * delta
				if (present(t)) t(i) = T0 * theta
			end do
			
			return
		end subroutine
	
		subroutine GetStandardAtmosphere_2D(z, p, t)
			real(kind=f), intent(in)              :: z(:, :)   !! Geometric Altitude (m)
			real(kind=f), optional, intent(out)   :: p(:, :)   !! pressure (Pa)
			real(kind=f), optional, intent(out)   :: t(:, :)   !! temperature (K)
			
			! Local variables
			real    :: sigma
			real    :: delta
			real    :: theta
			integer :: i, j
			integer :: NY, NZ
			
			NY = size(z, 1)
			NZ = size(z, 2)
			
			do i = 1, NY
				do j = 1, NZ
				
					! Get the scaling of the pressure and temperature at the altitude.
					call Atmosphere(real(z(i, j) / 1000._f), sigma, delta, theta)     !! Convert from m -> km
				
					if (present(p)) p(i, j) = p0 * delta
					if (present(t)) t(i, j) = T0 * theta
				end do
			end do
			
			return
		end subroutine
	
		subroutine GetStandardAtmosphere_3D(z, p, t)
			real(kind=f), intent(in)              :: z(:, :, :)   !! Geometric Altitude (m)
			real(kind=f), optional, intent(out)   :: p(:, :, :)   !! pressure (Pa)
			real(kind=f), optional, intent(out)   :: t(:, :, :)   !! temperature (K)
			
			! Local variables
			real    :: sigma
			real    :: delta
			real    :: theta
			integer :: i, j, k
			integer :: NX, NY, NZ
			
			NX = size(z, 1)
			NY = size(z, 2)
			NZ = size(z, 3)
			
			do i = 1, NX
				do j = 1, NY
					do k = 1, NZ
					
						! Get the scaling of the pressure and temperature at the altitude.
						call Atmosphere(real(z(i, j, k) / 1000._f), sigma, delta, theta)     !! Convert from m -> km
					
						if (present(p)) p(i, j, k) = p0 * delta
						if (present(t)) t(i, j, k) = T0 * theta
					end do
				end do
		  end do
			
			return
		end subroutine
	
		!+
		SUBROUTINE Atmosphere(alt, sigma, delta, theta)
		!   -------------------------------------------------------------------------
		! PURPOSE - Compute the properties of the 1976 standard atmosphere to 86 km.
		! AUTHOR - Ralph Carmichael, Public Domain Aeronautical Software
		! NOTE - If alt > 86, the values returned will not be correct, but they will
		!   not be too far removed from the correct values for density.
		!   The reference document does not use the terms pressure and temperature
		!   above 86 km.
		IMPLICIT NONE
		!============================================================================
		!     A R G U M E N T S                                                     |
		!============================================================================
			real,INTENT(IN)::  alt        ! geometric altitude, km.
			real,INTENT(OUT):: sigma      ! density/sea-level standard density
			real,INTENT(OUT):: delta      ! pressure/sea-level standard pressure
			real,INTENT(OUT):: theta      ! temperature/sea-level standard temperature
		!============================================================================
		!     L O C A L   C O N S T A N T S                                         |
		!============================================================================
			real,PARAMETER:: REARTH = 6369.0                 ! radius of the Earth (km)
			real,PARAMETER:: GMR = 34.163195                     ! hydrostatic constant
			INTEGER,PARAMETER:: NTAB=8       ! number of entries in the defining tables
		!============================================================================
		!     L O C A L   V A R I A B L E S                                         |
		!============================================================================
			INTEGER:: i,j,k                                                  ! counters
			real:: h                                       ! geopotential altitude (km)
			real:: tgrad, tbase      ! temperature gradient and base temp of this layer
			real:: tlocal                                           ! local temperature
			real:: deltah                             ! height above base of this layer
		!============================================================================
		!     L O C A L   A R R A Y S   ( 1 9 7 6   S T D.  A T M O S P H E R E )   |
		!============================================================================
			real,DIMENSION(NTAB),PARAMETER:: htab= &
															(/0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 71.0, 84.852/)
			real,DIMENSION(NTAB),PARAMETER:: ttab= &
							(/288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.946/)
			real,DIMENSION(NTAB),PARAMETER:: ptab= &
									 (/1.0, 2.233611E-1, 5.403295E-2, 8.5666784E-3, 1.0945601E-3, &
																				 6.6063531E-4, 3.9046834E-5, 3.68501E-6/)
			real,DIMENSION(NTAB),PARAMETER:: gtab= &
																		(/-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0, 0.0/)
		!----------------------------------------------------------------------------
			h=alt*REARTH/(alt+REARTH)      ! convert geometric to geopotential altitude
		
			i=1
			j=NTAB                                       ! setting up for binary search
			DO
				k=(i+j)/2                                              ! integer division
				IF (h < htab(k)) THEN
					j=k
				ELSE
					i=k
				END IF
				IF (j <= i+1) EXIT
			END DO
		
			tgrad=gtab(i)                                     ! i will be in 1...NTAB-1
			tbase=ttab(i)
			deltah=h-htab(i)
			tlocal=tbase+tgrad*deltah
			theta=tlocal/ttab(1)                                    ! temperature ratio
		
			IF (tgrad == 0.0) THEN                                     ! pressure ratio
				delta=ptab(i)*EXP(-GMR*deltah/tbase)
			ELSE
				delta=ptab(i)*(tbase/tlocal)**(GMR/tgrad)
			END IF
		
			sigma=delta/theta                                           ! density ratio
			RETURN
		END Subroutine Atmosphere   ! -----------------------------------------------
end module
