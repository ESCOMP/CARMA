! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"


module sulfate_utils
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod
 
  implicit none
 ! Declare the public methods.
  public wtpct_tabaz
  public sulfate_density
  public sulfate_surf_tens
  
  real(kind=f), public:: dnwtp(46), dnc0(46), dnc1(46)
   
  data dnwtp/0,1,5,10,20,25,30,35,40,41,45,50,53,55,56,60,65,66,70, &
     & 72,73,74,75,76,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,  &
     & 94,95,96,97,98,100/
     
   data dnc0/1,1.13185,1.17171,1.22164,1.3219,1.37209,1.42185,1.4705,&
     & 1.51767,1.52731,1.56584,1.61834,1.65191,1.6752,1.68708,1.7356,   &
     & 1.7997,1.81271,1.86696,1.89491,1.9092,1.92395,1.93904,1.95438,   &
     & 1.98574,2.00151,2.01703,2.03234,2.04716,2.06082,2.07363,2.08461, &
     & 2.09386,2.10143,2.10764,2.11283,2.11671,2.11938,2.12125,2.1219,  &
     & 2.12723,2.12654,2.12621,2.12561,2.12494,2.12093/
     
   data dnc1/0,-0.000435022,-0.000479481,-0.000531558,-0.000622448,  &
     & -0.000660866,-0.000693492,-0.000718251,-0.000732869,-0.000735755,&
     & -0.000744294,-0.000761493,-0.000774238,-0.00078392,-0.000788939, &
     & -0.00080946,-0.000839848,-0.000845825,-0.000874337,-0.000890074, &
     & -0.00089873,-0.000908778,-0.000920012,-0.000932184,-0.000959514, &
     & -0.000974043,-0.000988264,-0.00100258,-0.00101634,-0.00102762,   &
     & -0.00103757,-0.00104337,-0.00104563,-0.00104458,-0.00104144,     &
     & -0.00103719,-0.00103089,-0.00102262,-0.00101355,-0.00100249,     &
     & -0.00100934,-0.000998299,-0.000990961,-0.000985845,-0.000984529, &
     & -0.000989315/  
contains

  !!  This function calculates the weight % H2SO4 composition of 
  !!  sulfate aerosol, using Tabazadeh et. al. (GRL, 1931, 1997).
  !!  Rated for T=185-260K, activity=0.01-1.0
  !!
  !!  Argument list input:   
  !!    temp = Temperature (K)
  !!    h2o_mass = water vapor mass concentration (g/cm3)
  !!    h2o_vp = water eq. vaper pressure (dynes/cm2)
  !!
  !!  Output:
  !!    Weight % H2SO4 in H2O/H2SO4 particle
  !!
  !!  Include global constants and variables (BK=Boltzman constant,
  !!   AVG=Avogadro's constant)
  !!
  !! @author Jason English
  !! @ version Apr-2010
  function wtpct_tabaz(carma, temp, h2o_mass, h2o_vp, rc)
   
    real(kind=f)                         :: wtpct_tabaz
    type(carma_type), intent(in)         :: carma     !! the carma object
    real(kind=f), intent(in)             :: temp      !! temperature [K]
    real(kind=f), intent(in)             :: h2o_mass  !! water vapor mass concentration (g/cm3)
    real(kind=f), intent(in)             :: h2o_vp    !! water eq. vaper pressure (dynes/cm2) 
    integer, intent(inout)               :: rc        !! return code, negative indicates failure
      
    !  Declare variables for this routine only
    real(kind=f)     :: atab1,btab1,ctab1,dtab1,atab2,btab2,ctab2,dtab2
    real(kind=f)     :: h2o_num, p_h2o, vp_h2o
    real(kind=f)     :: contl, conth, contt, conwtp
    real(kind=f)     :: activ 
         
    ! Get number density of water (/cm3) from mass concentration (g/cm3)
    h2o_num=h2o_mass*AVG/gwtmol(1)

    !  Get partial pressure of water (dynes/cm2) from concentration (/cm3)
    ! Ideal gas law: P=nkT
    p_h2o=h2o_num*bk*temp

    !  Convert from dynes/cm2 to mb (hPa)
    p_h2o=p_h2o/1000.0     ! partial pressure
    vp_h2o=h2o_vp/1000.0   ! eq. vp

    !  Prevent a NaN calculation  
    !  In the upper thermosphere p_h2o can be very low and vp_h2o can be very high
    if (p_h2o.lt.1.e-10 .and. vp_h2o.gt.0.) p_h2o=1.e-10   
   
    !  Activity = water pp in mb / water eq. vp over pure water in mb
    activ = p_h2o/vp_h2o
 
    if (activ.lt.0.05) then
      activ = max(activ,1.e-6)    ! restrict minimum activity
      atab1 	= 12.37208932	
      btab1 	= -0.16125516114
      ctab1 	= -30.490657554
      dtab1 	= -2.1133114241	
      atab2 	= 13.455394705	
      btab2 	= -0.1921312255	
      ctab2 	= -34.285174607	
      dtab2 	= -1.7620073078
    elseif (activ.ge.0.05.and.activ.le.0.85) then
      atab1 	= 11.820654354	
      btab1 	= -0.20786404244
      ctab1 	= -4.807306373
      dtab1 	= -5.1727540348	
      atab2 	= 12.891938068	
      btab2 	= -0.23233847708
      ctab2 	= -6.4261237757	
      dtab2 	= -4.9005471319	
    elseif (activ.gt.0.85) then
      activ = min(activ,1.)      ! restrict maximum activity
      atab1 	= -180.06541028	
      btab1 	= -0.38601102592
      ctab1 	= -93.317846778
      dtab1 	= 273.88132245
      atab2 	= -176.95814097	
      btab2 	= -0.36257048154
      ctab2 	= -90.469744201
      dtab2 	= 267.45509988
    else
      if (do_print) write(LUNOPRT,*) 'invalid activity: activity,pp,vp=',activ, p_h2o
      rc = RC_ERROR
      return
    endif

    contl = atab1*(activ**btab1)+ctab1*activ+dtab1
    conth = atab2*(activ**btab2)+ctab2*activ+dtab2
      
!    temp=min(max(temp,185.),260.)  ! T range 185-260
    
    contt = contl + (conth-contl) * ((temp -190.)/70.)
    conwtp = (contt*98.) + 1000.

    wtpct_tabaz = (100.*contt*98.)/conwtp
    wtpct_tabaz = min(max(wtpct_tabaz,1.),100.) ! restrict between 1 and 100 %
      
    !  Note: restricting activity to 1.e-6 minimum allows for a maximum of
    !  98.5 wtpct at T=650K, 95.8 wtpct at T=300K, and 90.9 wtpct at 180K.
  
    return
  end function wtpct_tabaz
           
  function sulfate_density(carma, WTP,TEMP, rc)
    real(kind=f)                         :: sulfate_density
    type(carma_type), intent(in)         :: carma   !! the carma object
    real(kind=f), intent(in)             :: WTP     !! weight percent
    real(kind=f), intent(in)             :: TEMP    !! temperature 
    integer, intent(inout)               :: rc      !! return code, negative indicates failure
    
    ! Local declerations
    integer           :: I
    real(kind=f)      :: DEN1, DEN2
    real(kind=f)      :: FRAC

    if (wtp .LT. 0.0 .OR. wtp .GT. 100.0) then
      if (do_print) write(LUNOPRT,*)'sulfate_density: Illegal value for wtp:',wtp
      rc = RC_ERROR
      return
    endif

    I=1

    DO WHILE (WTP .GT. DNWTP(I))
     I=I+1
    END DO

    DEN2=DNC0(I)+DNC1(I)*TEMP

    IF (I.EQ.1 .OR. WTP.EQ.DNWTP(I)) THEN
      sulfate_density=DEN2
      RETURN
    ENDIF

    DEN1=DNC0(I-1)+DNC1(I-1)*TEMP
    FRAC=(DNWTP(I)-WTP)/(DNWTP(I)-DNWTP(I-1))
    sulfate_density=DEN1*FRAC+DEN2*(1.0-FRAC)

    RETURN
  END FUNCTION sulfate_density

  !!  Calculates surface tension (erg/cm2) of sulfate of 
  !!  different compositions as a linear function of temperature.
  !!
  !!  Argument list input:
  !!     WTP = aerosol composition in weight % H2SO4 (0-100)
  !!     TEMP = temperature in Kelvin
  !!
  !!  Output:
  !!     sulfate_density (g/cm3) [function name]
  !!
  !!  This function requires setup_sulfate_density to be run
  !!  first to read in the density coefficients DNC0 and DNC1
  !!  and the tabulated weight percents DNWTP.
  !!
  !! @author Mike Mills
  !! @version Jul-2001
  function sulfate_surf_tens(carma, wtp,temp, rc)
    real(kind=f)                         :: sulfate_surf_tens
    type(carma_type), intent(in)         :: carma   !! the carma object
    real(kind=f), intent(in)             :: wtp     !! weight percent
    real(kind=f), intent(in)             :: temp    !! temperature 
    integer, intent(inout)               :: rc      !! return code, negative indicates failure
    
    ! Local declerations
    integer           :: i  
    real(kind=f)      :: sig1, sig2
    real(kind=f)      :: frac
    real(kind=f)      :: stwtp(15), stc0(15), stc1(15)
    
    data stwtp/0, 23.8141, 38.0279, 40.6856, 45.335, 52.9305, 56.2735, &
       & 59.8557, 66.2364, 73.103, 79.432, 85.9195, 91.7444, 97.6687, 100/
    
    data stc0/117.564, 103.303, 101.796, 100.42, 98.4993, 91.8866,     &
       & 88.3033, 86.5546, 84.471, 81.2939, 79.3556, 75.608, 70.0777,  &
       & 63.7412, 61.4591 /
  
    data stc1/-0.153641, -0.0982007, -0.0872379, -0.0818509,           &
       & -0.0746702, -0.0522399, -0.0407773, -0.0357946, -0.0317062,   &
       & -0.025825, -0.0267212, -0.0269204, -0.0276187, -0.0302094,    &
       & -0.0303081 /
       
    if (wtp .lt. 0.0_f .OR. wtp .gt. 100.0_f) then
      if (do_print) write(LUNOPRT,*)'sulfate_surf_tens: Illegal value for wtp:',wtp
      if (do_print) write(LUNOPRT,*)'sulfate_surf_tens: temp=',temp
      rc = RC_ERROR
      return
    endif
  
    i=1
  
    do while (wtp.gt.stwtp(i))
     i=i+1
    end do
  
    sig2=stc0(i)+stc1(i)*temp
  
    if (I.eq.1 .or. wtp.eq.stwtp(i)) then
      sulfate_surf_tens=sig2
      return
    end if
  
    sig1=stc0(i-1)+stc1(i-1)*temp
    frac=(stwtp(i)-wtp)/(stwtp(i)-stwtp(i-1))
    sulfate_surf_tens=sig1*frac+sig2*(1.0_f-frac)
  
    return
  end function sulfate_surf_tens
            
end module sulfate_utils
