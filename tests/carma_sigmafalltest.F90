!! This code is to demonstrate the CARMA sedimentation routines
!! using a constant fall velocity. Sigma coordinates are used
!! by the test.
!!
!! Upon execution, a text file (carma_sigmafalltest.txt) is generated.
!! The text file can be read with the IDL procedure read_sigmafalltest.pro.
!!
!! @author Peter Colarco (based on Chuck Bardeen's code)
!! @version Feb-2009


program carma_sigmafalltest
  implicit none

  write(*,*) "Sedimentation Test (Sigma Coordinates)"

  call test_sedimentation_sigma()  
  
  write(*,*) "Done"
end program


subroutine test_sedimentation_sigma()
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
  integer, parameter    :: NZ           = 72
  integer, parameter    :: NZP1         = NZ+1
  integer, parameter    :: NELEM        = 1
  integer, parameter    :: NBIN         = 8
  integer, parameter    :: NGROUP       = 1
  integer, parameter    :: NSOLUTE      = 0
  integer, parameter    :: NGAS         = 0
  integer, parameter    :: NWAVE        = 0
  integer, parameter    :: nstep        = 100*6

  real(kind=f), parameter   :: dtime  = 1000._f
  real(kind=f), parameter   :: deltax = 100._f
  real(kind=f), parameter   :: deltay = 100._f
  
  integer, parameter        :: I_DUST       = 1
  integer, parameter        :: I_ICE        = 2

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
  
  real(kind=f), allocatable, target  :: mmr(:,:,:,:,:)
  
  real(kind=f), allocatable          :: lat(:,:)
  real(kind=f), allocatable          :: lon(:,:)
  
  integer               :: i
  integer               :: ix
  integer               :: iy
  integer               :: ixy
  integer               :: istep
  integer               :: ielem
  integer               :: ibin
  integer               :: ithread
  integer, parameter    :: lun = 42

  real(kind=f)          :: time
  real(kind=f)          :: rmin, rmrat, rho
  logical               :: do_explised = .false.
!  logical               :: do_explised = .true.
  real(kind=f)          :: vf_const = 2.0_f
!  real(kind=f)          :: vf_const = 0.0_f
  
  integer               :: omp_get_num_threads, omp_get_max_threads, &
                           omp_get_thread_num
  

  real(kind=f)          :: a72(73), b72(73), t72(72), ze(73)
  real(kind=f)          :: hyai66(67), hybi66(67), hyam66(66), hybm66(66)
  real(kind=f)          :: hyai125(126), hybi125(126), hyam125(125), hybm125(125)
  real(kind=f)          :: a(NZP1), b(NZP1), dz(NZ), zm(NZP1), rhoa(NZ)
  real(kind=f), parameter :: ps = 98139.8  ! GEOS-5, Omaha, NE, 20090101

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
      
  ! The WACCM 66 level hybrid coefficients.
  data hyai66 /&
    4.5005e-09, 7.4201e-09, 1.22337e-08, 2.017e-08, 3.32545e-08, &
    5.48275e-08, 9.0398e-08, 1.4904e-07, 2.4572e-07, 4.05125e-07, 6.6794e-07, &
    1.101265e-06, 1.81565e-06, 2.9935e-06, 4.963e-06, 8.150651e-06, &
    1.3477e-05, 2.2319e-05, 3.67965e-05, 6.0665e-05, 9.91565e-05, 0.00015739, &
    0.00023885, 0.0003452, 0.000475135, 0.000631805, 0.000829155, 0.00108274, &
    0.00140685, 0.00181885, 0.0023398, 0.00299505, 0.0038147, 0.00483445, &
    0.00609635, 0.00764935, 0.0095501, 0.011864, 0.0146655, 0.018038, &
    0.0220755, 0.0268825, 0.0325735, 0.039273, 0.0471145, 0.0562405, &
    0.0668005, 0.0789485, 0.07731271, 0.07590131, 0.07424086, 0.07228743, &
    0.06998932, 0.06728574, 0.06410509, 0.06036322, 0.05596111, 0.05078225, &
    0.0446896, 0.03752191, 0.02908949, 0.02084739, 0.01334443, 0.00708499, &
    0.00252136, 0, 0 /

  data hybi66 /&
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    0, 0.01505309, 0.03276228, 0.05359622, 0.07810627, 0.1069411, 0.1408637, &
    0.180772, 0.227722, 0.2829562, 0.3479364, 0.4243822, 0.5143168, &
    0.6201202, 0.7235355, 0.8176768, 0.8962153, 0.9534761, 0.9851122, 1 /

  data hyam66 /&
    5.9603e-09, 9.8269e-09, 1.620185e-08, 2.671225e-08, 4.4041e-08, &
    7.261275e-08, 1.19719e-07, 1.9738e-07, 3.254225e-07, 5.365325e-07, &
    8.846025e-07, 1.4584575e-06, 2.404575e-06, 3.97825e-06, 6.5568255e-06, &
    1.08138255e-05, 1.7898e-05, 2.955775e-05, 4.873075e-05, 7.991075e-05, &
    0.00012827325, 0.00019812, 0.000292025, 0.0004101675, 0.00055347, &
    0.00073048, 0.0009559475, 0.001244795, 0.00161285, 0.002079325, &
    0.002667425, 0.003404875, 0.004324575, 0.0054654, 0.00687285, &
    0.008599725, 0.01070705, 0.01326475, 0.01635175, 0.02005675, 0.024479, &
    0.029728, 0.03592325, 0.04319375, 0.0516775, 0.0615205, 0.0728745, &
    0.078130605, 0.07660701, 0.075071085, 0.073264145, 0.071138375, &
    0.06863753, 0.065695415, 0.062234155, 0.058162165, 0.05337168, &
    0.047735925, 0.041105755, 0.0333057, 0.02496844, 0.01709591, 0.01021471, &
    0.004803175, 0.00126068, 0 /

  data hybm66 /&
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    0.007526545, 0.023907685, 0.04317925, 0.065851245, 0.092523685, &
    0.1239024, 0.16081785, 0.204247, 0.2553391, 0.3154463, 0.3861593, &
    0.4693495, 0.5672185, 0.67182785, 0.77060615, 0.85694605, 0.9248457, &
    0.96929415, 0.9925561 /

  ! The WACCM 125 level hybrid coefficients.
  data hyai125 /&    
    4.5005e-09, 7.4247466e-09, 1.2236616e-08, 2.0143422e-08, &
    3.3006333e-08, 5.317267e-08, 8.160345e-08, 1.1600713e-07, 1.5343022e-07, &
    1.9154761e-07, 2.2843057e-07, 2.648671e-07, 3.0222616e-07, 3.3995556e-07, &
    3.7754795e-07, 4.1533548e-07, 4.5535556e-07, 4.9759899e-07, &
    5.4204611e-07, 5.8866705e-07, 6.374221e-07, 6.8881999e-07, 7.4381186e-07, &
    8.0260573e-07, 8.6541864e-07, 9.3247687e-07, 1.0040162e-06, 1.080282e-06, &
    1.1619754e-06, 1.2496277e-06, 1.3436572e-06, 1.4445101e-06, &
    1.5526629e-06, 1.6686237e-06, 1.7929346e-06, 1.9263614e-06, &
    2.0696043e-06, 2.2233774e-06, 2.3884455e-06, 2.5656287e-06, &
    2.7558058e-06, 2.9599187e-06, 3.1791492e-06, 3.4146519e-06, &
    3.6676372e-06, 3.9394056e-06, 4.2313546e-06, 4.544986e-06, 4.8819134e-06, &
    5.2440647e-06, 5.6334087e-06, 6.0520116e-06, 6.5020985e-06, &
    6.9860655e-06, 7.5064932e-06, 8.0661613e-06, 8.6701024e-06, &
    9.3223836e-06, 1.002711e-05, 1.0788756e-05, 1.1612197e-05, 1.2502751e-05, &
    1.3466215e-05, 1.4531221e-05, 1.5710984e-05, 1.7020447e-05, &
    1.8476809e-05, 2.0099909e-05, 2.1912687e-05, 2.4079993e-05, 2.673686e-05, &
    3.0029676e-05, 3.416049e-05, 3.9624245e-05, 4.7318986e-05, 5.8509432e-05, &
    7.3993635e-05, 9.5598711e-05, 0.00012447961, 0.00016290808, &
    0.00021346928, 0.00027992192, 0.00036715931, 0.00048162804, 0.000631805, &
    0.000829155, 0.00108274, 0.00140685, 0.00181885, 0.0023398, 0.00299505, &
    0.0038147, 0.00483445, 0.00609635, 0.00764935, 0.0095501, 0.011864, &
    0.0146655, 0.018038, 0.0220755, 0.0268825, 0.0325735, 0.039273, &
    0.0471145, 0.0562405, 0.0668005, 0.0789485, 0.07731271, 0.07590131, &
    0.07424086, 0.07228743, 0.06998932, 0.06728574, 0.06410509, 0.06036322, &
    0.05596111, 0.05078225, 0.0446896, 0.03752191, 0.02908949, 0.02084739, &
    0.01334443, 0.00708499, 0.00252136, 0, 0 /
    
  data hybi125 /&
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.01505309, 0.03276228, 0.05359622, &
    0.07810627, 0.1069411, 0.1408637, 0.180772, 0.227722, 0.2829562, &
    0.3479364, 0.4243822, 0.5143168, 0.6201202, 0.7235355, 0.8176768, &
    0.8962153, 0.9534761, 0.9851122, 1 /

  data hyam125 /&
    5.9626233e-09, 9.8306813e-09, 1.6190019e-08, 2.65748775e-08, &
    4.30895015e-08, 6.738806e-08, 9.880529e-08, 1.34718675e-07, &
    1.72488915e-07, 2.0998909e-07, 2.46648835e-07, 2.8354663e-07, &
    3.2109086e-07, 3.58751755e-07, 3.96441715e-07, 4.3534552e-07, &
    4.76477275e-07, 5.1982255e-07, 5.6535658e-07, 6.13044575e-07, &
    6.63121045e-07, 7.16315925e-07, 7.73208795e-07, 8.34012185e-07, &
    8.98947755e-07, 9.68246535e-07, 1.0421491e-06, 1.1211287e-06, &
    1.20580155e-06, 1.29664245e-06, 1.39408365e-06, 1.4985865e-06, &
    1.6106433e-06, 1.73077915e-06, 1.859648e-06, 1.99798285e-06, &
    2.14649085e-06, 2.30591145e-06, 2.4770371e-06, 2.66071725e-06, &
    2.85786225e-06, 3.06953395e-06, 3.29690055e-06, 3.54114455e-06, &
    3.8035214e-06, 4.0853801e-06, 4.3881703e-06, 4.7134497e-06, &
    5.06298905e-06, 5.4387367e-06, 5.84271015e-06, 6.27705505e-06, &
    6.744082e-06, 7.24627935e-06, 7.78632725e-06, 8.36813185e-06, &
    8.996243e-06, 9.6747468e-06, 1.0407933e-05, 1.12004765e-05, &
    1.2057474e-05, 1.2984483e-05, 1.3998718e-05, 1.51211025e-05, &
    1.63657155e-05, 1.7748628e-05, 1.9288359e-05, 2.1006298e-05, &
    2.299634e-05, 2.54084265e-05, 2.8383268e-05, 3.2095083e-05, &
    3.68923675e-05, 4.34716155e-05, 5.2914209e-05, 6.62515335e-05, &
    8.4796173e-05, 0.0001100391605, 0.000143693845, 0.00018818868, &
    0.0002466956, 0.000323540615, 0.000424393675, 0.00055671652, 0.00073048, &
    0.0009559475, 0.001244795, 0.00161285, 0.002079325, 0.002667425, &
    0.003404875, 0.004324575, 0.0054654, 0.00687285, 0.008599725, 0.01070705, &
    0.01326475, 0.01635175, 0.02005675, 0.024479, 0.029728, 0.03592325, &
    0.04319375, 0.0516775, 0.0615205, 0.0728745, 0.078130605, 0.07660701, &
    0.075071085, 0.073264145, 0.071138375, 0.06863753, 0.065695415, &
    0.062234155, 0.058162165, 0.05337168, 0.047735925, 0.041105755, &
    0.0333057, 0.02496844, 0.01709591, 0.01021471, 0.00480317499999999, &
    0.00126067999999999, 0 /

 data hybm125 /&
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8.35614355487735e-18, 0.007526545, &
    0.023907685, 0.04317925, 0.065851245, 0.092523685, 0.1239024, 0.16081785, &
    0.204247, 0.2553391, 0.3154463, 0.3861593, 0.4693495, 0.5672185, &
    0.67182785, 0.77060615, 0.85694605, 0.9248457, 0.96929415, 0.9925561 /
      
  ! Do the GEOS case.
  a = a72
  b = b72

  ! Do the WACCM 66 case.
!  a = hyai66 * 10000_f
!  b = hybi66
  
  ! Do the WACCM 125 case.
!  a = hyai125 * 10000_f
!  b = hybi125

!  do i = 1, NZP1
!    write(*,*) a(i), b(i)*p0, a(i) + b(i)*p0
!  end do

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


!  write(*,*) ""
!  write(*,*) "Sedimentation of Dust Particles"

  ! Open the output text file
  open(unit=lun,file="carma_sigmafalltest.txt",status="unknown")
  
  ! Allocate the arrays that we need for the model
  allocate(xc(NZ,NY,NX), dx(NZ,NY,NX), yc(NZ,NY,NX), dy(NZ,NY,NX), &
           zc(NZ,NY,NX), zl(NZP1,NY,NX), p(NZ,NY,NX), pl(NZP1,NY,NX), &
           t(NZ,NY,NX))
  allocate(mmr(NZ,NY,NX,NELEM,NBIN))
  allocate(lat(NY,NX), lon(NY,NX))  

  ! Define the particle-grid extent of the CARMA test
!  write(*,*) "  CARMA_Create(carma, ", NBIN,    ", ", NELEM, ", ", NGROUP, &
!                                 ", ", NSOLUTE, ", ", NGAS, ", rc, 6) ..."
  call CARMA_Create(carma, NBIN, NELEM, NGROUP, NSOLUTE, NGAS, NWAVE, rc, LUNOPRT=6)
  if (rc /=0) write(*, *) "    *** FAILED ***, rc=", rc
	carma_ptr => carma

  ! Define the group
!  write(*,*) "  Add Group(s) ..."
  rho = 2.65_f
!  rmrat = (100._f**3)**(1._f/(NBIN*1._f))
!  rmin = 1.e-5_f * ((1._f+rmrat)/2._f)**(1._f/3._f)
  rmrat = 2.0
  rmin = 7.5e-4_f
  call CARMAGROUP_Create(carma, 1, 'dust', rmin, rmrat, I_SPHERE, 1._f, .FALSE., rc)
  if (rc /=0) stop "    *** FAILED ***"
  
  ! Define the element
!  write(*,*) "  Add Element(s) ..."
  call CARMAELEMENT_Create(carma, 1, 1, "dust", rho, I_INVOLATILE, I_DUST, rc)
  if (rc /=0) stop "    *** FAILED ***"
  
!  write(*,*) "  CARMA_AddCoagulation(carma, 1, 1, 1, I_COLLEC_DATA, rc) ..."
!  call CARMA_AddCoagulation(carma, 1, 1, 1, I_COLLEC_DATA, rc) 
!  if (rc /=0) write(*, *) "    *** FAILED ***, rc=", rc
  
! Setup the CARMA processes to exercise
!  write(*,*) "  CARMA_Initialize(carma, rc, do_vtran=.TRUE., "// &
!               "vf_const=", vf_const,", do_explised=",do_explised,") ..."
  call CARMA_Initialize(carma, rc, do_vtran=.TRUE., vf_const=vf_const, do_explised=do_explised)
  if (rc /=0) write(*, *) "    *** FAILED ***, rc=", rc
  
! Print the Group Information
!  write(*,*)  ""
!  call dumpGroup(carma, rc)
!  if (rc /=0) write(*, *) "    *** FAILED ***, rc=", rc
  
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

  ! Layer edges
  do i = 1, NZP1
!   pl(i,:,:) = a(i)+b(i)*p0
!   zl(i,:,:) = pl(i,:,:)/p0
!   zl(i,:,:) = a(i) + b(i)*p0
   pl(i,:,:) = a72(i)+b72(i)*ps
   zl(i,:,:) = pl(i,:,:)/ps
  end do

  
  ! Vertical center
!  t = 270._f
  do i = 1, NZ
!    p(i,:,:) = exp((log(pl(i,:,:)) + log(pl(i+1,:,:)) ) / 2._f)
!    zc(i,:,:) = p(i,:,:) / p0
!    zc(i,:,:) = hyam66(i) * 10000_f + hybm66(i) * p0
!    zc(i,:,:) = hyam125(i) * 10000_f + hybm125(i) * p0
    t(i,:,:) = t72(i)
    p(i,:,:) = (pl(i,:,:) + pl(i+1,:,:) ) / 2._f
    zc(i,:,:) = p(i,:,:) / p0
    rhoa(i) = p(i,NY,NX)/287._f/t(i,NY,NX)
    dz(i) = ze(i)-ze(i+1)
  end do

!  zm(NZ) = (p0 - p(NZ,NX,NY)) /  (rhoa(NZ)*9.816_f)
  zm(NZ) = ze(NZP1)+dz(NZ)/2._f

  do i = NZ-1, 1, -1
    zm(i) = zm(i+1) + (p(i+1,NY,NX) - p(i,NY,NX)) /  (rhoa(i)*9.816_f)
  end do

!  write(*,'(a6, 3a12)') "level", "zc", "p", "t"
!  write(*,'(a6, 3a12)') "", "(m)", "(Pa)", "(K)"
!  do i = 1, NZ
!    write(*,'(i6,3f12.3)') i, zc(i,NY,NX), p(i,NY,NX), t(i,NY,NX)
!  end do


!  write(*,*) ""
!  write(*,'(a6, 2a12)') "level", "zl", "pl"
!  write(*,'(a6, 2a12)') "", "(m)", "(Pa)"
!  do i = 1, NZP1
!    write(*,'(i6,2f12.3)') i, zl(i,NY,NX), pl(i,NY,NX)
!  end do

  			
  ! Put a blob in the model first bin at 8 km
  mmr(:,:,:,:,:) = 0._f
  do i = 1, NZ
   mmr(i,:,:,1,1) = 1e-10_f * exp( - ( ( zm(i) - 8.e3_f)/3.e3_f) ** 2) / &
                     ( p(i,:,:) / 287._f / t(i,:,:))
  end do

!  write(*,*)  ""
!  write(*, '(a6, 4a12)') "level", "mmr(i,NY,NX,1)", "mmr(i,NY,NX,2)"
!  do i = 1, NZ
!	  write(*, '(i6, 4g12.3)') i, mmr(i,NY,NX,1,1), mmr(i,NY,NX,1,2)
!  end do
  
  ! Write output for the falltest
  write(lun,*) NZ
  do i = 1, NZ
   write(lun,'(i3,2f10.1)') &
    i, zm(i), dz(i)
  end do
  
  write(lun,*) 0
  do i = 1, NZ
   write(lun,'(i3,e10.3,e10.3)') &
    i, real(mmr(i,NY,NX,1,1)), real(mmr(i,NY,NX,1,1)*rhoa(i))
  end do

		
  ! Iterate the model over a few time steps.
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
                              I_HYBRID, I_CART, lat(iy,ix), lon(iy,ix), &
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

    ! Write output for the falltest
    write(lun,'(f12.0)') istep*dtime
    do i = 1, NZ
      write(lun,'(i3,e10.3,e10.3)') &
        i, real(mmr(i,NY,NX,1,1)), real(mmr(i,NY,NX,1,1)*rhoa(i))
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

  if (rc /=0) write(*, *) "    *** FAILED ***, rc=", rc

!  write(*,*)  ""
!  write(*,*) "  CARMA_Destroy() ..."
  call CARMA_Destroy(carma, rc)
  if (rc /=0) write(*, *) "    *** FAILED ***, rc=", rc  
end subroutine

