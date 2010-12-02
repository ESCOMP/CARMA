!! This module defines types used in the CARMA module. The types need to be defined here
!! to avoid circular references between different modules (e.g. carma_mod and
!! carmastate_mod).
!!
!! @version July-2009 
!! @author  Chuck Bardeen 
module carma_types_mod
  use carma_precision_mod
  use carma_constants_mod

  !! The CARMAELEMENT data type represents one of the components of a cloud or aerosol particle.
  !!
  !! The procedure for adding a variable to the CARMAELEMENT data type is:
  !!  - Add the variable as a scalar or allocatable in the type definition.
  !!  - If the new variable is dynamic,
  !!    - Allocate the variable in the appropriate create or initialization routine.
  !!    - Deallocate the variable in the approprate finalize and destroy routines.
  !!  - Add an alias for the variable to carma_globaer.h and associate it with the variable
  !!     in this typedef.  
  !!
  !! NOTE: While the carmaelement_type is public, routines outside of the CARMA module should not look
  !! at or manuipulate fields of this structure directly. There should be CARMAELEMENT_XXX methods
  !! to do anything that is needed with this structure, and use of these methods will allow
  !! the CARMAELEMENT data type structure to evolve without impacting code in the parent model.
  !! The contents of the structure had to be made public, since the CARMA microphysics
  !! routines are implemented in separate files outside of this model; however, logically
  !! they are part of the model and are the only routines outside of this module that should
  !! access fields of this structure directly.
  type, public :: carmaelement_type
  
    !   name          Name of the element
    !   shortname     Short name of the element
    !   rho           Mass density of particle element [g/cm^3]
    !   igroup        Group to which the element belongs
    !   itype         Particle type specification
    !   icomposition  Particle compound specification
    !   isolute       Index of solute for the particle element
    !
    character(len=CARMA_NAME_LEN)               :: name
    character(len=CARMA_SHORT_NAME_LEN)         :: shortname
    real(kind=f), allocatable, dimension(:)     :: rho          ! (NBIN)
    integer                                     :: igroup
    integer                                     :: itype
    integer                                     :: icomposition
    integer                                     :: isolute
  end type carmaelement_type


  !! The CARMAGAS data type represents a gas.
  !!
  !! The procedure for adding a variable to the CARMAGAS data type is:
  !!  - Add the variable as a scalar or allocatable in the type definition.
  !!  - If the new variable is dynamic,
  !!    - Allocate the variable in the appropriate create or initialization routine.
  !!    - Deallocate the variable in the approprate finalize and destroy routines.
  !!  - Add an alias for the variable to carma_globaer.h and associate it with the variable
  !!     in this typedef.  
  !!
  !! NOTE: While the carmagas_type is public, routines outside of the CARMA module should not look
  !! at or manuipulate fields of this structure directly. There should be CARMAGAS_XXX methods
  !! to do anything that is needed with this structure, and use of these methods will allow
  !! the CARMAGAS data type structure to evolve without impacting code in the parent model.
  !! The contents of the structure had to be made public, since the CARMA microphysics
  !! routines are implemented in separate files outside of this model; however, logically
  !! they are part of the model and are the only routines outside of this module that should
  !! access fields of this structure directly.
  type, public :: carmagas_type
  
    !   name        Name of the gas
    !   shortname   Short name of the gas
    !   wtmol       Molecular weight for the gas [g/mol]
    !   ivaprtn     vapor pressure routine for the gas
    !
    character(len=CARMA_NAME_LEN)               :: name
    character(len=CARMA_SHORT_NAME_LEN)         :: shortname
    real(kind=f)                                :: wtmol
    integer                                     :: ivaprtn
    integer                                     :: icomposition
  end type carmagas_type


  !! The CARMAGROUP data type represents a cloud or aerosol partcile.
  !!
  !! The procedure for adding a variable to the CARMAGROUP data type is:
  !!  - Add the variable as a scalar or allocatable in the type definition.
  !!  - If the new variable is dynamic,
  !!    - Allocate the variable in the appropriate create or initialization routine.
  !!    - Deallocate the variable in the approprate finalize and destroy routines.
  !!  - Add an alias for the variable to carma_globaer.h and associate it with the variable
  !!     in this typedef.  
  !!
  !! NOTE: While the carmagroup_type is public, routines outside of the CARMA module should not look
  !! at or manuipulate fields of this structure directly. There should be CARMAGROUP_XXX methods
  !! to do anything that is needed with this structure, and use of these methods will allow
  !! the CARMAGROUP data type structure to evolve without impacting code in the parent model.
  !! The contents of the structure had to be made public, since the CARMA microphysics
  !! routines are implemented in separate files outside of this model; however, logically
  !! they are part of the model and are the only routines outside of this module that should
  !! access fields of this structure directly.
  type, public :: carmagroup_type
  
    !   name        Name of the particle
    !   shortname   Short name of the particle
    !   cnsttype    constituent type [I_CNSTTYPE_PROGNOSTIC | I_CNSTTYPE_DIAGNOSTIC]
    !   maxbin      the last prognostic bin in the group
    !   nelem       Number of elements in group           
    !   ncore       Number of core elements (itype = 2) in group           
    !   ishape      Describes particle shape for group
    !   ienconc     Particle number conc. element for group
    !   imomelem    Scondary moment element for group
    !   icorelem    Core elements (itype = 2) in group           
    !   solfac      Solubility factor for wet deposition
    !   is_ice      If .true. then ice particle 
    !   is_cloud    If .true. then cloud particle
    !   do_mie      If .true. then do mie calculations 
    !   do_wetdep   If .true. then do wet deposition 
    !   grp_do_drydep If .true. then do dry deposition 
    !   grp_do_vtran If .true. then do sedimentation 
    !   scavcoef    Scavenging coefficient for wet deopistion (1/mm)
    !   if_sec_mom  If .true. then core second moment (itype = 3) used    {setupgrow}
    !   irhswell    Indicates method for swelling particles from RH
    !   irhswcomp   Indicates composition for swelling particles from RH
    !   rmin        Radius of particle in first bin [cm]
    !   rmassmin    Mass of particle in first bin [g]
    !   rmrat       Ratio of masses of particles in consecutive bins
    !   eshape      Ratio of particle length / diameter 
    !   r           Radius bins [cm]
    !   rmass       Mass bins [g]
    !   vol         Particle volume [cm^3]
    !   dr          Width of bins in radius space [cm]
    !   dm          Width of bins in mass space [g]
    !   rmassup     Upper bin boundary mass [g]
    !   rup         Upper bin boundary radius [cm]
    !   rlow        Lower bin boundary radius [cm]
    !   refidx      refractive index
    !   qext        extinction efficiency
    !   ssa         single scattering albedo
    !   asym        asymmetry factor
    !   ifallrtn    routine to use to calculate fall velocity  [I_FALLRTN_...]

    !
    character(len=CARMA_NAME_LEN)               :: name
    character(len=CARMA_SHORT_NAME_LEN)         :: shortname
    integer                                     :: cnsttype
    integer                                     :: maxbin
    integer                                     :: nelem
    integer                                     :: ncore
    integer                                     :: ishape
    integer                                     :: ienconc
    integer                                     :: imomelem
    real(kind=f)                                :: solfac
    real(kind=f)                                :: scavcoef
    logical                                     :: if_sec_mom
    logical                                     :: is_ice
    logical                                     :: is_cloud
    logical                                     :: do_mie
    logical                                     :: do_wetdep
    logical                                     :: grp_do_drydep
    logical                                     :: grp_do_vtran
    integer                                     :: irhswell
    integer                                     :: irhswcomp
    integer                                     :: ifallrtn
    real(kind=f)                                :: rmin
    real(kind=f)                                :: rmassmin
    real(kind=f)                                :: rmrat
    real(kind=f)                                :: eshape
    real(kind=f), allocatable, dimension(:)     :: r          ! (NBIN)
    real(kind=f), allocatable, dimension(:)     :: rmass      ! (NBIN)
    real(kind=f), allocatable, dimension(:)     :: vol        ! (NBIN)
    real(kind=f), allocatable, dimension(:)     :: dr         ! (NBIN)
    real(kind=f), allocatable, dimension(:)     :: dm         ! (NBIN)
    real(kind=f), allocatable, dimension(:)     :: rmassup    ! (NBIN)
    real(kind=f), allocatable, dimension(:)     :: rup        ! (NBIN)
    real(kind=f), allocatable, dimension(:)     :: rlow       ! (NBIN)
    complex(kind=f), allocatable, dimension(:)  :: refidx     ! (NWAVE)
    real(kind=f), allocatable, dimension(:,:)   :: qext       ! (NWAVE,NBIN)
    real(kind=f), allocatable, dimension(:,:)   :: ssa        ! (NWAVE,NBIN)
    real(kind=f), allocatable, dimension(:,:)   :: asym       ! (NWAVE,NBIN)
    integer, allocatable, dimension(:)          :: icorelem   ! (NELEM)
  end type carmagroup_type
  
  
  !! The CARMASOLUTE data type represents a gas.
  !!
  !! The procedure for adding a variable to the CARMASOLUTE data type is:
  !!  - Add the variable as a scalar or allocatable in the type definition.
  !!  - If the new variable is dynamic,
  !!    - Allocate the variable in the appropriate create or initialization routine.
  !!    - Deallocate the variable in the approprate finalize and destroy routines.
  !!  - Add an alias for the variable to carma_globaer.h and associate it with the variable
  !!     in this typedef.  
  !!
  !! NOTE: While the carmagas_type is public, routines outside of the CARMA module should not look
  !! at or manuipulate fields of this structure directly. There should be CARMASOLUTE_XXX methods
  !! to do anything that is needed with this structure, and use of these methods will allow
  !! the CARMASOLUTE data type structure to evolve without impacting code in the parent model.
  !! The contents of the structure had to be made public, since the CARMA microphysics
  !! routines are implemented in separate files outside of this model; however, logically
  !! they are part of the model and are the only routines outside of this module that should
  !! access fields of this structure directly.
  type, public :: carmasolute_type
  
    !   name        Name of the solute
    !   shortname   Short name of the solute
    !   ions        Number of ions solute dissociates into
    !   wtmol       Molecular weight of solute
    !   rho         Mass density of solute
    !
    character(len=CARMA_NAME_LEN)               :: name
    character(len=CARMA_SHORT_NAME_LEN)         :: shortname
    integer                                     :: ions
    real(kind=f)                                :: wtmol
    real(kind=f)                                :: rho
  end type carmasolute_type


  !! The CARMA data type replaces the common blocks that were used in the F77 version of
  !! CARMA. This allows the code to be written to allow for multiple threads to call CARMA
  !! routines simulataneously. This thread safety is necessary for to run CARMA under OPEN/MP.
  !!
  !! The procedure for adding a variable to the CARMA data type is:
  !!  - Add the variable as a scalar or allocatable in the type definition.
  !!  - If the new variable is dynamic,
  !!    - Allocate the variable in the appropriate create or initialization routine.
  !!    - Deallocate the variable in the approprate finalize and destroy routines.
  !!  - Add an alias for the variable to carma_globaer.h and associate it with the variable
  !!     in this typedef.  
  !!
  !! NOTE: While the carmatype is public, routines outside of the CARMA module should not look
  !! at or manuipulate fields of this structure directly. There should be CARMA_XXX methods
  !! to do anything that is needed with this structure, and use of these methods will allow
  !! the CARMA data type structure to evolve without impacting code in the parent model.
  !! The contents of the structure had to be made public, since the CARMA microphysics
  !! rountines are implemented in separate files outside of this model; however, logically
  !! they are part of the model and are the only routines outside of this module that should
  !! access fields of this structure directly.
  type, public :: carma_type
  
    ! Model Dimensions
    !
    !  NGROUP   number of particle groups
    !  NELEM    number of particle components (elements)
    !  NBIN     number of size bins per element
    !  NGAS     number of gases (may be 0)
    !  NSOLUTE  number of solutes (may be 0)
    !  NWAVE    number of wavelength bands (may be 0)
    !
    integer :: NGROUP
    integer :: NELEM
    integer :: NBIN
    integer :: NGAS
    integer :: NSOLUTE
    integer :: NWAVE

    ! Output logical unit numbers
    !
    ! NOTE: CARMA will not directly access files or keep track of file names. It is the
    ! parent model's responsibility to provide the logical unit number to be used for
    ! model output.
    !
    integer :: LUNOPRT   ! output print file

    ! Model startup control variables
    !
    !   do_print  .t. if print output is desired
    !
    logical            :: do_print
    
    
    ! Configuration Objects
    !
    ! These are all other objects that are parts of the CARMA model. This is
    ! an attempt to break up the large common block that has historically been
    ! the structure of CARMA so the code is easier to understand and to
    ! maintain.
    !
    !   element     Particle component    
    !   gas         Gas   
    !   group       Particle    
    !   solute      Element solute    
    !
    ! NOTE: In the future, it may make sense to create objects that represent
    ! the CARMA processes. This would encapsulate all the variables related to
    ! a particular process into one structure. Candidate processes include:
    ! transport, growth, nucleation, coagulation, ...
    !
    type(carmaelement_type), allocatable, dimension(:)    :: element    ! (NELEM)
    type(carmagas_type),     allocatable, dimension(:)    :: gas        ! (NGAS)
    type(carmagroup_type),   allocatable, dimension(:)    :: group      ! (NGROUP)
    type(carmasolute_type),  allocatable, dimension(:)    :: solute     ! (NSOLUTE)



    ! Model option & control variables
    !
    !   conmax      Minumum relative concentration to consider in varstep   {prestep}   
    !   icoag       Coagulation mapping array                           {setupcoag}
    !   icoagelem   Coagulation element mapping array                   {setupcoag}
    !   icoagelem_cm Coagulation element mapping array for second mom   {setupcoag}
    !   ifall       Fall velocity options                               {setupvfall}
    !   icoagop     Coagulation kernel options                          {setupckern}
    !   icollec     Gravitational collection options                      {setupckern}
    !   itbnd_pc    Top boundary condition flag for particles             {init}
    !   ibbnd_pc    Bottom boundary condition flag for particles          {init}
    !   do_vdiff    If .true. then do Brownian diffusion                  {init}
    !   do_coag     If .true. then do coagulation                         {init}
    !   do_detrain  If .true. then do detrainment                         {init}
    !   do_drydep   If .true. then do dry deposition                      {init}
    !   do_fixedinitIf .true. then do initialize from reference atm       {init}
    !   do_grow     If .true. then do condensational growth and evap.     {init}
    !   do_incloud  If .true. then do incloud growth and coagulation      {init}
    !   do_explised If .true. then do sedimentation with substepping      {init}
    !   do_print_init If .true. then do print initializtion info          {init}
    !   do_step     if .true. then varstepping succeeded                  {init}
    !   do_substep  if .true. then use substepping                        {init}
    !   do_thermo   if .true. then do solve thermodynamic equation        {init}
    !   do_vdiff    If .true. then do Brownian diffusion                  {init}
    !   do_vtran    If .true. then do vertical transport                  {init}
    !   do_cnst_rlh If .true. then uses constants for rlhe and rlhm       {setupgrow}
    !   igrowgas    Gas that condenses into a particle element            {setupgrow}
    !   inucgas     Gas that nucleates a particle group                   {setupnuc}
    !   if_nuc      Nucleation conditional array                          {setupaer}
    !   inucproc    Nucleation conditional array                          {setupaer}
    !   nnuc2elem   Number of elements that nucleate to element           {setupnuc}
    !   inuc2elem   Nucleation transfers particles into element inuc2elem {setupnuc}
    !   ievp2elem   Total evap. transfers particles into group ievp2elem  {setupnuc}
    !   ievp2bin    Total evap. transfers particles into bin ievp2bin     {setupnuc}
    !   inuc2bin    Nucleation transfers particles into bin inuc2bin      {setupnuc}
    !   maxsubsteps Maximum number of time substeps allowed
    !   minsubsteps Maximum number of time substeps allowed
    !   maxretries  Maximum number of substepping retries allowed
    !
    logical                                       :: do_vdiff
    logical                                       :: do_drydep
    logical                                       :: do_coag
    logical                                       :: do_detrain
    logical                                       :: do_fixedinit
    logical                                       :: do_grow
    logical                                       :: do_incloud
    logical                                       :: do_vtran
    logical                                       :: do_explised
    logical                                       :: do_print_init
    logical                                       :: do_step
    logical                                       :: do_substep
    logical                                       :: do_thermo
    logical                                       :: do_cnst_rlh
    logical, allocatable, dimension(:,:)          :: if_nuc       !(NELEM,NELEM)
    real(kind=f)                                  :: conmax
    integer                                       :: maxsubsteps 
    integer                                       :: minsubsteps 
    integer                                       :: maxretries 
    integer                                       :: ifall
    integer                                       :: icoagop
    integer                                       :: icollec
    integer                                       :: itbnd_pc
    integer                                       :: ibbnd_pc
    integer, allocatable, dimension(:)            :: inucgas      ! NGROUP
    integer, allocatable, dimension(:)            :: igrowgas     ! NELEM
    integer, allocatable, dimension(:)            :: nnuc2elem    ! NELEM
    integer, allocatable, dimension(:)            :: ievp2elem    ! NELEM
    integer, allocatable, dimension(:)            :: nnucelem     ! NELEM
    integer, allocatable, dimension(:,:)          :: icoag        ! (NGROUP,NGROUP)
    integer, allocatable, dimension(:,:)          :: inucproc     ! (NELEM,NELEM)
    integer, allocatable, dimension(:,:)          :: inuc2elem    ! (NELEM,NELEM)
    integer, allocatable, dimension(:,:)          :: icoagelem    ! (NELEM,NGROUP)
    integer, allocatable, dimension(:,:)          :: icoagelem_cm ! (NELEM,NGROUP)
    integer, allocatable, dimension(:,:)          :: inucelem     ! (NELEM,NELEM*NGROUP)
    integer, allocatable, dimension(:,:,:)        :: inuc2bin     ! (NBIN,NGROUP,NGROUP)
    integer, allocatable, dimension(:,:,:)        :: ievp2bin     ! (NBIN,NGROUP,NGROUP)
    integer, allocatable, dimension(:,:,:)        :: nnucbin      ! (NGROUP,NBIN,NGROUP)
    integer, allocatable, dimension(:,:,:,:)      :: inucbin      ! (NBIN*NGROUP,NGROUP,NBIN,NGROUP)
  

    ! Particle bin structure
    !
    !   diffmass  Difference between <rmass> values
    !
    real(kind=f), allocatable, dimension(:,:,:,:)   :: diffmass   ! (NBIN,NGROUP,NBIN,NGROUP)

    !  Coagulation kernels and bin pair mapping
    !
    !   ck0           Constant coagulation kernel           {setupaer}
    !   grav_e_coll0  Constant value for collection effic.  {setupaer}
    !   volx          Coagulation subdivision variable      {setupcoag}
    !   ilow          Bin pairs for coagulation production  {setupcoag}
    !   jlow          Bin pairs for coagulation production  {setupcoag}
    !   iup           Bin pairs for coagulation production  {setupcoag}
    !   jup           Bin pairs for coagulation production  {setupcoag}
    !   npairl        Bin pair indices                      {setupcoag}
    !   npairu        Bin pair indices                      {setupcoag}
    !   kbin          lower bin for coagulation             {setupcoag}
    !
    real(kind=f)                                        :: ck0
    real(kind=f)                                        :: grav_e_coll0
    real(kind=f), allocatable, dimension(:,:,:,:,:)     :: volx    ! (NGROUP,NGROUP,NGROUP,NBIN,NBIN)
    integer, allocatable, dimension(:,:,:)              :: ilow    ! (NGROUP,NBIN,NBIN*NBIN)
    integer, allocatable, dimension(:,:,:)              :: jlow    ! (NGROUP,NBIN,NBIN*NBIN)
    integer, allocatable, dimension(:,:,:)              :: iup     ! (NGROUP,NBIN,NBIN*NBIN)
    integer, allocatable, dimension(:,:,:)              :: jup     ! (NGROUP,NBIN,NBIN*NBIN)
    integer, allocatable, dimension(:,:)                :: npairl  ! (NGROUP,NBIN)
    integer, allocatable, dimension(:,:)                :: npairu  ! (NGROUP,NBIN)
    integer, allocatable, dimension(:,:,:,:,:)          :: kbin    ! (NGROUP,NGROUP,NGROUP,NBIN,NBIN)

    !  Coagulation group pair mapping
    !
    !   iglow      Group pairs for coagulation production  {setupcoag}
    !   jglow      Group pairs for coagulation production  {setupcoag}
    !   igup       Group pairs for coagulation production  {setupcoag}
    !   jgup       Group pairs for coagulation production  {setupcoag}
    !
    integer, allocatable, dimension(:,:,:) :: iglow  ! (NGROUP,NBIN,NBIN*NBIN)
    integer, allocatable, dimension(:,:,:) :: jglow  ! (NGROUP,NBIN,NBIN*NBIN)
    integer, allocatable, dimension(:,:,:) :: igup   ! (NGROUP,NBIN,NBIN*NBIN)
    integer, allocatable, dimension(:,:,:) :: jgup   ! (NGROUP,NBIN,NBIN*NBIN)

    !  Particle fall velocities
    !
    !   vf_const  Constant vertical fall velocity when ifall=0          {setupaer}
    !
    real(kind=f)                                        :: vf_const


    ! Condensational growth parameters
    !
    ! NOTE: Some of these variables are used for storing intermediate values in
    ! the calculations. They may no longer be necessary, when the code is
    ! implemented as F90 and values as passed as parameters between subroutines.
    !
    !   rlh_nuc   Latent heat released by nucleation [cm^2/s^2]       {setupaer}
    !   pratt     Terms in PPM advection scheme for condensation      {setupgkern}
    !   prat
    !   pden1
    !   palr
    real(kind=f), allocatable, dimension(:,:)        :: rlh_nuc    ! (NELEM,NELEM)
    real(kind=f), allocatable, dimension(:,:,:)      :: pratt      ! (3,NBIN,NGROUP)
    real(kind=f), allocatable, dimension(:,:,:)      :: prat       ! (4,NBIN,NGROUP)
    real(kind=f), allocatable, dimension(:,:)        :: pden1      ! (NBIN,NGROUP)
    real(kind=f), allocatable, dimension(:,:)        :: palr       ! (4,NGROUP)
    
    ! Optical Properties
    !   wave      Bin-center wavelengths [micron]
    !
    real(kind=f), allocatable, dimension(:)          :: wave       ! (NWAVE)
  end type carma_type  
 
 
  !! The cstate data type replaces portions of the common blocks that were used
  !! in the F77 version of CARMA. This allows the code to be written to allow for
  !! multiple threads to call CARMA routines simulataneously. This thread safety is
  !! necessary for to run CARMA under OPEN/MP.
  !!
  !! The procedure for adding a variable to the cstate data type is:
  !!  - Add the variable as a scalar or allocatable in the type definition.
  !!  - If the new variable is dynamic,
  !!    - Allocate the variable in the create routine.
  !!    - Deallocate the variable in the destroy routines.
  !!  - Add an alias for the variable to cstate.h and associate it with the
  !!      variable in this typedef.  
  !!
  !! NOTE: While the carmastate_type is public, routines outside of the CARMA module
  !! should not look at or manuipulate fields of this structure directly. There should
  !! be CARMASTATE_XXX methods to do anything that is needed with this structure, and
  !! use of these methods will allow the cstate data type structure to evolve without
  !! impacting code in the parent model. The contents of the structure had to be made
  !! public, since the CARMA microphysics rountines are implemented in separate files
  !! outside of this model; however, logically they are part of the model and are the
  !! only routines outside of this module that should access fields of this structure
  !! directly.
  type, public :: carmastate_type
  
    ! Parent CARMA object
    type(carma_type), pointer                   :: carma

    ! Model Dimensions
    !
    !  NZ       number of grid points in the column
    !  NZP1     NZ+1
    !  NGROUP   number of particle groups
    !  NELEM    number of particle components (elements)
    !  NBIN     number of size bins per element
    !  NGAS     number of gases (may be 0)
    !
    integer :: NZ
    integer :: NZP1
    
    ! Model option & control variables
    !
    !   time        Simulation time at end of current timestep [s]
    !   dtime       Substep Timestep size [s]
    !   dtime_orig  Original Timestep size [s]
    !   nretries    Number of substepping retries attempted
    real(kind=f)                                  :: time
    real(kind=f)                                  :: dtime
    real(kind=f)                                  :: dtime_orig
    real(kind=f)                                  :: nretries

    !   max_nretry  Maximum number of retries in a step
    !   nstep       Total number of steps taken
    !   nsubstep    Total number of substeps taken
    !   nretry      Total number of retries taken
    integer                                       :: max_nsubstep
    real(kind=f)                                  :: max_nretry
    real(kind=f)                                  :: nstep
    integer                                       :: nsubstep
    real(kind=f)                                  :: nretry
    
    real(kind=f), allocatable, dimension(:)       :: zsubsteps   ! (NZ)


    ! Model Grid
    !
    !  igridv     flag to specify desired vertical grid coord system    {initatm}
    !  igridh     flag to specify desired horizontal grid coord system  {initatm}
    !  xmet       Horizontal ds/dx (ds is metric distance)              {initatm}
    !  ymet       Horizontal ds/dy (ds is metric distance)              {initatm}
    !  zmet       Vertical ds/dz (ds is metric distance)                {initatm}
    !  zmetl      Vertical ds/dz at edges (ds is metric distance)       {initatm}
    !  xc         Horizontal position at center of box                  {initatm}
    !  yc         Horizontal position at center of box                  {initatm}
    !  zc         Altitude at layer mid-point                           {initatm}
    !  dx         Horizontal grid spacing                               {initatm}
    !  dy         Horizontal grid spacing                               {initatm}
    !  dz         Thickness of vertical layers                          {initatm}
    !  zl         Altitude at top of layer                              {initatm}
    !  lon        Longitude [deg] at xc, yc                             {initatm}
    !  lat        Latitude [deg] at xc, yc                              {initatm}
    !
    integer :: igridv
    integer :: igridh
    real(kind=f), allocatable, dimension(:)     :: xmet   ! (NZ)
    real(kind=f), allocatable, dimension(:)     :: ymet   ! (NZ)
    real(kind=f), allocatable, dimension(:)     :: zmet   ! (NZ)
    real(kind=f), allocatable, dimension(:)     :: zmetl  ! (NZP1)
    real(kind=f), allocatable, dimension(:)     :: xc     ! (NZ)
    real(kind=f), allocatable, dimension(:)     :: yc     ! (NZ)
    real(kind=f), allocatable, dimension(:)     :: zc     ! (NZ)
    real(kind=f), allocatable, dimension(:)     :: dx     ! (NZ)
    real(kind=f), allocatable, dimension(:)     :: dy     ! (NZ)
    real(kind=f), allocatable, dimension(:)     :: dz     ! (NZ)
    real(kind=f), allocatable, dimension(:)     :: zl     ! (NZP1)
    real(kind=f)                                :: lon
    real(kind=f)                                :: lat

    ! Particle bin structure
    !
    !   rhop      Mass density of particle groups [g/cm^3]
    !   r_wet     Wet particle radius from RH swelling [cm]             {setupvfall}
    !   rhop_wet  Wet Mass density of particle groups [g/cm^3]
    !
    real(kind=f), allocatable, dimension(:,:,:) :: rhop       ! (NZ,NBIN,NGROUP)
    real(kind=f), allocatable, dimension(:,:,:) :: rhop_wet   ! (NZ,NBIN,NGROUP)
    real(kind=f), allocatable, dimension(:,:,:) :: r_wet      ! (NZ,NBIN,NGROUP)

    ! Primary model state variables
    !
    !  pc          Particle concentration [/x_units/y_units/z_units]  {initaer}
    !  pcd         Detrained particle concentration [/x_units/y_units/z_units]  {initaer}
    !  pc_surf     Particles on surface [/cm2]                        {initaer}
    !  gc          Gas concentration [g/x_units/y_units/z_units]      {initgas}
    !  cldfrc      Cloud fraction [fraction]
    !  rhcrit      Relative humidity for onset of liquid clouds [fraction]
    !
    real(kind=f), allocatable, dimension(:,:,:) :: pc         ! (NZ,NBIN,NELEM)
    real(kind=f), allocatable, dimension(:,:,:) :: pcd        ! (NZ,NBIN,NELEM)
    real(kind=f), allocatable, dimension(:,:)   :: pc_surf    ! (NBIN,NELEM)
    real(kind=f), allocatable, dimension(:,:)   :: gc         ! (NZ,NGAS)
    real(kind=f), allocatable, dimension(:)     :: cldfrc     ! (NZ)
    real(kind=f), allocatable, dimension(:)     :: rhcrit     ! (NZ)

    ! Secondary model variables
    !
    ! NOTE: Some of these variables are used for storing intermediate values in
    ! the calculations. They may no longer be necessary, when the code is
    ! implemented as F90 and values as passed as parameters between subroutines.
    !
    !   pcl         Particle concentration at beginning of time-step
    !   pconmax     Maximum particle concentration for each grid point
    !   gcl         Gas concentration at beginning of time-step
    !   d_gc        Change in gas concentration due to transport
    !   d_t        Change in temperature due to transport
    !   dpc_sed     Change in particle concentration due to sedimentation
    !   coaglg      Total particle loss rate due to coagulation for group
    !   coagpe      Particle production due to coagulation 
    !   rnuclg      Total particle loss rate due to nucleation for group
    !   rnucpe      Particle production due to nucleation 
    !   pc_nucl     Particles produced due to nucleation (for the whole step, not just the substep)
    !   growlg      Total particle loss rate due to growth for group
    !   growle      Partial particle loss rate due to growth for element 
    !   growpe      Particle production due to growth 
    !   evaplg      Total particle loss rate due to evaporation for group
    !   evapls      Partial particle loss rate due to evaporation for element
    !   evappe      Particle production due to evaporation
    !   coreavg     Average total core mass in bin
    !   coresig     logarithm^2 of std dev of core distribution
    !   evdrop      Particle production of droplet number
    !   evcore      Particle production of core elements
    !   gasprod     Gas production term
    !   rlheat      Latent heating rate [deg_K/s]   
    !   ftoppart    Downward particle flux across top boundary of model
    !   fbotpart    Upward flux particle across bottom boundary of model
    !   pc_topbnd   Particle concentration assumed just above the top boundary
    !   pc_botbnd   Particle concentration assumed just below the bottom boundary
    !   cmf         Core mass fraction in a droplet 
    !   totevap     .true. if droplets are totally evaporating to CN
    !   too_small   .true. if cores are smaller than smallest CN
    !   too_big     .true. if cores are larger than largest CN
    !   nuc_small   .true. if cores are smaller than smallest nucleated CN
    !
    real(kind=f), allocatable, dimension(:,:,:)   :: pcl        ! (NZ,NBIN,NELEM
    real(kind=f), allocatable, dimension(:,:)     :: gcl        ! (NZ,NGAS)
    real(kind=f), allocatable, dimension(:,:)     :: d_gc       ! (NZ,NGAS)
    real(kind=f), allocatable, dimension(:)       :: d_t        ! (NZ)
    real(kind=f), allocatable, dimension(:,:)     :: dpc_sed    ! (NBIN,NELEM)
    real(kind=f), allocatable, dimension(:,:)     :: pconmax    ! (NZ,NGROUP)
    real(kind=f), allocatable, dimension(:,:,:)   :: coaglg     ! (NZ,NBIN,NGROUP)
    real(kind=f), allocatable, dimension(:,:,:)   :: coagpe     ! (NZ,NBIN,NELEM)
    real(kind=f), allocatable, dimension(:,:,:)   :: rnuclg     ! (NBIN,NGROUP,NGROUP)
    real(kind=f), allocatable, dimension(:,:)     :: rnucpe     ! (NBIN,NELEM)
    real(kind=f), allocatable, dimension(:,:,:)   :: pc_nucl    ! (NZ,NBIN,NELEM)
    real(kind=f), allocatable, dimension(:,:)     :: growpe     ! (NBIN,NELEM)
    real(kind=f), allocatable, dimension(:,:)     :: evappe     ! (NBIN,NELEM)
    real(kind=f)                                  :: coreavg
    real(kind=f)                                  :: coresig
    real(kind=f)                                  :: evdrop
    real(kind=f), allocatable, dimension(:)       :: evcore     ! (NELEM)
    real(kind=f), allocatable, dimension(:,:)     :: growlg     ! (NBIN,NGROUP)
    real(kind=f), allocatable, dimension(:,:)     :: evaplg     ! (NBIN,NGROUP)
    real(kind=f), allocatable, dimension(:)       :: gasprod    ! (NGAS)
    real(kind=f)                                  :: rlheat
    real(kind=f), allocatable, dimension(:,:)     :: ftoppart   ! (NBIN,NELEM)
    real(kind=f), allocatable, dimension(:,:)     :: fbotpart   ! (NBIN,NELEM)
    real(kind=f), allocatable, dimension(:,:)     :: pc_topbnd  ! (NBIN,NELEM)
    real(kind=f), allocatable, dimension(:,:)     :: pc_botbnd  ! (NBIN,NELEM)
    real(kind=f), allocatable, dimension(:,:)     :: cmf        ! (NBIN,NGROUP)
    logical, allocatable, dimension(:,:)          :: totevap    ! (NBIN,NGROUP)
    logical                                       :: too_small
    logical                                       :: too_big
    logical                                       :: nuc_small
 
    !  Coagulation kernels and bin pair mapping
    !
    !   ckernel       Coagulation kernels [cm^3/s]          {setupckern}
    !   pkernel       Coagulation production variables      {setupcoag}
    !
    real(kind=f), allocatable, dimension(:,:,:,:,:) :: ckernel ! (NZ,NBIN,NBIN,NGROUP,NGROUP)
    real(kind=f), allocatable, dimension(:,:,:,:,:,:,:) :: pkernel ! (NZ,NBIN,NBIN,NGROUP,NGROUP,NGROUP,6)

    !  Particle fall velocities and diffusivities
    !
    !   bpm       Corrects for non-sphericity and non-continuum effects {setupvfall}
    !   vf        Fall velocities at layer mid-pt                       {setupvfall}
    !   re        Reynolds' number based on <vfall>                     {setupvfall}
    !   dkz       Vert Brownian diffusion coef at layer boundary [z_units^2/s] {setupbdif}
    !   vd        Particle dry deposition velocity  [z_units/s]         {setupvdry}
    !
    real(kind=f), allocatable, dimension(:,:,:)     :: bpm        ! (NZ,NBIN,NGROUP)
    real(kind=f), allocatable, dimension(:,:,:)     :: vf         ! (NZP1,NBIN,NGROUP)
    real(kind=f), allocatable, dimension(:,:,:)     :: re         ! (NZ,NBIN,NGROUP)
    real(kind=f), allocatable, dimension(:,:,:)     :: dkz        ! (NZP1,NBIN,NGROUP)
    real(kind=f), allocatable, dimension(:,:)       :: vd         ! (NBIN,NGROUP)
    
    ! Atmospheric Structure
    !
    !  rhoa      Air density at layer mid-pt [g/x_units/y_units/z_units]  {initatm}
    !  rhoa_wet  Wet Air density averaged over grid box [g/x_units/y_units/z_units] {initatm}
    !  t         Air temperature at layer mid-pt [deg_K]                  {initatm}
    !  p         Atmospheric pressure at layer mid-pt [dyne/cm^2]         {initatm}
    !  pl        Atmospheric pressure at layer edge [dyne/cm^2]           {initatm}
    !  rmu       Air viscosity at layer mid-pt [g/cm/s]                   {initatm}
    !  thcond    Thermal conductivity of dry air [erg/cm/sec/deg_K]       {initatm}
    !  told      Temperature at beginning of time-step
    !  relhum    Hacked in relative humidity from hostmodel
    !
    real(kind=f), allocatable, dimension(:)     :: rhoa       ! (NZ)
    real(kind=f), allocatable, dimension(:)     :: rhoa_wet   ! (NZ)
    real(kind=f), allocatable, dimension(:)     :: t          ! (NZ)
    real(kind=f), allocatable, dimension(:)     :: p          ! (NZ)
    real(kind=f), allocatable, dimension(:)     :: pl         ! (NZP1)
    real(kind=f), allocatable, dimension(:)     :: rmu        ! (NZ)
    real(kind=f), allocatable, dimension(:)     :: thcond     ! (NZ)
    real(kind=f), allocatable, dimension(:)     :: told       ! (NZ)
    real(kind=f), allocatable, dimension(:)     :: relhum     ! (NZ)

    ! Condensational growth parameters
    !
    ! NOTE: Some of these variables are used for storing intermediate values in
    ! the calculations. They may no longer be necessary, when the code is
    ! implemented as F90 and values as passed as parameters between subroutines.
    !
    !   gwtmol    Molecular weight for gases [g/mol]                  {setupgrow}
    !   diffus    Diffusivity of gas in air [cm^2/s]                  {setupgrow}
    !   rlhe      Latent heat of evaporation for gas [cm^2/s^2]       {setupgrow}
    !   rlhm      Latent heat of ice melting for gas [cm^2/s^2]       {setupgrow}
    !   pvapl     Saturation vapor pressure over water [dyne/cm^2]    {vaporp}   
    !   pvapi     Saturation vapor pressure over ice [dyne/cm^2]      {vaporp}   
    !   surfctwa  Surface tension of water-air interface              {setupgkern}
    !   surfctiw  Surface tension of water-ice interface              {setupgkern}
    !   surfctia  Surface tension of ice-air interface                {setupgkern}
    !   akelvin   Exponential arg. in curvature term for growth       {setupgkern}
    !   akelvini  Curvature term for ice                              {setupgkern}
    !   ft        Ventilation factor                                  {setupgkern}
    !   gro       Growth kernel [UNITS?]                              {setupgkern}
    !   gro1      Growth kernel conduction term [UNITS?]              {setupgkern}
    !   gro2      Growth kernel radiation term [UNITS?]               {setupgkern}
    !   supsatl   Supersaturation of vapor w.r.t. liquid water [dimless]
    !   supsati   Supersaturation of vapor w.r.t. ice [dimless]                  
    !   supsatlold Supersaturation (liquid) before time-step    {prestep}
    !   supsatiold Supersaturation (ice) before time-step    {prestep}
    !   scrit     Critical supersaturation for nucleation [dimless]   {setupnuc}
    !   qrad      Particle heating rate [deg_K/s]
    real(kind=f), allocatable, dimension(:,:)    :: diffus     ! (NZ,NGAS)
    real(kind=f), allocatable, dimension(:,:)    :: rlhe       ! (NZ,NGAS)
    real(kind=f), allocatable, dimension(:,:)    :: rlhm       ! (NZ,NGAS)
    real(kind=f), allocatable, dimension(:,:)    :: pvapl      ! (NZ,NGAS)
    real(kind=f), allocatable, dimension(:,:)    :: pvapi      ! (NZ,NGAS)
    real(kind=f), allocatable, dimension(:)      :: surfctwa   ! (NZ)
    real(kind=f), allocatable, dimension(:)      :: surfctiw   ! (NZ)
    real(kind=f), allocatable, dimension(:)      :: surfctia   ! (NZ)
    real(kind=f), allocatable, dimension(:,:)    :: akelvin    ! (NZ,NGAS)
    real(kind=f), allocatable, dimension(:,:)    :: akelvini   ! (NZ,NGAS)
    real(kind=f), allocatable, dimension(:,:,:)  :: ft         ! (NZ,NBIN,NGROUP)
    real(kind=f), allocatable, dimension(:,:,:)  :: gro        ! (NZ,NBIN,NGROUP)
    real(kind=f), allocatable, dimension(:,:,:)  :: gro1       ! (NZ,NBIN,NGROUP)
    real(kind=f), allocatable, dimension(:,:)    :: gro2       ! (NZ,NGROUP)
    real(kind=f), allocatable, dimension(:,:)    :: supsatl    ! (NZ,NGAS)
    real(kind=f), allocatable, dimension(:,:)    :: supsati    ! (NZ,NGAS)
    real(kind=f), allocatable, dimension(:,:)    :: supsatlold ! (NZ,NGAS)
    real(kind=f), allocatable, dimension(:,:)    :: supsatiold ! (NZ,NGAS)
    real(kind=f), allocatable, dimension(:,:,:)  :: scrit      ! (NZ,NBIN,NGROUP)
    real(kind=f), allocatable, dimension(:,:,:)  :: qrad       ! (NZ,NBIN,NGROUP)
  end type carmastate_type   
end module