! CARMA Type aliases
! ---------------------
! This file containts shortcut names that map the variable names that
! were traditionally used in the common blocks by the Fortran 77 version
! of CARMA (globeaer.h), to the corresponding structure members in the
! Fortran 90 version of CARMA. This allows the older code to be
! converted to F90 with minimal changes, but without adding any
! processing overhead.
! ---------------------------------------------

! NOTE: Using macros causes two limitations:
!
!   1) You can not have another #define as a parameter to a macro. This causes
!      multiple expansions of the parameter. To prevent this, assign the parameter
!      to a varaible and use the variable in the macro.
!
!   2) You can not have comments on the same line as a macro. Put comments on the
!       line before the one with the macro.



#define NZ            cstate%NZ
#define NZP1          cstate%NZP1
#define NGAS          carma%NGAS
#define NBIN          carma%NBIN
#define NGROUP        carma%NGROUP
#define NELEM         carma%NELEM
#define NSOLUTE       carma%NSOLUTE

!  Model logical units for I/O
#define LUNOPRT       carma%LUNOPRT

!  Model startup control variables
#define do_print      carma%do_print

!  Gridding Information
#define igridv        cstate%igridv
#define igridh        cstate%igridh
#define xmet          cstate%xmet
#define ymet          cstate%ymet
#define zmet          cstate%zmet
#define zmetl         cstate%zmetl
#define xc            cstate%xc
#define yc            cstate%yc
#define zc            cstate%zc
#define dx            cstate%dx
#define dy            cstate%dy
#define dz            cstate%dz
#define zl            cstate%zl
#define lon           cstate%lon
#define lat           cstate%lat

! Element object
#define elemname(ielem)       carma%element(ielem)%name
#define rhoelem(ibin, ielem)  carma%element(ielem)%rho(ibin)
#define igelem(ielem)         carma%element(ielem)%igroup
#define itype(ielem)          carma%element(ielem)%itype
#define icomp(ielem)          carma%element(ielem)%icomposition
#define isolelem(ielem)       carma%element(ielem)%isolute

! Gas object
#define gasname(igas)         carma%gas(igas)%name
#define gwtmol(igas)          carma%gas(igas)%wtmol
#define ivaprtn(igas)         carma%gas(igas)%ivaprtn
#define igcomp(igas)          carma%gas(igas)%icomposition

! Group object
#define groupname(igroup)       carma%group(igroup)%name
#define nelemg(igroup)          carma%group(igroup)%nelem
#define ncore(igroup)           carma%group(igroup)%ncore
#define ishape(igroup)          carma%group(igroup)%ishape
#define ienconc(igroup)         carma%group(igroup)%ienconc
#define imomelem(igroup)        carma%group(igroup)%imomelem
#define solfac(igroup)          carma%group(igroup)%solface
#define scavcoef(igroup)        carma%group(igroup)%scavcoef
#define if_sec_mom(igroup)      carma%group(igroup)%if_sec_mom
#define is_grp_ice(igroup)      carma%group(igroup)%is_ice
#define is_grp_cloud(igroup)    carma%group(igroup)%is_cloud
#define grp_do_vtran(igroup)    carma%group(igroup)%grp_do_vtran
#define irhswell(igroup)        carma%group(igroup)%irhswell
#define irhswcomp(igroup)       carma%group(igroup)%irhswcomp
#define rmrat(igroup)           carma%group(igroup)%rmrat
#define eshape(igroup)          carma%group(igroup)%eshape
#define r(ibin,igroup)          carma%group(igroup)%r(ibin)
#define rmass(ibin,igroup)      carma%group(igroup)%rmass(ibin) 
#define vol(ibin,igroup)        carma%group(igroup)%vol(ibin)
#define dr(ibin,igroup)         carma%group(igroup)%dr(ibin)
#define dm(ibin,igroup)         carma%group(igroup)%dm(ibin)
#define rmassup(ibin,igroup)    carma%group(igroup)%rmassup(ibin)
#define rmin(igroup)            carma%group(igroup)%rmin
#define rmassmin(igroup)        carma%group(igroup)%rmassmin
#define rup(ibin,igroup)        carma%group(igroup)%rup(ibin)
#define rlow(ibin,igroup)       carma%group(igroup)%rlow(ibin)
#define icorelem(icore,igroup)  carma%group(igroup)%icorelem(icore)
#define ifallrtn(igroup)        carma%group(igroup)%ifallrtn

! Solute object
#define solname(isolute)      carma%solute(isolute)%name
#define sol_ions(isolute)     carma%solute(isolute)%ions
#define solwtmol(isolute)     carma%solute(isolute)%wtmol
#define rhosol(isolute)       carma%solute(isolute)%rho

!  Model option & control variables
#define do_cnst_rlh   carma%do_cnst_rlh
#define do_coag       carma%do_coag
#define do_detrain    carma%do_detrain
#define do_fixedinit  carma%do_fixedinit
#define do_grow       carma%do_grow
#define do_explised   carma%do_explised
#define do_incloud    carma%do_incloud
#define do_print_init carma%do_print_init
#define do_step       carma%do_step
#define do_substep    carma%do_substep
#define do_thermo     carma%do_thermo
#define do_vdiff      carma%do_vdiff
#define do_vtran      carma%do_vtran
#define if_nuc        carma%if_nuc
#define time          cstate%time
#define dtime         cstate%dtime
#define dtime_orig    cstate%dtime_orig
#define nretries      cstate%nretries
#define dtmin         carma%dtmin
#define dtmax         carma%dtmax
#define conmax        carma%conmax
#define maxsubsteps   carma%maxsubsteps
#define minsubsteps   carma%minsubsteps
#define maxretries    carma%maxretries
#define ifall         carma%ifall
#define icoagop       carma%icoagop
#define icollec       carma%icollec
#define itbnd_pc      carma%itbnd_pc
#define ibbnd_pc      carma%ibbnd_pc
#define inucgas       carma%inucgas
#define igrowgas      carma%igrowgas
#define nnuc2elem     carma%nnuc2elem
#define ievp2elem     carma%ievp2elem
#define nnucelem      carma%nnucelem
#define inucproc      carma%inucproc
#define inuc2elem     carma%inuc2elem
#define inucelem      carma%inucelem
#define inuc2bin      carma%inuc2bin
#define ievp2bin      carma%ievp2bin
#define nnucbin       carma%nnucbin
#define inucbin       carma%inucbin

#define max_nsubstep  cstate%max_nsubstep
#define max_nretry    cstate%max_nretry
#define nstep         cstate%nstep
#define nsubstep      cstate%nsubstep
#define nretry        cstate%nretry
#define zsubsteps     cstate%zsubsteps
  
!  Particle grid structure
#define diffmass      carma%diffmass
#define rhop          cstate%rhop
#define r_wet         cstate%r_wet
#define rhop_wet      cstate%rhop_wet

!  Atmospheric structure
#define rhoa          cstate%rhoa
#define rhoa_wet      cstate%rhoa_wet
#define t             cstate%t
#define p             cstate%p
#define pl            cstate%pl
#define relhum        cstate%relhum
#define told          cstate%told
#define rmu           cstate%rmu
#define thcond        cstate%thcond
#define dkz           cstate%dkz

! Model primary vars
#define pc            cstate%pc
#define pcd           cstate%pcd
#define pc_surf       cstate%pc_surf
#define gc            cstate%gc
#define cldfrc        cstate%cldfrc
#define rhcrit        cstate%rhcrit

!  Model secondary variables
#define pcl           cstate%pcl
#define gcl           cstate%gcl
#define d_gc          cstate%d_gc
#define d_t           cstate%d_t
#define dpc_sed       cstate%dpc_sed
#define pconmax       cstate%pconmax
#define coaglg        cstate%coaglg
#define coagpe        cstate%coagpe
#define rnuclg        cstate%rnuclg
#define rnucpe        cstate%rnucpe
#define pc_nucl       cstate%pc_nucl
#define growpe        cstate%growpe
#define evappe        cstate%evappe
#define coreavg       cstate%coreavg
#define coresig       cstate%coresig
#define evdrop        cstate%evdrop
#define evcore        cstate%evcore
#define growlg        cstate%growlg
#define evaplg        cstate%evaplg
#define gasprod       cstate%gasprod
#define rlheat        cstate%rlheat
#define cmf           cstate%cmf
#define totevap       cstate%totevap
#define pc_topbnd     cstate%pc_topbnd
#define pc_botbnd     cstate%pc_botbnd
#define ftoppart      cstate%ftoppart
#define fbotpart      cstate%fbotpart
#define cmf     		  cstate%cmf
#define totevap       cstate%totevap
#define too_small     cstate%too_small
#define too_big       cstate%too_big
#define nuc_small     cstate%nuc_small

!  Coagulation kernels and bin pair mapping
#define ck0           carma%ck0
#define grav_e_coll0  carma%grav_e_coll0
#define icoag         carma%icoag
#define icoagelem     carma%icoagelem
#define icoagelem_cm  carma%icoagelem_cm
#define kbin          carma%kbin

#define ckernel       cstate%ckernel
#define pkernel       cstate%pkernel

#define volx          carma%volx
#define ilow          carma%ilow
#define jlow          carma%jlow
#define iup           carma%iup
#define jup           carma%jup
#define npairl        carma%npairl
#define npairu        carma%npairu

!   Coagulation group pair mapping
#define iglow         carma%iglow
#define jglow         carma%jglow
#define igup          carma%igup
#define jgup          carma%jgup

!  Particle fall velocities, transport rates, and coagulation kernels
#define bpm           cstate%bpm
#define vf            cstate%vf
#define re            cstate%re
#define vf_const      carma%vf_const

!  Condensational growth parameters
#define diffus        cstate%diffus
#define rlhe          cstate%rlhe
#define rlhm          cstate%rlhm
#define pvapl         cstate%pvapl
#define pvapi         cstate%pvapi
#define surfctwa      cstate%surfctwa
#define surfctiw      cstate%surfctiw
#define surfctia      cstate%surfctia
#define akelvin       cstate%akelvin
#define akelvini      cstate%akelvini
#define ft            cstate%ft
#define gro           cstate%gro
#define gro1          cstate%gro1
#define gro2          cstate%gro2
#define supsatl       cstate%supsatl
#define supsati       cstate%supsati
#define supsatlold    cstate%supsatlold
#define supsatiold    cstate%supsatiold
#define scrit         cstate%scrit
#define rlh_nuc       carma%rlh_nuc
#define qrad          cstate%qrad
#define pratt         carma%pratt
#define prat          carma%prat
#define pden1         carma%pden1
#define palr          carma%palr
