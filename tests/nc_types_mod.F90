module nc_types_mod
        implicit none

        ! All variables and procedures are private except those explicitly declared to be public.
        private
        public nctype_create
        public nc_type 
   
        type :: ncatm_type   
                logical :: f_nc_p
                logical :: f_nc_t
                logical :: f_nc_dt
                logical :: f_nc_z
                logical :: f_nc_rhoa
                logical :: f_nc_rlh
        end type ncatm_type


        type :: ncgas_type
                logical :: f_nc_mmr
                logical :: f_nc_sl
                logical :: f_nc_si
                logical :: f_nc_ei
                logical :: f_nc_el
                logical :: f_nc_wt
        end type ncgas_type

        type :: ncgroup_type
                logical :: f_nc_r 
                logical :: f_nc_rmass
                logical :: f_nc_vol
                logical :: f_nc_dr
                logical :: f_nc_dm 
                logical :: f_nc_rup
                logical :: f_nc_rlw
                logical :: f_nc_rrat
                logical :: f_nc_rprat
                logical :: f_nc_arat
                logical :: f_nc_df
                logical :: f_nc_nmon

                logical :: f_nc_nd
                logical :: f_nc_ad
                logical :: f_nc_md
                logical :: f_nc_re
                logical :: f_nc_rew
                logical :: f_nc_rm
                logical :: f_nc_jn
                logical :: f_nc_mr
                logical :: f_nc_pa
                logical :: f_nc_ar
                logical :: f_nc_vm
                logical :: f_nc_odvis
                logical :: f_nc_exvis

                logical :: f_nc_wr_bin
                logical :: f_nc_nd_bin
                logical :: f_nc_ro_bin
                logical :: f_nc_mr_bin
                logical :: f_nc_vd_bin
        end type ncgroup_type

        type :: ncelem_type
                logical :: f_nc_mr_elm_bin
        end type ncelem_type

        type :: nc_type
                type(ncatm_type),     allocatable                  :: f_nc_atm
                type(ncgas_type),     allocatable, dimension(:)    :: f_nc_gas       
                type(ncelem_type),    allocatable, dimension(:)    :: f_nc_elem       
                type(ncgroup_type),   allocatable, dimension(:)    :: f_nc_group    
        
                contains
                procedure :: nctype_create
        
        end type nc_type


        contains 
        subroutine nctype_create(ncflgs,NGAS, NELEM, NGROUP)
                class(nc_type), intent(out) :: ncflgs  
                integer, intent(in)         :: NGAS
                integer, intent(in)         :: NELEM
                integer, intent(in)         :: NGROUP

                integer :: igas, ielem, igroup

                ! Allocate tables
                allocate( &
                        ncflgs%f_nc_atm, &
                        ncflgs%f_nc_gas(NGAS), &
                        ncflgs%f_nc_elem(NELEM), &
                        ncflgs%f_nc_group(NGROUP) &
                        )

                ! Atm flags 
                        ncflgs%f_nc_atm%f_nc_p    = .FALSE.
                        ncflgs%f_nc_atm%f_nc_t    = .FALSE.
                        ncflgs%f_nc_atm%f_nc_dt   = .FALSE.
                        ncflgs%f_nc_atm%f_nc_z    = .FALSE.
                        ncflgs%f_nc_atm%f_nc_rhoa = .FALSE.
                        ncflgs%f_nc_atm%f_nc_rlh  = .FALSE.


                ! GAS flags
                if (NGAS > 0) then
                        do igas = 1, NGAS
                                ncflgs%f_nc_gas(igas)%f_nc_mmr = .FALSE.
                                ncflgs%f_nc_gas(igas)%f_nc_sl  = .FALSE.
                                ncflgs%f_nc_gas(igas)%f_nc_si  = .FALSE.
                                ncflgs%f_nc_gas(igas)%f_nc_ei  = .FALSE.
                                ncflgs%f_nc_gas(igas)%f_nc_el  = .FALSE.
                                ncflgs%f_nc_gas(igas)%f_nc_wt  = .FALSE.
                        end do
                end if

                ! GROUP flags
                if (NGROUP > 0) then
                        do igroup = 1, NGROUP
                                ncflgs%f_nc_group(igroup)%f_nc_r       = .FALSE.
                                ncflgs%f_nc_group(igroup)%f_nc_rmass   = .FALSE.
                                ncflgs%f_nc_group(igroup)%f_nc_vol     = .FALSE.
                                ncflgs%f_nc_group(igroup)%f_nc_dr      = .FALSE.
                                ncflgs%f_nc_group(igroup)%f_nc_dm      = .FALSE. 
                                ncflgs%f_nc_group(igroup)%f_nc_rup     = .FALSE.
                                ncflgs%f_nc_group(igroup)%f_nc_rlw     = .FALSE.
                                ncflgs%f_nc_group(igroup)%f_nc_rrat    = .FALSE.
                                ncflgs%f_nc_group(igroup)%f_nc_rprat   = .FALSE.
                                ncflgs%f_nc_group(igroup)%f_nc_arat    = .FALSE.
                                ncflgs%f_nc_group(igroup)%f_nc_df      = .FALSE.
                                ncflgs%f_nc_group(igroup)%f_nc_nmon    = .FALSE.

                                ncflgs%f_nc_group(igroup)%f_nc_nd = .FALSE.
                                ncflgs%f_nc_group(igroup)%f_nc_ad = .FALSE.
                                ncflgs%f_nc_group(igroup)%f_nc_md = .FALSE.
                                ncflgs%f_nc_group(igroup)%f_nc_re = .FALSE.
                                ncflgs%f_nc_group(igroup)%f_nc_rew = .FALSE.
                                ncflgs%f_nc_group(igroup)%f_nc_rm = .FALSE.
                                ncflgs%f_nc_group(igroup)%f_nc_jn = .FALSE.
                                ncflgs%f_nc_group(igroup)%f_nc_mr = .FALSE.
                                ncflgs%f_nc_group(igroup)%f_nc_pa = .FALSE.
                                ncflgs%f_nc_group(igroup)%f_nc_ar = .FALSE.
                                ncflgs%f_nc_group(igroup)%f_nc_vm = .FALSE.
                                ncflgs%f_nc_group(igroup)%f_nc_exvis = .FALSE.
                                ncflgs%f_nc_group(igroup)%f_nc_odvis = .FALSE. 
                                
                                ncflgs%f_nc_group(igroup)%f_nc_wr_bin = .FALSE. 
                                ncflgs%f_nc_group(igroup)%f_nc_nd_bin = .FALSE.
                                ncflgs%f_nc_group(igroup)%f_nc_ro_bin = .FALSE.
                                ncflgs%f_nc_group(igroup)%f_nc_mr_bin = .FALSE.
                                ncflgs%f_nc_group(igroup)%f_nc_vd_bin = .FALSE. 
                        end do
                end if

                if (NELEM > 0) then
                        do ielem = 1, NELEM
                                ncflgs%f_nc_elem(ielem)%f_nc_mr_elm_bin = .FALSE.
                        end do
                end if
        end subroutine nctype_create

end module nc_types_mod
