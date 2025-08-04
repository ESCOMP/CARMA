! 

module test2nc_mod
        use ncio_mod
        use carma_types_mod
        use carma_enums_mod
        use nc_types_mod

        implicit none

        public ncdef
        public ncwrt_olp
        public ncwrt_ilp
        public ncdef_optics
        public ncwrt_optics
        contains

        subroutine ncdef(ncflgs, filename_out, fill_value, nlev, nlat, nlon, nbin, ngas, nelem, ngroup, &
                 sname_gas, sname_elm, sname_grp, &  
                 outid)
                
                 ! Input 
                 type(nc_type), intent(in)     :: ncflgs
                 real(kind=f), intent(in)      :: fill_value
                 integer, intent(in)           :: nlev, nlat, nlon, nbin, ngas, nelem, ngroup
                 character(LEN=*), intent(in)  :: filename_out
                 character(len=80), intent(in), optional :: sname_gas(ngas), sname_elm(nelem), sname_grp(ngroup) 
        
                ! Output
                integer, intent(out)             :: outid
                 
                ! Working variables
                integer                 :: i,ib
                character(len=80)       :: sname,vname
 
                 ! File creation and definition of dimensions
                 outid = nccreate(filename_out, "time", "time", "days since 0000-01-01 00:00:00", nlev, nlat, nlon, nbin)
                
                 ! Atmospheric structure  variables
                 if (ncflgs%f_nc_atm%f_nc_p) call ncdefreal(outid, nlev, nlat, nlon, "P", "Pressure", "Pa", fill_value)
                 if (ncflgs%f_nc_atm%f_nc_t) call ncdefreal(outid, nlev, nlat, nlon, "T", "Temperature", "K", fill_value)
                 if (ncflgs%f_nc_atm%f_nc_dt) call ncdefreal(outid, nlev, nlat, nlon, "deltaT", "T-T(1)", "K", fill_value)
                 if (ncflgs%f_nc_atm%f_nc_z) call ncdefreal(outid, nlev, nlat, nlon, "Z", "Geopotential Height", "m", fill_value)
                 if (ncflgs%f_nc_atm%f_nc_rhoa) call ncdefreal(outid, nlev, nlat, nlon, "RHOA", "Air density", "kg m-3", fill_value)
                 if (ncflgs%f_nc_atm%f_nc_rlh) call ncdefreal(outid, nlev, nlat, nlon, "RLHEAT", "Latent heat", "K/s", fill_value)

                 ! Gases variables
                 do i = 1, ngas
                        sname = sname_gas(i)
                        if (ncflgs%f_nc_gas(i)%f_nc_mmr) call ncdefreal(outid, nlev, nlat, nlon, sname, "Mass mixing ratio", "kg/kg", fill_value)
                        if (ncflgs%f_nc_gas(i)%f_nc_si) call ncdefreal(outid, nlev, nlat, nlon, trim(sname)//"SI", "Saturation ratio of vapor wrt ice", "", fill_value)
                        if (ncflgs%f_nc_gas(i)%f_nc_sl) call ncdefreal(outid, nlev, nlat, nlon, trim(sname)//"SL", "Saturation ratio of gas wrt liquid", "", fill_value)
                        if (ncflgs%f_nc_gas(i)%f_nc_ei) call ncdefreal(outid, nlev, nlat, nlon, trim(sname)//"EI", "Saturation mixing ratio of vapor wrt ice", "", fill_value)
                        if (ncflgs%f_nc_gas(i)%f_nc_el) call ncdefreal(outid, nlev, nlat, nlon, trim(sname)//"EL", "Saturation mixing ratio of vapor wrt liquid", "", fill_value)
                        if (ncflgs%f_nc_gas(i)%f_nc_wt) call ncdefreal(outid, nlev, nlat, nlon, trim(sname)//"WT", "Weight percent aerosol composition", "", fill_value)
                 end do

                 ! Element variables
                 do i = 1, nelem
                        sname = sname_elm(i)
                        do ib = 1, nbin
                                write(vname, '(A, I2.2)') trim(sname), ib
                                if (ncflgs%f_nc_elem(i)%f_nc_mr_elm_bin) call ncdefreal(outid, nlev, nlat, nlon, trim(vname)//"MR", "Mass mixing ratio", "kg/kg", fill_value)
                        end do
                 end do

                 ! Group variables 
                 do i = 1, ngroup
                        sname = sname_grp(i)
                        if (ncflgs%f_nc_group(i)%f_nc_r)       call ncdefbin(outid, nbin, trim(sname)//"_r", "Radius bin", "cm")
                        if (ncflgs%f_nc_group(i)%f_nc_rmass)   call ncdefbin(outid, nbin, trim(sname)//"_rmass", "Mass bins", "g")
                        if (ncflgs%f_nc_group(i)%f_nc_vol)     call ncdefbin(outid, nbin, trim(sname)//"_vol", "Particle volume", "cm3")
                        if (ncflgs%f_nc_group(i)%f_nc_dr)      call ncdefbin(outid, nbin, trim(sname)//"_dr", "Width of bins in radius space", "cm")
                        if (ncflgs%f_nc_group(i)%f_nc_dm)      call ncdefbin(outid, nbin, trim(sname)//"_dm", "Width of bins in mass space", "g") 
                        if (ncflgs%f_nc_group(i)%f_nc_rup)     call ncdefbin(outid, nbin, trim(sname)//"_rup", "Upper bin boundary radius", "cm")
                        if (ncflgs%f_nc_group(i)%f_nc_rlw)     call ncdefbin(outid, nbin, trim(sname)//"_rlow", "Lower bin boundary radius", "cm")
                        if (ncflgs%f_nc_group(i)%f_nc_rrat)    call ncdefbin(outid, nbin, trim(sname)//"_rrat", "Ratio of maximum diameter to diameter of equivalent sphere", "1")
                        if (ncflgs%f_nc_group(i)%f_nc_rprat)   call ncdefbin(outid, nbin, trim(sname)//"_rprat", "Ratio of mobility diameter of a porous particle to diameter of equivlent sphere", "")
                        if (ncflgs%f_nc_group(i)%f_nc_arat)    call ncdefbin(outid, nbin, trim(sname)//"_arat", "Ratio of projected area to projected area of containing sphere", "")
                        if (ncflgs%f_nc_group(i)%f_nc_df)      call ncdefbin(outid, nbin, trim(sname)//"_df", "Fractal dimension", "")
                        if (ncflgs%f_nc_group(i)%f_nc_nmon)    call ncdefbin(outid, nbin, trim(sname)//"_nmon", "Number of monomers", "")
                 
                        if (ncflgs%f_nc_group(i)%f_nc_nd)  call ncdefreal(outid, nlev, nlat, nlon, trim(sname)//"ND", "Number density", "#/cm3", fill_value)
                        if (ncflgs%f_nc_group(i)%f_nc_ad)  call ncdefreal(outid, nlev, nlat, nlon, trim(sname)//"AD", "Surface area density", "um2/cm3", fill_value)
                        if (ncflgs%f_nc_group(i)%f_nc_md)  call ncdefreal(outid, nlev, nlat, nlon, trim(sname)//"MD", "Mass density", "g/cm3", fill_value)
                        if (ncflgs%f_nc_group(i)%f_nc_re)  call ncdefreal(outid, nlev, nlat, nlon, trim(sname)//"RE", "Dry effective radius", "um", fill_value)
                        if (ncflgs%f_nc_group(i)%f_nc_rew) call ncdefreal(outid, nlev, nlat, nlon, trim(sname)//"REw", "Wet effective radius", "um", fill_value)
                        if (ncflgs%f_nc_group(i)%f_nc_rm)  call ncdefreal(outid, nlev, nlat, nlon, trim(sname)//"RM", "Mitchell effective radius", "um", fill_value)
                        if (ncflgs%f_nc_group(i)%f_nc_jn)  call ncdefreal(outid, nlev, nlat, nlon, trim(sname)//"JN", "Nucleation rate", "#/cm3/s", fill_value)
                        if (ncflgs%f_nc_group(i)%f_nc_mr)  call ncdefreal(outid, nlev, nlat, nlon, trim(sname)//"MR", "Mass mixing ratio", "kg/kg", fill_value)
                        if (ncflgs%f_nc_group(i)%f_nc_pa)  call ncdefreal(outid, nlev, nlat, nlon, trim(sname)//"PA", "Projected area", "1/cm", fill_value)
                        if (ncflgs%f_nc_group(i)%f_nc_ar)  call ncdefreal(outid, nlev, nlat, nlon, trim(sname)//"AR", "Area ratio", "", fill_value)
                        if (ncflgs%f_nc_group(i)%f_nc_vm)  call ncdefreal(outid, nlev, nlat, nlon, trim(sname)//"VM", "Fall velocity", "cm/s", fill_value) 
                        if (ncflgs%f_nc_group(i)%f_nc_exvis) call ncdefreal(outid, nlev, nlat, nlon, trim(sname)//"EX_VIS", "Extinction 550 nm", "km-1", fill_value)
                        if (ncflgs%f_nc_group(i)%f_nc_odvis) call ncdefreal(outid, nlev, nlat, nlon, trim(sname)//"OD_VIS", "Optical depth 550 nm", "", fill_value)


                        do ib = 1, nbin
                                write(vname, '(A, I2.2)') trim(sname), ib
                                if (ncflgs%f_nc_group(i)%f_nc_vd_bin) call ncdefreal(outid, nlev, nlat, nlon, trim(vname)//"VD", "Dry deposition velocity", "cm/s", fill_value)
                                if (ncflgs%f_nc_group(i)%f_nc_nd_bin) call ncdefreal(outid, nlev, nlat, nlon, trim(vname)//"ND", "Number density", "#/cm3", fill_value)
                                if (ncflgs%f_nc_group(i)%f_nc_wr_bin) call ncdefreal(outid, nlev, nlat, nlon, trim(vname)//"WR", "Wet radius", "cm", fill_value)
                                if (ncflgs%f_nc_group(i)%f_nc_ro_bin) call ncdefreal(outid, nlev, nlat, nlon, trim(vname)//"RO", "Wet particle density", "g/cm3", fill_value)
                                if (ncflgs%f_nc_group(i)%f_nc_mr_bin) call ncdefreal(outid, nlev, nlat, nlon, trim(vname)//"MR", "Mass mixing ratio", "kg/kg", fill_value)
                        end do

                end do
               
               
                call ncdefend(outid) 
        end subroutine ncdef


        subroutine ncwrt_olp(ncflgs, outid, nlev, nlat, nlon, ngas, nelem, ngroup, nbin, & 
                       lev, lat, lon, bins, sname_grp, &
                       r, rmass, vol, dr, dm, rup, rlow, rrat, arat, rprat,  df, nmon)
                
                ! Input 
                type(nc_type), intent(in)     :: ncflgs
                integer, intent(in)           :: outid, nlev, nlat, nlon, nbin, ngas, nelem, ngroup
                integer, intent(in)           :: bins(nbin)
                real(kind=f), intent(in)      :: lev(nlev), lat(nlat), lon(nlon)
                character(len=80), intent(in)    :: sname_grp(ngroup)
                
                real(kind=f), intent(in), optional  :: r(nbin,ngroup), rmass(nbin,ngroup), vol(nbin,ngroup)
                real(kind=f), intent(in), optional  :: dr(nbin,ngroup), dm(nbin,ngroup)
                real(kind=f), intent(in), optional  :: rup(nbin,ngroup), rlow(nbin,ngroup), rrat(nbin,ngroup)
                real(kind=f), intent(in), optional  :: rprat(nbin,ngroup), arat(nbin,ngroup)
                real(kind=f), intent(in), optional  :: df(nbin,ngroup), nmon(nbin,ngroup)


                ! Working variables
                integer                 :: i,ib
                character(len=80)       :: sname 

               
                call ncwrtdims(outid, nlev, lev, nlat, lat, nlon, lon, nbin, bins)
                
                do i = 1, ngroup
                        sname = sname_grp(i) 
                        if (ncflgs%f_nc_group(i)%f_nc_r)        call ncwrtbin(outid, trim(sname)//"_r", nbin, r(:,i) )
                        if (ncflgs%f_nc_group(i)%f_nc_rmass)    call ncwrtbin(outid, trim(sname)//"_rmass", nbin, rmass(:,i) )
                        if (ncflgs%f_nc_group(i)%f_nc_vol)      call ncwrtbin(outid, trim(sname)//"_vol", nbin, vol(:,i) )
                        if (ncflgs%f_nc_group(i)%f_nc_dr)       call ncwrtbin(outid, trim(sname)//"_dr", nbin, dr(:,i) )
                        if (ncflgs%f_nc_group(i)%f_nc_dm)       call ncwrtbin(outid, trim(sname)//"_dm", nbin, dm(:,i) ) 
                        if (ncflgs%f_nc_group(i)%f_nc_rup)      call ncwrtbin(outid, trim(sname)//"_rup", nbin, rup(:,i) )
                        if (ncflgs%f_nc_group(i)%f_nc_rlw)      call ncwrtbin(outid, trim(sname)//"_rlow", nbin, rlow(:,i) )
                        if (ncflgs%f_nc_group(i)%f_nc_rrat)     call ncwrtbin(outid, trim(sname)//"_rrat", nbin, rrat(:,i) )
                        if (ncflgs%f_nc_group(i)%f_nc_rprat)    call ncwrtbin(outid, trim(sname)//"_rprat", nbin, rprat(:,i) )
                        if (ncflgs%f_nc_group(i)%f_nc_arat)     call ncwrtbin(outid, trim(sname)//"_arat", nbin, arat(:,i) )
                        if (ncflgs%f_nc_group(i)%f_nc_df)       call ncwrtbin(outid, trim(sname)//"_df", nbin, df(:,i) )
                        if (ncflgs%f_nc_group(i)%f_nc_nmon)     call ncwrtbin(outid, trim(sname)//"_nmon", nbin, nmon(:,i) )  
                end do
        end subroutine ncwrt_olp


        subroutine ncwrt_ilp(ncflgs, outid, nlev, nlat, nlon, nbin, ngas, nelem, ngroup, istep,  &
                        lev, lat, lon, bins, time, secsperday, date, &
                        sname_gas, sname_elm, sname_grp, &
                        p, t, z, dT, rhoa, rlheat, &
                        mmr, si, sl, ei, el, wt, &
                        nd, ad, md, re, rew, rm, jn, mr, pa, ar, vm, ex_vis, od_vis, &
                        vd_bin, nd_bin, wr_bin, ro_bin, mr_bin, mr_elm_bin)


                ! Input
                type(nc_type), intent(in)               :: ncflgs
                integer, intent(in)                     :: outid, istep, date, nlev, nlat, nlon, nbin, ngas, nelem, ngroup, lat, lon
                integer, intent(in)                     :: bins(nbin)
                real(kind=f), intent(in)                :: time, secsperday
                real(kind=f), intent(in)                :: lev(nlev)
                character(len=80), intent(in), optional :: sname_gas(ngas), sname_elm(nelem), sname_grp(ngroup)

                
                real(kind=f), intent(in), optional :: p(nlev), t(nlev), z(nlev), dT(nlev), rhoa(nlev), rlheat(nlev)
                real(kind=f), intent(in), optional :: mmr(nlev,ngas), si(nlev,ngas), sl(nlev,ngas), ei(nlev,ngas)
                real(kind=f), intent(in), optional :: el(nlev,ngas), wt(nlev,ngas)
                real(kind=f), intent(in), optional :: nd(nlev,ngroup), ad(nlev,ngroup), md(nlev,ngroup), re(nlev,ngroup), rew(nlev,ngroup)
                real(kind=f), intent(in), optional :: rm(nlev,ngroup), jn(nlev,ngroup), mr(nlev,ngroup), pa(nlev,ngroup), ar(nlev,ngroup)
                real(kind=f), intent(in), optional :: vm(nlev,ngroup), od_vis(nlev,ngroup), ex_vis(nlev,ngroup)

                real(kind=f), intent(in), optional :: vd_bin(nlev,ngroup,nbin), nd_bin(nlev,ngroup,nbin), wr_bin(nlev,ngroup,nbin)
                real(kind=f), intent(in), optional :: ro_bin(nlev,ngroup,nbin), mr_bin(nlev,ngroup,nbin)

                real(kind=f), intent(in), optional :: mr_elm_bin(nlev,nelem,nbin)
                
                ! Working variables
                integer                 :: i,ib
                character(len=80)       :: sname, vname

                ! Time variable
                call ncwrttime(outid, istep, time / secsperday, date)
                
                ! Atmospheric variables
                if (ncflgs%f_nc_atm%f_nc_p)    call ncwrtreal(outid, "P", istep, nlev, nlat, nlon, lat, lon, p)
                if (ncflgs%f_nc_atm%f_nc_t)    call ncwrtreal(outid, "T", istep, nlev, nlat, nlon, lat, lon, t)
                if (ncflgs%f_nc_atm%f_nc_z)    call ncwrtreal(outid, "Z", istep, nlev, nlat, nlon, lat, lon, z)
                if (ncflgs%f_nc_atm%f_nc_dt)   call ncwrtreal(outid, "deltaT", istep, nlev, nlat, nlon, lat, lon, dT)
                if (ncflgs%f_nc_atm%f_nc_rhoa) call ncwrtreal(outid, "RHOA", istep, nlev, nlat, nlon, lat, lon, rhoa)
                if (ncflgs%f_nc_atm%f_nc_rlh)  call ncwrtreal(outid, "RLHEAT", istep, nlev, nlat, nlon, lat, lon, rlheat)

                ! Gases
                do i = 1, ngas
                        sname = sname_gas(i)
                        if (ncflgs%f_nc_gas(i)%f_nc_mmr) call ncwrtreal(outid, sname, istep, nlev, nlat, nlon, lat, lon,  mmr(:,i))
                        if (ncflgs%f_nc_gas(i)%f_nc_si)  call ncwrtreal(outid, trim(sname)//"SI", istep, nlev, nlat, nlon, lat, lon, si(:,i))
                        if (ncflgs%f_nc_gas(i)%f_nc_sl)  call ncwrtreal(outid, trim(sname)//"SL", istep, nlev, nlat, nlon, lat, lon, sl(:,i))
                        if (ncflgs%f_nc_gas(i)%f_nc_ei)  call ncwrtreal(outid, trim(sname)//"EI", istep, nlev, nlat, nlon, lat, lon, ei(:,i))
                        if (ncflgs%f_nc_gas(i)%f_nc_el)  call ncwrtreal(outid, trim(sname)//"EL", istep, nlev, nlat, nlon, lat, lon, el(:,i))
                        if (ncflgs%f_nc_gas(i)%f_nc_wt)  call ncwrtreal(outid, trim(sname)//"WT", istep, nlev, nlat, nlon, lat, lon, wt(:,i))
                end do

                 ! Element variables
                 do i = 1, nelem
                        sname = sname_elm(i)
                        do ib = 1, nbin
                                write(vname, '(A, I2.2)') trim(sname), ib
                               if (ncflgs%f_nc_elem(i)%f_nc_mr_elm_bin) call ncwrtreal(outid, trim(vname)//"MR", istep, nlev, nlat, nlon, lat, lon, mr_elm_bin(:,i,ib)) 
                        end do
                 end do

               ! Groups variables 
               do i = 1, ngroup
                        sname = sname_grp(i)
                        if (ncflgs%f_nc_group(i)%f_nc_nd)    call ncwrtreal(outid, trim(sname)//"ND", istep, nlev, nlat, nlon, lat, lon, nd(:,i))
                        if (ncflgs%f_nc_group(i)%f_nc_ad)    call ncwrtreal(outid, trim(sname)//"AD", istep, nlev, nlat, nlon, lat, lon, ad(:,i))
                        if (ncflgs%f_nc_group(i)%f_nc_md)    call ncwrtreal(outid, trim(sname)//"MD", istep, nlev, nlat, nlon, lat, lon, md(:,i))
                        if (ncflgs%f_nc_group(i)%f_nc_re)    call ncwrtreal(outid, trim(sname)//"RE", istep, nlev, nlat, nlon, lat, lon, re(:,i))
                        if (ncflgs%f_nc_group(i)%f_nc_rew)   call ncwrtreal(outid, trim(sname)//"REw", istep, nlev, nlat, nlon, lat, lon, rew(:,i))
                        if (ncflgs%f_nc_group(i)%f_nc_rm)    call ncwrtreal(outid, trim(sname)//"RM", istep, nlev, nlat, nlon, lat, lon, rm(:,i))
                        if (ncflgs%f_nc_group(i)%f_nc_jn)    call ncwrtreal(outid, trim(sname)//"JN", istep, nlev, nlat, nlon, lat, lon, jn(:,i))
                        if (ncflgs%f_nc_group(i)%f_nc_mr)    call ncwrtreal(outid, trim(sname)//"MR", istep, nlev, nlat, nlon, lat, lon, mr(:,i))
                        if (ncflgs%f_nc_group(i)%f_nc_pa)    call ncwrtreal(outid, trim(sname)//"PA", istep, nlev, nlat, nlon, lat, lon, pa(:,i))
                        if (ncflgs%f_nc_group(i)%f_nc_ar)    call ncwrtreal(outid, trim(sname)//"AR", istep, nlev, nlat, nlon, lat, lon, ar(:,i))
                        if (ncflgs%f_nc_group(i)%f_nc_vm)    call ncwrtreal(outid, trim(sname)//"VM", istep, nlev, nlat, nlon, lat, lon, vm(:,i)) 
                        if (ncflgs%f_nc_group(i)%f_nc_exvis) call ncwrtreal(outid, trim(sname)//"EX_VIS", istep, nlev, nlat, nlon, lat, lon, ex_vis(:,i))
                        if (ncflgs%f_nc_group(i)%f_nc_odvis) call ncwrtreal(outid, trim(sname)//"OD_VIS", istep, nlev, nlat, nlon, lat, lon, od_vis(:,i))

                        do ib = 1, nbin
 
                                write(vname, '(A, I2.2)') trim(sname), ib
                                if (ncflgs%f_nc_group(i)%f_nc_vd_bin) call ncwrtreal(outid, trim(vname)//"VD", istep, nlev, nlat, nlon, lat, lon, vd_bin(:,i,ib))
                                if (ncflgs%f_nc_group(i)%f_nc_nd_bin) call ncwrtreal(outid, trim(vname)//"ND", istep, nlev, nlat, nlon, lat, lon, nd_bin(:,i,ib))
                                if (ncflgs%f_nc_group(i)%f_nc_wr_bin) call ncwrtreal(outid, trim(vname)//"WR", istep, nlev, nlat, nlon, lat, lon, wr_bin(:,i,ib))
                                if (ncflgs%f_nc_group(i)%f_nc_ro_bin) call ncwrtreal(outid, trim(vname)//"RO", istep, nlev, nlat, nlon, lat, lon, ro_bin(:,i,ib))
                                if (ncflgs%f_nc_group(i)%f_nc_mr_bin) call ncwrtreal(outid, trim(vname)//"MR", istep, nlev, nlat, nlon, lat, lon, mr_bin(:,i,ib))
                        end do
               end do

       end subroutine ncwrt_ilp

       subroutine ncdef_optics( filename_out, nwave, nbin, nrun, sname_run, outid)

                ! Input

                integer, intent(in)           :: nwave, nbin, nrun
                character(LEN=*), intent(in)  :: filename_out
                character(len=80), intent(in) :: sname_run(nrun)

                ! Output
                integer, intent(out)             :: outid

                ! Working variables
                integer                 :: i
                character(len=80)       :: sname


                ! File creation and definition of dimensions
                outid = nccreate_optics(filename_out, nwave, nbin)
                
                call ncdefreal_optics(outid, nwave, 1, "refidx_r", "Refractive index - real", "")
                call ncdefreal_optics(outid, nwave, 1, "refidx_i", "Refractive index - imaginary", "")
                 do i = 1, nrun
                        sname = sname_run(i)
                        if (nbin.gt.1) then
                                call ncdefbin(outid, nbin, trim(sname)//"_r", "Bin radius", "cm") 
                                call ncdefbin(outid, nbin, trim(sname)//"_rmass", "Mass bins", "g")
                                call ncdefbin(outid, nbin, trim(sname)//"_vol", "Particle volume", "cm3")
                                call ncdefbin(outid, nbin, trim(sname)//"_dr", "Width of bins in radius space", "cm")
                                call ncdefbin(outid, nbin, trim(sname)//"_dm", "Width of bins in mass space", "g")
                                call ncdefbin(outid, nbin, trim(sname)//"_rrat", "Ratio of maximum diameter to diameter of equivalent sphere", "1") 
                                call ncdefbin(outid, nbin, trim(sname)//"_rprat", "Ratio of mobility diameter of a porous particle to diameter of equivlent sphere", "1") 
                                call ncdefbin(outid, nbin, trim(sname)//"_arat", "Ratio of projected area to projected area of containing sphere", "1")
                                call ncdefbin(outid, nbin, trim(sname)//"_nmon", "number of monomers", "")
                                call ncdefbin(outid, nbin, trim(sname)//"_df", "fractal dimension", "")
                                call ncdefbin(outid, nbin, trim(sname)//"_rg", "Bin radius for agglomerates", "cm")
                                call ncdefbin(outid, nbin, trim(sname)//"_rp", "Bin mobility radius for agglomerates", "cm")
                        end if
                        call ncdefreal_optics(outid, nwave, nbin, trim(sname)//"_qext", "Extinction efficiency", "")
                        call ncdefreal_optics(outid, nwave, nbin, trim(sname)//"_asy", "Asymmetry factor", "")
                        call ncdefreal_optics(outid, nwave, nbin, trim(sname)//"_ssa", "Single scattering albedo", "")
                end do
                
                call ncdefend(outid) 
                
        end subroutine ncdef_optics


        subroutine ncwrt_optics(outid, nwave, nbin, nrun, wave, bins, sname_run,  & 
                r_refidx, i_refidx, r, rmass, vol, dr, dm, rrat, rprat, arat, nmon, rg, rp, df, qext, asy, ssa)

        ! Input
        integer, intent(in)           :: outid,nwave, nbin, nrun
        integer, intent(in)           :: bins(nbin)
        real(kind=f), intent(in)      :: wave(nwave), r_refidx(nwave), i_refidx(nwave)
        character(len=80), intent(in) :: sname_run(nrun)

        real(kind=f), intent(in), optional :: r(nbin,nrun), rmass(nbin,nrun), vol(nbin,nrun), dr(nbin,nrun), dm(nbin,nrun) 
        real(kind=f), intent(in), optional :: rrat(nbin,nrun), rprat(nbin,nrun), arat(nbin,nrun), nmon(nbin,nrun), df(nbin,nrun)
        real(kind=f), intent(in), optional :: rg(nbin,nrun), rp(nbin,nrun)
        real(kind=f), intent(in), optional :: qext(nwave, nbin, nrun), asy(nwave, nbin, nrun), ssa(nwave, nbin, nrun)

        ! Working variables
        integer                 :: i
        character(len=80)       :: sname

        call ncwrtdims_optics(outid, nwave, wave, nbin, bins)
        call ncwrtreal_optics(outid, "refidx_r", nwave, 1, r_refidx(:))
        call ncwrtreal_optics(outid, "refidx_i", nwave, 1, i_refidx(:))

        do i = 1, nrun
                sname = sname_run(i)
                if (present(r))     call ncwrtbin(outid, trim(sname)//"_r", nbin, r(:,i)) 
                if (present(rmass)) call ncwrtbin(outid, trim(sname)//"_rmass", nbin, rmass(:,i) )
                if (present(vol))   call ncwrtbin(outid, trim(sname)//"_vol", nbin, vol(:,i) )
                if (present(dr))    call ncwrtbin(outid, trim(sname)//"_dr", nbin, dr(:,i) )
                if (present(dm))    call ncwrtbin(outid, trim(sname)//"_dm", nbin, dm(:,i) ) 
                if (present(rrat))  call ncwrtbin(outid, trim(sname)//"_rrat", nbin, rrat(:,i) )
                if (present(rprat))  call ncwrtbin(outid, trim(sname)//"_rprat", nbin, rprat(:,i) )
                if (present(arat))  call ncwrtbin(outid, trim(sname)//"_arat", nbin, arat(:,i) )
                if (present(nmon))  call ncwrtbin(outid, trim(sname)//"_nmon", nbin, nmon(:,i) )
                if (present(rg))    call ncwrtbin(outid, trim(sname)//"_rg", nbin, rg(:,i))
                if (present(rp))    call ncwrtbin(outid, trim(sname)//"_rp", nbin, rp(:,i))
                if (present(df))    call ncwrtbin(outid, trim(sname)//"_df", nbin, df(:,i) )
                
                if (present(qext))  call ncwrtreal_optics(outid, trim(sname)//"_qext", nwave, nbin, qext(:,:,i))
                if (present(asy))   call ncwrtreal_optics(outid, trim(sname)//"_asy",  nwave, nbin, asy(:,:,i))
                if (present(ssa))   call ncwrtreal_optics(outid, trim(sname)//"_ssa",  nwave, nbin, ssa(:,:,i))
        end do
        end subroutine ncwrt_optics

end module
