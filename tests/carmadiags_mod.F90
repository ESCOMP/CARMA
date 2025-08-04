module carmadiags_mod
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmaelement_mod
  use carmagroup_mod
  use carmastate_mod

  implicit none
  public carmadiags

  contains

  subroutine carmadiags(cstate, carma, nlev, nbin, nelem, ngroup, nwave, &
                dz, qext, idx_wave, &
                nd, ad, md, re, rew, rm, jn, mr, pa, ar, vm, ex, od,&
                wr_bin, nd_bin, ro_bin, mr_bin, vd_bin, mr_elm_bin)
 
  ! Input 
  type(carma_type), intent(in)       :: carma
  type(carmastate_type), intent(in)  :: cstate
  integer, intent(in)                :: nlev, nbin, nelem, ngroup, nwave
  integer, intent(in), optional      :: idx_wave
  real(kind=f), intent(in)           :: dz 
  real(kind=f), intent(in), optional :: qext(nwave, nbin, ngroup)

  ! Output
  real(kind=f), intent(out), optional :: nd(nlev, ngroup) 
  real(kind=f), intent(out), optional :: ad(nlev, ngroup) 
  real(kind=f), intent(out), optional :: md(nlev, ngroup)
  real(kind=f), intent(out), optional :: re(nlev, ngroup)
  real(kind=f), intent(out), optional :: rew(nlev, ngroup)
  real(kind=f), intent(out), optional :: rm(nlev, ngroup)
  real(kind=f), intent(out), optional :: jn(nlev, ngroup)
  real(kind=f), intent(out), optional :: mr(nlev, ngroup)
  real(kind=f), intent(out), optional :: pa(nlev, ngroup)
  real(kind=f), intent(out), optional :: ar(nlev, ngroup)
  real(kind=f), intent(out), optional :: vm(nlev, ngroup) 
  real(kind=f), intent(out), optional :: ex(nlev, ngroup)
  real(kind=f), intent(out), optional :: od(nlev, ngroup)

  real(kind=f), intent(out), optional :: wr_bin(nlev, ngroup, nbin)
  real(kind=f), intent(out), optional :: nd_bin(nlev, ngroup, nbin)
  real(kind=f), intent(out), optional :: ro_bin(nlev, ngroup, nbin)
  real(kind=f), intent(out), optional :: mr_bin(nlev, ngroup, nbin)
  real(kind=f), intent(out), optional :: vd_bin(nlev, ngroup, nbin)
  real(kind=f), intent(out), optional :: mr_elm_bin(nlev, nelem, nbin)
  
  ! Local Variables
  integer      :: igroup, ielem, ibin, ienconc, rc, cnsttype, maxbin
  real(kind=f) :: r(nbin), rmass(nbin), rrat(nbin), arat(nbin)
  logical      :: is_cloud, is_ice, grp_do_drydep
  real(kind=f) :: vqext(nbin, ngroup)

  real(kind=f)           :: mmr(nlev)             !! the bin mass mixing ratio [kg/kg]
  real(kind=f)           :: totalmmr(nlev)        !! mmr of the entire particle (kg/m3)
  real(kind=f)           :: numberDensity(nlev)   !! number density [#/cm3]
  real(kind=f)           :: nucleationRate(nlev)  !! nucleation rate [1/cm3/s]
  real(kind=f)           :: rwet(nlev)           !! wet particle radius [cm]
  real(kind=f)           :: rhop_wet(nlev)        !! wet particle density [g/cm3]
  real(kind=f)           :: dd                    !! particle sedimentation mass flux to surface [kg/m2/s] 
  real(kind=f)           :: vd                    !! deposition velocity [cm/s]
  real(kind=f)           :: vf(nlev+1)            !! fall velocity [cm/s]
  real(kind=f)           :: dtpart(nlev)          !! delta particle temperature [K]
  real(kind=f)           :: re2(nlev, ngroup), re3(nlev, ngroup)
  real(kind=f)           :: rew2(nlev, ngroup), rew3(nlev, ngroup)

  
  real(kind=f) :: vnd(nlev, ngroup)
  real(kind=f) :: vad(nlev, ngroup) 
  real(kind=f) :: vmd(nlev, ngroup)
  real(kind=f) :: vre(nlev, ngroup)
  real(kind=f) :: vrew(nlev, ngroup)
  real(kind=f) :: vrm(nlev, ngroup)
  real(kind=f) :: vjn(nlev, ngroup)
  real(kind=f) :: vmr(nlev, ngroup)
  real(kind=f) :: vpa(nlev, ngroup)
  real(kind=f) :: var(nlev, ngroup)
  real(kind=f) :: vvm(nlev, ngroup) 
  real(kind=f) :: vex(nlev, ngroup)
  real(kind=f) :: vod(nlev, ngroup)

  real(kind=f) :: vwr_bin(nlev, ngroup, nbin)
  real(kind=f) :: vnd_bin(nlev, ngroup, nbin)
  real(kind=f) :: vro_bin(nlev, ngroup, nbin)
  real(kind=f) :: vmr_bin(nlev, ngroup, nbin)
  real(kind=f) :: vvd_bin(nlev, ngroup, nbin)
  real(kind=f) :: vmr_elm_bin(nlev, nelem, nbin)
  
  ! Initialization
  re2(:, :) = 0._f
  re3(:, :) = 0._f
  rew2(:, :) = 0._f
  rew3(:, :) = 0._f

  vnd(:, :)   = 0._f
  vad(:, :)   = 0._f
  vmd(:, :)   = 0._f
  vre(:, :)   = 0._f
  vrew(:, :)  = 0._f
  vrm(:, :)   = 0._f
  vjn(:, :)   = 0._f
  vmr(:, :)   = 0._f
  vpa(:, :)   = 0._f
  var(:, :)   = 0._f
  vvm(:, :)   = 0._f 
  vex(:, :)   = 0._f
  vod(:, :)   = 0._f 
  vwr_bin(:, :, :) = 0._f
  vnd_bin(:, :, :) = 0._f
  vro_bin(:, :, :) = 0._f 
  vvd_bin(:, :, :) = 0._f
  vmr_elm_bin(:, :, :) = 0._f
  
 
  vqext(:, :) = 2.0_f
  if (present(qext)) vqext(:, :) = qext(idx_wave, :, :) 


  ! Calculation
  do ielem = 1, nelem
       
       call CARMAELEMENT_Get(carma, ielem, rc, igroup=igroup)
       if (rc /=0) stop "    *** CARMAELEMENT_Get FAILED ***"

       call CARMAGROUP_Get(carma, igroup, rc, ienconc=ienconc, cnsttype=cnsttype, r=r, rmass=rmass, maxbin=maxbin, &
               is_cloud=is_cloud, is_ice=is_ice, do_drydep=grp_do_drydep, rrat=rrat, arat=arat)
    
       do ibin = 1, nbin
           call CARMASTATE_GetBin(cstate, ielem, ibin, mmr, rc, &
                        numberDensity=numberDensity, nucleationRate=nucleationRate, r_wet=rwet, rhop_wet=rhop_wet, &
                        sedimentationflux=dd, vd=vd, vf=vf, dtpart=dtpart, totalmmr=totalmmr)
          
           if (rwet(1) == 0.0_f) rwet(:)=r(ibin)
           
           if (rc /=0) stop "    *** CARMASTATE_GetBin FAILED ***"
                   
           vmr_elm_bin(:, ielem, ibin) = mmr(:) 
           if (ielem == ienconc) then
                
                vnd(:, igroup)   = vnd(:, igroup)  + numberDensity(:)             
                re2(:, igroup)   = re2(:, igroup)  + numberDensity(:) * ((r(ibin)*rrat(ibin))**2)
                re3(:, igroup)   = re3(:, igroup)  + numberDensity(:) * ((r(ibin)*rrat(ibin))**3)
                rew2(:, igroup)  = rew2(:, igroup) + numberDensity(:) * ((rwet*rrat(ibin))**2)
                rew3(:, igroup)  = rew3(:, igroup) + numberDensity(:) * ((rwet*rrat(ibin))**3) 
                vad(:, igroup)   = vad(:, igroup)  + numberDensity(:) * 4.0_f * PI * (rwet**2) * 1.0e8_f
                vmd(:, igroup)   = vmd(:, igroup)  + numberDensity(:) * rmass(ibin)
                vmr(:, igroup)   = vmr(:, igroup)  + totalmmr(:)
                vpa(:, igroup)   = vpa(:, igroup)  + numberDensity(:) * PI * ((rwet * rrat(ibin))**2) * arat(ibin) 
                vvm(:, igroup)   = vvm(:, igroup)  + numberDensity(:) * rmass(ibin) * vf(2:) 
                vex(:, igroup)   = vex(:, igroup)  + numberDensity(:) * vqext(ibin,igroup) * PI * (r(ibin)**2) * 1e5_f 
                vod(:, igroup)   = vod(:, igroup)  + numberDensity(:) * vqext(ibin,igroup) * PI * (r(ibin)**2) * dz * 100._f
     
                vmr_bin(:, igroup, ibin) = totalmmr(:)
                vnd_bin(:, igroup, ibin) = numberDensity(:) 
                vro_bin(:, igroup, ibin) = rhop_wet(:)
                vwr_bin(:, igroup, ibin) = rwet(:)* 1e4_f
 
                vjn(:, igroup)  = vjn(:, igroup)  + nucleationRate(:)
  
                if (cnsttype == I_CNSTTYPE_PROGNOSTIC) then
                        if (ibin <= maxbin) then
                                if (grp_do_drydep) then
                                        vvd_bin(:, igroup, ibin) = - vd / 100._f
                                end if
                        end if
                end if

        end if
     end do   ! loop over bins

      if (ielem == ienconc) then
        where (re2(:, igroup) > 0.0_f)
                vre(:, igroup)  = (re3(:, igroup) / re2(:, igroup)) * 1e4_f 
                vrew(:, igroup) = (rew3(:, igroup) / rew2(:, igroup)) * 1e4_f
                vrm(:, igroup)  = (3._f / 4._f) * (vmd(:, igroup)  / (0.917_f * vpa(:, igroup))) * 1e4_f
                var(:, igroup)  = vpa(:, igroup) / PI / rew2(:, igroup)
        end where

        where (vmd(:, igroup) > 0.0_f)
                vvm(:, igroup) = vvm(:, igroup) / vmd(:, igroup)
        end where
     end if

  end do      ! loop over elements


  if (present(nd))     nd(:, :)  = vnd(:, :)
  if (present(ad))     ad(:, :)  = vad(:, :) 
  if (present(md))     md(:, :)  = vmd(:, :)
  if (present(re))    re(:, :)  = vre(:, :) 
  if (present(rew))    rew(:, :) = vrew(:, :) 
  if (present(rm))     rm(:, :)  = vrm(:, :)
  if (present(jn))     jn(:, :)  = vjn(:, :)
  if (present(mr))     mr(:, :)  = vmr(:, :)
  if (present(pa))     pa(:, :)  = vpa(:, :)
  if (present(ar))     ar(:, :)  = var(:, :)
  if (present(vm))     vm(:, :)  = vvm(:, :)
  if (present(ex))     ex(:, :)  = vex(:, :)
  if (present(od))     od(:, :)  = vod(:, :)
  if (present(wr_bin)) wr_bin(:, :, :) = vwr_bin(:, :, :)
  if (present(nd_bin)) nd_bin(:, :, :) = vnd_bin(:, :, :)
  if (present(ro_bin)) ro_bin(:, :, :) = vro_bin(:, :, :)
  if (present(mr_bin)) mr_bin(:, :, :) = vmr_bin(:, :, :)
  if (present(vd_bin)) vd_bin(:, :, :) = vvd_bin(:, :, :)
  if (present(mr_elm_bin))  mr_elm_bin(:, :, :)  = vmr_elm_bin(:, :, :)
end subroutine carmadiags

end module
