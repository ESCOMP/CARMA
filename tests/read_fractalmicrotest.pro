  openr, lun, 'carma_fractalmicrotest.txt', /get_lun

  ; Read in bin structure
  readf,lun, nbin, rmon
  
  df = fltarr(nbin)
  r = fltarr(nbin)
  rrat = fltarr(nbin)
  rprat = fltarr(nbin)  
  nmon = fltarr(nbin)

  data= fltarr(6)
  for ib = 0, nbin-1 do begin
    readf, lun, data
    df[ib] = data[1]
    r[ib] = data[2]
    rrat[ib] = data[3]
    rprat[ib] = data[4]
    nmon[ib] = data[5]
  endfor
  
  ; Write bin parameters out
  ;print, "ib,  df,  r_sphere,   rrat,    rprat"
  ;for ib = 0, nbin-1 do begin
  ;  print, ib, df(ib), r(ib), rrat(ib), rprat(ib)
  ;endfor


  ; Read in the vertical grid.
  readf, lun, nz

  z    = fltarr(nz)
  dz   = fltarr(nz)

  data = fltarr(3)
  
  for iz = 0, nz-1 do begin
    readf, lun, data
    
    z[iz]  = data[1]
    dz[iz] = data[2]
  endfor

  ; Read in the particles for each time step.
  pc_f_  = fltarr(nz)   ; particle concentration #/cm-3
  q_f_    = fltarr(nz)  ; mass concentration g/cm-3
  pc_s_  = fltarr(nz)   ; particle concentration #/cm-3
  q_s_    = fltarr(nz)  ; mass concentration g/cm-3
  data2 = fltarr(5)

  nt = 0
  while(not(eof(lun))) do begin
    readf, lun, t1
    for iz = 0, nz-1 do begin
      ;for ib = 0, nbin-1 do begin
        readf, lun, data2
        pc_f_[iz]  = data2[1] 
        q_f_[iz]    = data2[2]
        pc_s_[iz]  = data2[3] 
        q_s_[iz]    = data2[4]
     ;endfor
   endfor 
  
    if (nt eq 0) then begin
      time = t1
      pc_f = pc_f_
      q_f    = q_f_
      pc_s = pc_s_
      q_s    = q_s_
    endif else begin
      time = [time,t1]
      pc_f  = [pc_f,pc_f_]
      q_f    = [q_f,q_f_]
      pc_s  = [pc_s,pc_s_]
      q_s    = [q_s,q_s_]
    endelse
    nt = nt+1
  endwhile
  
  free_lun, lun
  

  pc_f = reform(pc_f,nz,nt)
  q_f = reform(q_f,nz,nt)
  pc_s = reform(pc_s,nz,nt)
  q_s = reform(q_s,nz,nt)
  z = z/1000.

  ;Calculate the column mass, which should be conserved.
  zmass_f    = fltarr(nz,nt)
  pctot_f    = fltarr(nz,nt)
  zmass_s    = fltarr(nz,nt)
  pctot_s    = fltarr(nz,nt)

  for it = 0, nt-1 do begin
    for iz = 0, nz-1 do begin
      zmass_f[iz,it] = total(q_f[iz,it])
      pctot_f[iz,it] = total(pc_f[iz,it])
      zmass_s[iz,it] = total(q_s[iz,it])
      pctot_s[iz,it] = total(pc_s[iz,it])
    endfor
  endfor



  ;plotting to do
  ;plot effective r,rf,rm versus altitude
  ; try to duplicate plot RE_RB_RF

  r_eff_f = fltarr(nz,nt)
  r_eff_s = fltarr(nz,nt)
  rf_eff_f = fltarr(nz,nt)
  rm_eff_f = fltarr(nz,nt)

;  for it=0,nt-1 do begin
;    for iz=0,nz-1 do begin
;      rtemp3_f = 0. & rftemp3_f = 0. & rmtemp3_f = 0.
;      rtemp2_f = 0. & rftemp2_f = 0. & rmtemp2_f = 0.
;      rtemp3_s = 0.
;      rtemp2_s = 0.
;      for ib=0,nbin-1 do begin
;        rtemp3_f = rtemp3_f + 4.0/3.0*!pi*r(ib)^3.*pc_f(iz,it,ib) * 1.0e12
;        rtemp2_f = rtemp2_f + 4.0*!pi*r(ib)^2.*pc_f(iz,it,ib) * 1.0e8
;
;        rftemp3_f = rftemp3_f + 4.0/3.0*!pi*(r(ib)*rrat(ib))^3.*pc_f(iz,it,ib) * 1.0e12
;        rftemp2_f = rftemp2_f + 4.0*!pi*(r(ib)*rrat(ib))^2.*pc_f(iz,it,ib) * 1.0e8
;
;        rmtemp3_f = rmtemp3_f + 4.0/3.0*!pi*(r(ib)*rprat(ib))^3.*pc_f(iz,it,ib) * 1.0e12
;        rmtemp2_f = rmtemp2_f + 4.0*!pi*(r(ib)*rprat(ib))^2.*pc_f(iz,it,ib) * 1.0e8
;
;        rtemp3_s = rtemp3_s + 4.0/3.0*!pi*r(ib)^3.*pc_s(iz,it,ib) * 1.0e12
;        rtemp2_s = rtemp2_s + 4.0*!pi*r(ib)^2.*pc_s(iz,it,ib) * 1.0e8
;
;      endfor
;    
;      r_eff_f(iz,it) = 3.0*rtemp3_f/rtemp2_f
;      rf_eff_f(iz,it) = 3.0*rftemp3_f/rftemp2_f
;      rm_eff_f(iz,it) = 3.0*rmtemp3_f/rmtemp2_f
;      r_eff_s(iz,it) = 3.0*rtemp3_s/rtemp2_s
;    endfor
;  endfor

  !p.multi = [0,1,1]
  loadct, 39

  for it = 0, nt-1 do begin
;    plot, r_eff_f[*,nt-1], z[*], yrange=[0,110], xrange=[1.0e-3,1.0e3],xstyle=1,/xlog,/nodata, $
;         title = 'time = '+string(time[it])+' seconds', charsize=1.2,$
;         xtitle='Effective Radius (um)', ytitle = 'Altitude [km]', thick=3, $
;         xtickname=['10!U-3!N','10!U-2!N','10!U-1!N','10!U0!N','10!U1!N','10!U2!N','10!U3!N']
;    oplot, r_eff_s[*,it], z[*], lin=0, thick=3
;    oplot, r_eff_f[*,it], z[*], lin=2, thick=3, color=225
;    oplot, rf_eff_f[*,it], z[*], lin=3, thick=3, color=125
;    ;oplot, rm_eff_f[*,it], z[*], lin=1, thick=3
;    



 
    ;plot, pctot_f[*,0], z[*], yrange=[0,110], xrange=[1.0e-3,1.0e3],/xlog,xstyle=1,/nodata, $
    ;     title = 'time = '+string(time[it])+' seconds', charsize=1.5, $
    ;     xtitle='Particle Concentration [# cm-3]', ytitle = 'Altitude [km]', thick=3
    ;oplot, pctot_f[*,it], z[*], lin=2, thick=3
    ;oplot, pctot_s[*,it], z[*], lin=0, thick=3
 

    plot, zmass_f[*,0], z[*], yrange=[0,80], xrange=[0,2.0e-10],xstyle=1, /nodata, $
         title = 'time = '+string(time[it])+' seconds', charsize=1.2,$
         xtitle='Particle Concentration [g cm-3]', ytitle = 'Altitude [km]', thick=1

    oplot, zmass_s[*,it], z[*], lin=0, thick=3
    oplot, zmass_f[*,it], z[*], lin=0, thick=3, color=225

    xyouts, 0.15,0.89, "Fractal", color=225, /NORMAL, charsize=1.5
    xyouts, 0.15,0.84, "Spherical", /NORMAL, charsize=1.5

    wait, 0.01
  endfor
  

end


