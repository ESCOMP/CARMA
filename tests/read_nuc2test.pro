  openr, lun, 'carma_nuc2test.txt', /get_lun

  ; Read in the sizes.
  readf, lun, ngroup, nelem, nbin, ngas
  
  r     = fltarr(ngroup, nbin)
  rmass = fltarr(ngroup, nbin)

  data = fltarr(4)
  
  for igroup = 0, ngroup-1 do begin
    for ibin = 0, nbin-1 do begin
      readf, lun, data
    
      r[igroup, ibin]     = data[2]
      rmass[igroup, ibin] = data[3]
    endfor
  endfor

  ; Read in the particles for each time step.
  mmr_     = fltarr(nelem, nbin)
  mmrgas_  = fltarr(ngas)
  satliq_  = fltarr(ngas)
  satice_  = fltarr(ngas)

  data = fltarr(3)
  datag = fltarr(4)

  nt = 0
  while(not(eof(lun))) do begin
    readf, lun, t1
    for ielem = 0, nelem-1 do begin
      for ibin = 0, nbin-1 do begin
        readf, lun, data
      
        mmr_[ielem, ibin]  = data[2]
      endfor
    endfor
   
    if (nt eq 0) then begin
      time = t1
      mmr  = mmr_
    endif else begin
      time = [time,t1]
      mmr  = [mmr,mmr_]
    endelse
   
    for igas = 0, ngas-1 do begin
      readf, lun, datag
      
      mmrgas_[igas]  = datag[1]
      satliq_[igas]  = datag[2]
      satice_[igas]  = datag[3]
    endfor
   
    if (nt eq 0) then begin
      mmrgas  = mmrgas_
      satliq  = satliq_
      satice  = satice_
    endif else begin
      mmrgas  = [mmrgas,mmrgas_]
      satliq  = [satliq,satliq_]
      satice  = [satice,satice_]
    endelse

    nt = nt+1
  endwhile
  
  free_lun, lun

  mmr    = reform(mmr,nelem,nt,nbin)
  mmrgas = reform(mmrgas,nt,ngas)
  satliq = reform(satliq,nt,ngas)
  satice = reform(satice,nt,ngas)
  
  
  !p.multi = [0,1,4]
  loadct, 39

  ;Calculate the column mass, which should be conserved.
  mmrelem   = fltarr(nt,nelem)
  mmrtotal   = fltarr(nt)
  
  for ielem = 0, nelem-1 do begin
    for it = 0, nt-1 do begin
      mmrelem[it,ielem] = total(mmr[ielem,it,*])
      mmrtotal[it] = total(mmrelem[it,*]) + total(mmrgas[it, *])
    endfor
  endfor
  
 mmr[where(mmr le 0.)]           = !Values.F_NAN
; mmrelem[where(mmrelem le 0.)]   = !Values.F_NAN
; mmrtotal[where(mmrtotal le 0.)] = !Values.F_NAN
satliq[0,*] = !Values.F_NAN
satice[0,*] = !Values.F_NAN
 
  for it = 0, nt-1 do begin
    plot, r[0,*]*1000., mmr[0,0,*], yrange=[1e-20, 1e-10], xrange=[min(r[0,*])*1000., max(r[0,*])*1000.], $
         title = 'Sulfate, time = '+string(time[it])+' seconds', $
         xtitle='Radius [nm]', ytitle = 'MMR [kg/kg]', thick=4, $
         /XLOG, /YLOG, charsize=2.0
 
    ; Add a legend
;    plots, [1.5e-10,1.75e-10], 4.2, thick=3, color=66, lin=0
;    plots, [1.5e-10,1.75e-10], 2.7, thick=3, color=66, lin=1
;    plots, [1.5e-10,1.75e-10], 1.2, thick=3, color=66, lin=2
;    xyouts, 1.8e-10, 4.2, 'None', color=66
;    xyouts, 1.8e-10, 2.7, 'Fitzgerald', color=66
;    xyouts, 1.8e-10, 1.2, 'Gerber', color=66

    oplot, r[0,*]*1000., mmr[0,it,*], thick=4, color=66

    plot, r[1,*], mmr[1,0,*], yrange=[1e-15, 1e-5], xrange=[min(r[1,*]), max(r[1,*])], $
         title = 'H2O and Ice', $
         xtitle='Radius [um]', ytitle = 'MMR [kg/kg]', thick=4, $
         /XLOG, /YLOG, charsize=2.0

    for ielem = 1, nelem-1 do begin
      oplot, r[1,*], mmr[ielem,it,*], lin=ielem, thick=4, color=66
    endfor

    for igas = 0, ngas-1 do begin
      oplot, [min(r[1,*]), max(r[1,*])], [mmrgas[it, igas], mmrgas[it, igas]], thick=4, color=96, lin=igas
    endfor

    ; Show the mmr evolution.
    plot, mmrtotal[*], xtitle = 'Time Step', ytitle = 'mmr [kg/kg]', thick=4, $
        title = 'Total mmr evolution', charsize=2.0, $
        yrange=[min([min(mmrtotal), min(mmrgas), min(mmrelem)]), max([max(mmrtotal), 1.5*max(mmrgas), max(mmrelem)])] 

    for ielem = 0, nelem-1 do begin
      oplot, mmrelem[*,ielem], thick=4, lin=ielem
    endfor
    
    for igas = 0, ngas-1 do begin
      oplot, mmrgas[*,igas], thick=4, lin=igas
    endfor

    oplot, mmrtotal[0:it], thick=4, color=26

    for ielem = 0, nelem-1 do begin
      oplot, mmrelem[0:it,ielem], thick=4, color=66, lin=ielem
    endfor
    
   for igas = 0, ngas-1 do begin
      oplot, mmrgas[0:it, igas], thick=4, color=96, lin=igas
    endfor

    ; Show the saturation evolution.
    plot, satice[*], xtitle = 'Time Step', ytitle = 's', thick=4, $
        title = 'Gas Saturation Ratio', $
        yrange=[0, 2], charsize=2.0 
        
    oplot, [0, nt], [1., 1.], thick=2

    for igas = 0, ngas-1 do begin
      oplot, satliq[*,igas], thick=4, lin=igas
      oplot, satice[*,igas], thick=4, lin=igas
    endfor

   for igas = 0, ngas-1 do begin
      oplot, satliq[0:it, igas], thick=4, color=196, lin=igas
      oplot, satice[0:it, igas], thick=4, color=66, lin=igas
    endfor

    wait, 20. / nt
  endfor
  
;  plot, q[*,0,0], z[*], yrange=[0,15], xrange=[0,2.e-10], $
;       title = 'Falling History', $
;       xtitle='Particle Concentration [g cm-3]', ytitle = 'Altitude [km]', thick=3
;    oplot, [0,2.5e-10], [8,8], thick=1, lin=1
;  xyouts, 1.5e-10, 9.2, 't = 0 sec'

  ; Add a legend
;  plots, [1.5e-10,1.75e-10], 4.2, thick=3, lin=0
;  plots, [1.5e-10,1.75e-10], 2.7, thick=3, lin=1
;  plots, [1.5e-10,1.75e-10], 1.2, thick=3, lin=2
;  xyouts, 1.8e-10, 4.2, 'None'
;  xyouts, 1.8e-10, 2.7, 'Fitzgerald'
;  xyouts, 1.8e-10, 1.2, 'Gerber'

;  it = 200
;  for ielem = 0, nelem-1 do begin
;    oplot, q[*,it,ielem], z[*], color=86, thick=3, line=ielem
;  endfor
;  oplot, [0,2.5e-10], [4,4], color=86, thick=1, lin=1
;  xyouts, 1.5e-10, 7.7, 't = 200000 sec', color=86
  
;  it = 400
;  for ielem = 0, nelem-1 do begin
;    oplot, q[*,it,ielem], z[*], color=126, thick=3, line=ielem
;  endfor
;  oplot, [0,2.5e-10], [0,0], color=126, thick=1, lin=1
;  xyouts, 1.5e-10, 6.2, 't = 400000 sec', color=126
  
;  plot, mass[0,*], xtitle = 'Time Step', ytitle = 'Column Mass [g cm-2]', thick=6, $
;        title = 'Total mass evolution'
;  for ielem = 0, nelem-1 do begin
;    oplot, mass[ielem,*], thick=6, color=66, lin=ielem
;  endfor

end
