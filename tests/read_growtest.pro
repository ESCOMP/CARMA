  openr, lun, 'carma_growtest.txt', /get_lun

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
  
  t_ = 0.
  rlheat_ = 0.

  nt = 0L
  while(not(eof(lun))) do begin
    readf, lun, t1

    readf, lun, t_, rlheat_

    for ielem = 0, nelem-1 do begin
      for ibin = 0, nbin-1 do begin
        readf, lun, data
      
        mmr_[ielem, ibin]  = data[2]
      endfor
    endfor
   
    if (nt eq 0) then begin
      time   = t1
      mmr    = mmr_
      t      = t_
      rlheat = rlheat_
    endif else begin
      time   = [time,t1]
      mmr    = [mmr,mmr_]
      t      = [t,t_]
      rlheat = [rlheat,rlheat_]
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

  !p.multi = [0,1,5]
  loadct, 39

  ;Calculate the column mass, which should be conserved.
  mmrelem   = fltarr(nt,nelem)
  mmrtotal   = fltarr(nt)
  
  for ielem = 0, nelem-1 do begin
    for it = 0L, nt-1 do begin
      mmrelem[it,ielem] = total(mmr[ielem,it,*])
      mmrtotal[it] = total(mmrelem[it,*]) + total(mmrgas[it, *])
    endfor
  endfor
  
 mmr[where(mmr le 0.)]           = !Values.F_NAN
; mmrelem[where(mmrelem le 0.)]   = !Values.F_NAN
; mmrtotal[where(mmrtotal le 0.)] = !Values.F_NAN
satliq[0,*] = !Values.F_NAN
satice[0,*] = !Values.F_NAN
 
  for it = 0L, nt-1 do begin
    plot, r[*], mmr[0,0,*], yrange=[1e-30, 10*max(mmrtotal)], $
         title = 'time = '+string(time[it])+' seconds', $
         xtitle='Radius [um]', ytitle = 'MMR [kg/kg]', thick=3, $
         /XLOG, /YLOG, charsize=2.0
 
    ; Add a legend
    plots, [60,62], 1e-10, thick=3, lin=0, color=66
    plots, [60,62], 1e-5, thick=3, lin=0, color=96
    xyouts, 63, 1e-10, 'Ice', color=66
    xyouts, 63, 1e-5, 'Water Vapor', color=96

    for ielem = 0, nelem-1 do begin
      oplot, r[*], mmr[ielem,it,*], lin=ielem, thick=3, color=66
    endfor

   for igas = 0, ngas-1 do begin
      oplot, [min(r), max(r)], [mmrgas[it, igas], mmrgas[it, igas]], thick=3, color=96, lin=igas
    endfor

    ; Show the mmr evolution.
    plot, mmrtotal[*], xtitle = 'Time Step', ytitle = 'mmr [kg/kg]', thick=3, $
        title = 'Total mmr evolution', charsize=2.0, $
        yrange=[min([min(mmrtotal), min(mmrgas), min(mmrelem)]), max([max(mmrtotal), 1.5*max(mmrgas), max(mmrelem)])] 

    ; Add a legend
    plots, [2,2.5], 5.5e-6, thick=3, lin=0, color=66
    plots, [4,4.5], 5.5e-6, thick=3, lin=0, color=96
    plots, [9,9.5], 5.5e-6, thick=3, lin=0, color=26
    xyouts, 2.7, 5.25e-6, 'Ice', color=66
    xyouts, 4.7, 5.25e-6, 'Water Vapor', color=96
    xyouts, 9.7, 5.25e-6, 'Total Water', color=26


    for ielem = 0, nelem-1 do begin
      oplot, mmrelem[*,ielem], thick=3, lin=ielem
    endfor
     for igas = 0, ngas-1 do begin
      oplot, mmrgas[*,igas], thick=3, lin=igas
    endfor

    oplot, mmrtotal[0:it], thick=3, color=26

    for ielem = 0, nelem-1 do begin
      oplot, mmrelem[0:it,ielem], thick=3, color=66, lin=ielem
    endfor
    
   for igas = 0, ngas-1 do begin
      oplot, mmrgas[0:it, igas], thick=3, color=96, lin=igas
    endfor

    ; Show the saturation evolution.
    plot, satice[*], xtitle = 'Time Step', ytitle = 's', thick=3, $
        title = 'Gas Saturation Ratio', $
        yrange=[0, 5], charsize=2.0 
        
    ; Add a legend
    plots, [2,2.5], 4.5, thick=3, lin=0, color=66
    plots, [2,2.5], 3.5, thick=3, lin=0, color=196
    xyouts, 2.7, 4.25, 'Sat Ice', color=66
    xyouts, 2.7, 3.25, 'Sat Liq', color=196

    oplot, [0, nt], [1., 1.], thick=3

    for igas = 0, ngas-1 do begin
      oplot, satliq[*,igas], thick=3, lin=igas
      oplot, satice[*,igas], thick=3, lin=igas
    endfor

   for igas = 0, ngas-1 do begin
      oplot, satliq[0:it, igas], thick=3, color=196, lin=igas
      oplot, satice[0:it, igas], thick=3, color=66, lin=igas
    endfor

    ; Show the temperature evolution.
    plot, t[*], xtitle = 'Time Step', ytitle = 'dT (K)', thick=3, $
        title = 'Delta Temperature', $
        yrange=[0., max(t)], charsize=2.0 
        
    oplot, t[0:it], thick=3, lin=0, color=66

    ; Show the latent heat.
    plot, rlheat[*], xtitle = 'Time Step', ytitle = 'LH (K/s)', thick=3, $
        title = 'Latent Heat', $
        yrange=[min(rlheat), max(rlheat)], charsize=2.0 
        
    oplot, rlheat[0:it], thick=3, lin=0, color=66

    wait, 2. / nt
  endfor
  
end
