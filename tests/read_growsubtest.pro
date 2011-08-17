  openr, lun, 'carma_growsubtest.txt', /get_lun

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
  nsubstep_ = 0
  nretry_ = 0

  nt = 0
  while(not(eof(lun))) do begin
    readf, lun, t1
    
    readf, lun, nsubstep_, nretry_, t_
    
    for ielem = 0, nelem-1 do begin
      for ibin = 0, nbin-1 do begin
        readf, lun, data
      
        mmr_[ielem, ibin]  = data[2]
      endfor
    endfor
   
    if (nt eq 0) then begin
      time = t1
      mmr  = mmr_
      nsubstep = nsubstep_
      nretry   = nretry_
      t    = t_
    endif else begin
      time = [time,t1]
      mmr  = [mmr,mmr_]
      nsubstep = [nsubstep,nsubstep_]
      nretry   = [nretry,nretry_]
      t        = [t,t_]
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
  
  nsubstep[0] = !values.f_nan
  nretry[0]   = !values.f_nan
  

  !p.multi = [0,2,3]
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
    plot, r[*], mmr[0,0,*], yrange=[1e-30, 10*max(mmrtotal)], $
         title = 'time = '+string(time[it])+' seconds', $
         xtitle='Radius [um]', ytitle = 'MMR [kg/kg]', thick=6, $
         /XLOG, /YLOG, charsize=2.0
 
    ; Add a legend
    plots, [60,62], 1e-10, thick=3, lin=0, color=66
    plots, [60,62], 1e-5, thick=3, lin=0, color=96
    xyouts, 63, 1e-10, 'Ice', color=66
    xyouts, 63, 1e-5, 'Water Vapor', color=96

    for ielem = 0, nelem-1 do begin
      oplot, r[*], mmr[ielem,it,*], lin=ielem, thick=6, color=66
    endfor

   for igas = 0, ngas-1 do begin
      oplot, [min(r), max(r)], [mmrgas[it, igas], mmrgas[it, igas]], thick=6, color=96, lin=igas
    endfor

    ; Show the mmr evolution.
    plot, mmrtotal[*], xtitle = 'Time Step', ytitle = 'mmr [kg/kg]', thick=6, $
        title = 'Total mmr evolution', charsize=2.0, $
        yrange=[min([min(mmrtotal), min(mmrgas), min(mmrelem)]), max([max(mmrtotal), 1.5*max(mmrgas), max(mmrelem)])] 

    ; Add a legend
    plots, [.25,.5], 5.5e-6, thick=3, lin=0, color=66
    plots, [1.5,1.75], 5.5e-6, thick=3, lin=0, color=96
    plots, [3.25,3.5], 5.5e-6, thick=3, lin=0, color=26
    xyouts, .75, 5.25e-6, 'Ice', color=66
    xyouts, 2., 5.25e-6, 'Water Vapor', color=96
    xyouts, 3.75, 5.25e-6, 'Total Water', color=26


    for ielem = 0, nelem-1 do begin
      oplot, mmrelem[*,ielem], thick=6, lin=ielem
    endfor
     for igas = 0, ngas-1 do begin
      oplot, mmrgas[*,igas], thick=6, lin=igas
    endfor

    oplot, mmrtotal[0:it], thick=6, color=26

    for ielem = 0, nelem-1 do begin
      oplot, mmrelem[0:it,ielem], thick=6, color=66, lin=ielem
    endfor
    
    for igas = 0, ngas-1 do begin
      oplot, mmrgas[0:it, igas], thick=6, color=96, lin=igas
    endfor
    

    ; Show the saturation evolution.
    plot, satice[*], xtitle = 'Time Step', ytitle = 's', thick=6, $
        title = 'Gas Saturation Ratio', $
        yrange=[0, 2], charsize=2.0 
        
    ; Add a legend
    plots, [1,1.25], 1.75, thick=3, lin=0, color=66
    plots, [3.5,3.75], 1.75, thick=3, lin=0, color=196
    xyouts, 1.5, 1.7, 'Sat Ice', color=66
    xyouts, 4, 1.7, 'Sat Liq', color=196

    oplot, [0, nt], [1., 1.], thick=3

    for igas = 0, ngas-1 do begin
      oplot, satliq[*,igas], thick=6, lin=igas
      oplot, satice[*,igas], thick=6, lin=igas
    endfor

    for igas = 0, ngas-1 do begin
      oplot, satliq[0:it, igas], thick=6, color=196, lin=igas
      oplot, satice[0:it, igas], thick=6, color=66, lin=igas
    endfor
    

    ; Show the temperature evolution.
    plot, t[*], xtitle = 'Time Step', ytitle = 'dT (K)', thick=6, $
        title = 'Delta Temperature', $
        yrange=[0., max(t)], charsize=2.0 
        
    oplot, t[0:it], thick=6, lin=0, color=66

    

    ; Show the substepping evolution.
    plot, nsubstep[0:*], xtitle = 'Time Step', ytitle = 'Nsubsteps', thick=6, $
        title = 'Number of Substeps', $
        yrange=[0., 1.2*max(nsubstep)], charsize=2.0
        
    oplot, nsubstep[0:it], thick=6, lin=0, color=66

    ; Show the retry evolution.
    plot, nretry[0:*], xtitle = 'Time Step', ytitle = 'Nretry', thick=6, $
        title = 'Number of Retries', $
        charsize=2.0 
        
    oplot, nretry[0:it], thick=6, lin=0, color=66

    wait, 10. / nt
  endfor
  
end
