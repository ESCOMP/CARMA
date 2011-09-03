  openr, lun, 'carma_sulfatetest.txt', /get_lun

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
  mmrgas = reform(mmrgas,ngas,nt)
  satliq = reform(satliq,ngas, nt)
  satice = reform(satice,ngas, nt)
   
  nsubstep[0] = !values.f_nan
  nretry[0]   = !values.f_nan

  !p.multi = [0,1,5]
  loadct, 39

  ;Calculate the column mass, which should be conserved.
  mmrelem   = fltarr(nt,nelem)
  mmrtotal   = fltarr(nt)
  
  for ielem = 0, nelem-1 do begin
    for it = 0, nt-1 do begin
      mmrelem[it,ielem] = total(mmr[ielem,it,*])
      mmrtotal[it] = total(mmrelem[it,*]) + mmrgas[1, it]
    endfor
  endfor
  
 mmr[where(mmr le 0.)]           = !Values.F_NAN
;mmrelem[where(mmrelem le 0.)]   = !Values.F_NAN
;mmrtotal[where(mmrtotal le 0.)] = !Values.F_NAN
; satliq[where(satliq le 0.)]     = !Values.F_NAN
; satice[where(satice le 0.)]     = !Values.F_NAN
 satliq[*, 0] = !Values.F_NAN
 satice[*, 0] = !Values.F_NAN

 
  for it = 0, nt-1 do begin
    ; ======plot 1 ==============  
    plot, r[*], mmr[0,0,*], yrange=[1e-35, 10*max(mmrtotal)], $
         title = 'time = '+string(time[it])+' seconds', $
         xtitle='Radius [um]', ytitle = 'MMR [kg/kg]', thick=6, $
         /XLOG, /YLOG, charsize=2.0
 
    ; Add a legend
;    plots, [1.5e-10,1.75e-10], 1.0, thick=3, color=66, lin=0
;    plots, [1.5e-10,1.75e-10], 1.0, thick=3, color=66, lin=1
;    plots, [1.5e-10,1.75e-10], 1.2, thick=3, color=66, lin=2
;    xyouts, 1.8e-10, 4.2, 'sulfate mmr', color=66
;    xyouts, 1.8e-10, 2.7, 'gas mmr', color=96
;    xyouts, 1.8e-10, 1.2, 'Gerber', color=66


    oplot, r[*], mmr[0,it,*], thick=6, color=66                                  ; Sulfates, darkblue

    oplot, [min(r), max(r)], [mmrgas[1, it], mmrgas[1, it]], thick=6, color=96   ; H2SO4 gas, light blue

    
    ; ======plot 2 ============== 
    ; Show the mmr evolution.
    plot, mmrtotal[*], xtitle = 'Time Step', ytitle = 'mmr [kg/kg]', thick=6, $
        title = 'Total mmr evolution', charsize=2.0, $
        xrange=[0,nt-1], $
        yrange=[0., 1.5*max(mmrtotal[*])] 

    ; Add a legend
    plots, [nt/25.,  nt/25.+1.],  1.4*max(mmrtotal[*]), thick=3, lin=0, color=66
    plots, [5.*nt/50.,5.*nt/50.+1.],  1.4*max(mmrtotal[*]), thick=3, lin=0, color=96
    plots, [11.*nt/50.,11.*nt/50.+1.], 1.4*max(mmrtotal[*]), thick=3, lin=0, color=26    
    xyouts, nt/25.+1.5, 1.3*max(mmrtotal[*]), 'Sulfate', color=66
    xyouts, 5.*nt/50.+1.5, 1.3*max(mmrtotal[*]), 'Sulfate Gas', color=96
    xyouts, 11.*nt/50.+1.5, 1.3*max(mmrtotal[*]), 'Total H2SO4', color=26
    
    oplot, mmrelem[*,0], thick=6                     ; Sulfates, trajectory

     
    oplot, mmrgas[1,*], thick=6                      ; H2SO4 gas, trajectory    
   

    oplot, mmrtotal[0:it], thick=6, color=26         ; Total H2SO4, purple

    
    oplot, mmrelem[0:it,0], thick=6, color=66        ; Sulfates, dark blue
    
       
    oplot, mmrgas[1, 0:it], thick=6, color=96        ; H2SO4 gas, light blue
   
    ; ======plot 3 ==============
    ; Show the saturation evolution.
    plot, satliq[*, 2], xtitle = 'Time Step', ytitle = 's', thick=6, $
        title = 'Gas Saturation Ratio', $
        xrange=[0,nt-1], yrange=[0,max(satliq, /NAN)], charsize=2.0 
        
    oplot, satliq[1, *], thick=6
    oplot, satice[1, *], thick=6

    oplot, satliq[1, 0:it], thick=6, color=196     ; liquid, yellow

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

    wait, 15. / nt
  endfor

end
