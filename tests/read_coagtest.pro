  openr, lun, 'carma_coagtest.txt', /get_lun
  readf, lun, nbin, nelem, ngroup
  
;  ibin_ = intarr(nbin)
  r    = fltarr(nbin,ngroup)      ; radius (cm)
  dr   = fltarr(nbin,ngroup)      ; delta radius (cm)
  data = fltarr(3)
  
  for igroup = 0, ngroup-1 do begin
    for ibin = 0, nbin-1 do begin
      readf, lun, data
;     ibin_[ibin] = fix( data[0] )
      r[ibin,igroup]   = data[1]
      dr[ibin,igroup]   = data[2]
    endfor
  endfor

  data2 = fltarr(2)
  pc_   = fltarr(nbin,nelem)      ; number   (#/cm^-3) in a bin
  
  
  nt = 0
  while(not(eof(lun))) do begin
    readf, lun, t1
    
    for ielem = 0, nelem-1 do begin
      for ibin = 0, nbin-1 do begin
        readf, lun, data2
;    ibin_[ibin] = fix( data2[0] )
        pc_[ibin,ielem]    = data2[1]
      endfor
    endfor
    
;   ibin = [ibin,ibin_]
    if (t1 eq 0.) then begin
      pc = pc_
      time = t1
    endif else begin
      pc = [pc,pc_]
      time = [time,t1]
    endelse
   
    nt = nt+1
  endwhile
  
  free_lun, lun

  pc = reform(pc,nbin,nt,nelem)
  

; Map units for this test
;  4 elements: 1 - BC, 2 - OC, 3 - OC shell, 4 - BC core
;  Elements 1 - 3 are PC; element 4 is COREMASS
;  For now, all elements have same size and mass dimensions
  d  = 2.*r*1e4
  dd = 2.*dr*1e4
  
  dnlogd = alog(10.) * pc
  
  for ibin = 0, nbin-1 do begin
    dnlogd[ibin,*,0] = dnlogd[ibin,*,0]*d[ibin,0]/dd[ibin,0]
  endfor
  
  ; Can't display 0 in log coordinates
  dnlogd[where(dnlogd le 0.)] = !Values.F_NAN
  
  

; The data has been manipulated and converted to match Figure 2 from
; Jacobson et al., Atmospheric Environment 28, 1327-1338, 1994
; Initial monodisperse distribution is slightly different -- probably because
; I don't know what bin edges were used to calculate dr
;
; - JAS, August 15, 2007

  plot, [ 0.005, 0.1 ], [ 1., 1.e8 ], /nodata $
      , xtitle = 'Diameter (um)' $
      , /xlog, xstyle = 1 $
      , ytitle = 'dN / d log D (cm-3)' $
      , /ylog, ystyle = 1 $
      , charsize = 1.25
  oplot, d, dnlogd[*,0,0], psym=2, symsize=1.25 
  oplot, d, dnlogd[*,nt-1,0], lin=2

end
