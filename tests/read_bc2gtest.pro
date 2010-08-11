  openr, lun, 'carma_bc2gtest.txt', /get_lun
  readf, lun, nbin, nelem, ngroup
  
;  ibin_ = intarr(nbin)
  r    = fltarr(nbin,ngroup)      ; radius (cm)
  dr   = fltarr(nbin,ngroup)      ; delta radius (cm)
  data = fltarr(ngroup+1)
  
  for igroup = 0, ngroup-1 do begin
    for ibin = 0, nbin-1 do begin
      readf, lun, data
;     ibin_[ibin] = fix( data[0] )
      r[ibin,igroup]   = data[1]
      dr[ibin,igroup]   = data[2]
    endfor
  endfor

  data2 = fltarr(2)
  pc_   = fltarr(nbin,nelem)      ; mass   (g/cm^-3) in a bin
  
  
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
;  3 elements: 1 - OC, 2 - BC shell, 3 - OC core
;  Elements 1 - 2 are PC; element 3 is COREMASS
;  For now, all elements have same size and mass dimensions
  d  = 2.*r*1e4
  dd = 2.*dr*1e4
  
  dnlogd = alog(10.) * pc
  for ibin = 0, nbin-1 do begin
    dnlogd[ibin,*,0] = dnlogd[ibin,*,0]*d[ibin,0]/dd[ibin,0]
    dnlogd[ibin,*,1] = dnlogd[ibin,*,1]*d[ibin,0]/dd[ibin,0]
  endfor
  
  ;  BCOC
  massbc   = total(pc[*,*,1],1) - total(pc[*,*,2],1) 
  massoc   = total(pc[*,*,0],1) + total(pc[*,*,2],1) 
  mass     = massbc + massoc
  
  massmbc  = total(pc[*,*,1],1) - total(pc[*,*,2],1) 
  massmoc  = total(pc[*,*,2],1) 
  massmix  = massmbc + massmoc


  ; Plot
  !p.font=0
  !p.multi = [0,1,2]
  loadct, 39
  
  ; Can't display 0 in log coordinates
  dnlogd[where(dnlogd le 0.)] = !Values.F_NAN
  
  for it = 0, nt-1 do begin
    plot, d[*,0], dnlogd[*,0,0], $
     /xlog, /ylog, /nodata, xtitle = 'Diameter [um]', ytitle = 'dM/dlogD [g/cm3]', $ 
     xrange=[.005,.3], yrange=[1.e-21,1.e-11], xstyle=1
     
    oplot, d[*,0], dnlogd[*,0,0], color=254, psym=2
    oplot, d[*,1], dnlogd[*,0,1], color=176, psym=4

    oplot, d[*,1], dnlogd[*,it,0], thick=6, color=254, lin=0
    oplot, d[*,0], dnlogd[*,it,1], thick=3, color=176, lin=0
    oplot, d[*,0], dnlogd[*,it,2], thick=3, color=84, lin=0

    xyouts, .1, 1.e-12, 'OC (pure group)', color=254
    xyouts, .1, 1.e-13, 'OC (core of mixed group)', color=84
    xyouts, .1, 1.e-14, 'BC (shell of mixed group)', color=176

    plot, [0.,1.], [mass[it], mass[it]], $
     /nodata, ytitle = 'Total Mass Density [g/cm3]', $
     yrange=[0.,2.5e-13], xstyle=1
    oplot, [0.,1.], [mass[it], mass[it]], thick=3, color=66, lin=0
    oplot, [0.,1.], [massoc[it], massoc[it]], thick=9, color=84, lin=0
    oplot, [0.,1.], [massbc[it], massbc[it]], thick=3, color=176, lin=0
    oplot, [0.,1.], [massmix[it], massmix[it]], thick=3, color=66, lin=2
    oplot, [0.,1.], [massmoc[it], massmoc[it]], thick=3, color=84, lin=2
    oplot, [0.,1.], [massmbc[it], massmbc[it]], thick=6, color=176, lin=2

        fac = 1./80.
        plots, fac*[50,55], 2.0e-13, thick=3, color=66
        plots, fac*[50,55], 1.8e-13, thick=3, color=66, lin=2
        xyouts, fac*57, 1.95e-13, 'Total Mass', color=66
        xyouts, fac*57, 1.75e-13, 'Mixed Group Mass', color=66

        plots, fac*[50,55], 0.9e-13, thick=3, color=176
        plots, fac*[50,55], 0.7e-13, thick=9, color=84
        plots, fac*[50,55], 0.5e-13, thick=6, color=176, lin=2
        plots, fac*[50,55], 0.3e-13, thick=3, color=84, lin=2
        xyouts, fac*57, 0.85e-13, 'Total BC', color=176
        xyouts, fac*57, 0.65e-13, 'Total OC', color=84
        xyouts, fac*57, 0.45e-13, 'Mixed BC', color=176
        xyouts, fac*57, 0.25e-13, 'Mixed OC', color=84
    
    wait, .1
  endfor
  
  ; At the end, show the mass evolution.
	plot, d[*,0], dnlogd[*,0,0], $
	 /xlog, /ylog, /nodata, xtitle = 'Diameter [um]', ytitle = 'dM/dlogD [g/cm3]', $ 
	 xrange=[.005,.3], yrange=[1.e-21,1.e-11], xstyle=1
	 
	oplot, d[*,0], dnlogd[*,0,0], color=254, psym=2
	oplot, d[*,1], dnlogd[*,0,1], color=176, psym=4

	oplot, d[*,1], dnlogd[*,it-1,0], thick=6, color=254, lin=0
	oplot, d[*,0], dnlogd[*,it-1,1], thick=3, color=176, lin=0
	oplot, d[*,0], dnlogd[*,it-1,2], thick=3, color=84, lin=0

        xyouts, .1, 1.e-12, 'OC (pure group)', color=254
        xyouts, .1, 1.e-13, 'OC (core of mixed group)', color=84
        xyouts, .1, 1.e-14, 'BC (shell of mixed group)', color=176
    

	plot, mass, xtitle = 'Time Step', ytitle = 'Total Mass Density [g/cm3]'
	oplot, mass, thick=3, color=66, lin=0
	oplot, massoc, thick=9, color=84, lin=0
	oplot, massbc, thick=3, color=176, lin=0
	oplot, massmix, thick=3, color=66, lin=2
	oplot, massmoc, thick=3, color=84, lin=2
	oplot, massmbc, thick=6, color=176, lin=2

        plots, [50,55], 2.0e-13, thick=3, color=66
        plots, [50,55], 1.8e-13, thick=3, color=66, lin=2
        xyouts, 57, 1.95e-13, 'Total Mass', color=66
        xyouts, 57, 1.75e-13, 'Mixed Group Mass', color=66

        plots, [50,55], 0.9e-13, thick=3, color=176
        plots, [50,55], 0.7e-13, thick=9, color=84
        plots, [50,55], 0.5e-13, thick=6, color=176, lin=2
        plots, [50,55], 0.3e-13, thick=3, color=84, lin=2
        xyouts, 57, 0.85e-13, 'Total BC', color=176
        xyouts, 57, 0.65e-13, 'Total OC', color=84
        xyouts, 57, 0.45e-13, 'Mixed BC', color=176
        xyouts, 57, 0.25e-13, 'Mixed OC', color=84


end
