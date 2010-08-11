  openr, lun, 'carma_mietest.txt', /get_lun


  ; Read in the wavelengths.
  readf, lun, NGROUP, NWAVE, NBIN

  wave   = fltarr(NWAVE)

  data = fltarr(2)
  
  for iwave = 0, NWAVE-1 do begin
    readf, lun, data
    
    wave[iwave]  = data[1]
  endfor


  ; Read in the radius, refractive index and optical properties.
  r      = fltarr(NBIN,NGROUP)
  refidx = fltarr(NWAVE,NGROUP)
  
  qext   = fltarr(NWAVE,NGROUP,NBIN)
  ssa    = fltarr(NWAVE,NGROUP,NBIN)
  asym   = fltarr(NWAVE,NGROUP,NBIN)

  for igroup = 0, NGROUP-1 do begin
    data = fltarr(2)
    
    for ibin = 0, NBIN-1 do begin
      readf, lun, data
      r(ibin, igroup) = data[1]
    endfor
  
    data = fltarr(3)
    
    for iwave = 0, NWAVE-1 do begin
      readf, lun, data
    endfor

    data = fltarr(5)
    
    for iwave = 0, NWAVE-1 do begin
      for ibin = 0, NBIN-1 do begin
        readf, lun, data
      print, iwave, ibin, data

        qext(iwave, igroup, ibin) = data[2]
        ssa(iwave, igroup, ibin) = data[3]
        asym(iwave, igroup, ibin) = data[4]
      endfor
    endfor
  endfor

  free_lun, lun


  ; Plot qext, ssa and asym
  !p.multi=[0,1,3]
  loadct, 39

  for iwave = 0, NWAVE-1 do begin
    plot, r[*,0]*1e4, qext[iwave,0,*], yrange=[1e-5,5], xrange=[0.01,15], $
         title = 'Extinction Efficiency, wavelength = '+string(wave[iwave]*1e4)+' (um)', $
         xtitle='Radius [um]', ytitle = 'Qext', thick=3, /XLOG, /YLOG, charsize=2.0

    plot, r[*,0]*1e4, ssa[iwave,0,*], yrange=[0,1], xrange=[0.01,15], $
         title = 'Single Scattering Albedo', $
         xtitle='Radius [um]', ytitle = 'w', thick=3, /XLOG, charsize=2.0

    plot, r[*,0]*1e4, asym[iwave,0,*], yrange=[-1,1], xrange=[0.01,15], $
         title = 'Asymmetry Factor', $
         xtitle='Radius [um]', ytitle = 'g', thick=3, /XLOG, charsize=2.0

    wait, .1
  endfor
  
  plot, r[*,0]*1e4, qext[0,0,*], yrange=[1e-5,5], xrange=[0.01,15], $
       title = 'Extinction Efficiency, wavelength = '+string(wave[0]*1e4)+' to '+string(wave[NWAVE-1]*1e4)+' (um)', $
       xtitle='Radius [um]', ytitle = 'Qext', thick=2, /XLOG, /YLOG, charsize=2.0
  for iwave = 0, NWAVE-1 do begin
    oplot, r[*,0]*1e4, qext[iwave,0,*], thick=2, color=30+4*iwave
  endfor

  plot, r[*,0]*1e4, ssa[0,0,*], yrange=[0,1], xrange=[0.01,15], $
       title = 'Single Scattering Albedo', $
       xtitle='Radius [um]', ytitle = 'w', thick=2, /XLOG, charsize=2.0
  for iwave = 0, NWAVE-1 do begin
    oplot, r[*,0]*1e4, ssa[iwave,0,*], thick=2, color=30+4*iwave
  endfor

  plot, r[*,0]*1e4, asym[0,0,*], yrange=[-1,1], xrange=[0.01,15], $
       title = 'Asymmetry Factor', $
       xtitle='Radius [um]', ytitle = 'g', thick=2, /XLOG, charsize=2.0
  for iwave = 0, NWAVE-1 do begin
    oplot, r[*,0]*1e4, asym[iwave,0,*], thick=2, color=30+4*iwave
  endfor
end
