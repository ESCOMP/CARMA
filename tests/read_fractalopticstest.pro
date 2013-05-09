  openr, lun, 'carma_fractalopticstest.txt', /get_lun


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
  
  qext_fractal   = fltarr(NWAVE,NGROUP,NBIN)
  ssa_fractal    = fltarr(NWAVE,NGROUP,NBIN)
  asym_fractal   = fltarr(NWAVE,NGROUP,NBIN)

  qext_sphere   = fltarr(NWAVE,NGROUP,NBIN)
  ssa_sphere    = fltarr(NWAVE,NGROUP,NBIN)
  asym_sphere   = fltarr(NWAVE,NGROUP,NBIN)

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

    data = fltarr(8)
    
    for iwave = 0, NWAVE-1 do begin
      for ibin = 0, NBIN-1 do begin
        readf, lun, data
        ;print, wave(iwave)*1.0e4,iwave, ibin, data

        qext_fractal(iwave, igroup, ibin) = data[2]
        ssa_fractal(iwave, igroup, ibin) = data[3]
        asym_fractal(iwave, igroup, ibin) = data[4]
        qext_sphere(iwave, igroup, ibin) = data[5]
        ssa_sphere(iwave, igroup, ibin) = data[6]
        asym_sphere(iwave, igroup, ibin) = data[7]
      endfor
    endfor
  endfor

  free_lun, lun
  wave(*) = wave(*) * 1.0e4

  ; Plot qext, ssa and asym
  !p.multi=[0,1,3]
  loadct, 39

  plot, wave(*), qext_fractal(*,0,0), yrange=[0,26], xrange=[0.0,1.0], $
       title = 'Extinction Efficiency', $
       xtitle='Wavelength [um]', ytitle = 'Qext', thick=2, charsize=2.0, /nodata
  oplot, wave(*), qext_fractal(*,0,0), linestyle=0,color=200, thick=2.0
  oplot, wave(*), qext_sphere(*,0,0), linestyle=0, thick=2.0
  xyouts, 0.66,0.92, "For Titan hydrocarbon hazes",/NORMAL, charsize=1.1
  xyouts, 0.75,0.89, "Fractal", color=200, /NORMAL
  xyouts, 0.75,0.86, "Sphere", /NORMAL


  plot, wave(*), ssa_fractal[*,0,0], yrange=[0,1], xrange=[0.0,1.0], $
       title = 'Single Scattering Albedo', $
       xtitle='Wavelength [um]', ytitle = 'w', thick=2,  charsize=2.0, linestyle=2, /nodata
  oplot, wave(*), ssa_fractal(*,0,0), linestyle=0, color=200, thick=2.0
  oplot, wave(*), ssa_sphere(*,0,0), linestyle=0, thick=2.0

  plot, wave(*), asym_fractal[*,0,0], yrange=[-1.05,1.05], xrange=[-0.0,1.0], ystyle=1, $
       title = 'Asymmetry Factor', $
       xtitle='Wavelength [um]', ytitle = 'g', thick=2, charsize=2.0, linestyle=2, /nodata
  oplot, wave(*), asym_fractal(*,0,0),linestyle=0, color=200, thick=2.0
  oplot, wave(*), asym_sphere(*,0,0),linestyle=0, thick=2.0

end
