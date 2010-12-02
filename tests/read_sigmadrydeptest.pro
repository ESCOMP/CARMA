  openr, lun, 'carma_sigmadrydeptest.txt', /get_lun

  ; Read in the vertical grid.
  readf, lun, nz, nelem
  z    = fltarr(nz)
  dz   = fltarr(nz)

  data = fltarr(3)
  
  for iz = 0, nz-1 do begin
    readf, lun, data
    
    z[iz]  = data[1]
    dz[iz] = data[2]
  endfor

  ; Read in the particles for each time step.
  mmr_  = fltarr(nz, nelem)
  q_    = fltarr(nz, nelem)

  data = fltarr(4)

  nt = 0
  while(not(eof(lun))) do begin
    readf, lun, t1
   
    for ielem = 0, nelem-1 do begin
      for iz = 0, nz-1 do begin
        readf, lun, data
      
        mmr_[iz, ielem]  = data[2]
        q_[iz, ielem]    = data[3]
      endfor
    endfor
   
    if (nt eq 0) then begin
      time = t1
      mmr  = mmr_
      q    = q_
    endif else begin
      time = [time,t1]
      mmr  = [mmr,mmr_]
      q    = [q,q_]
    endelse
   
    nt = nt+1
  endwhile
  
  free_lun, lun

  mmr = reform(mmr,nz,nt,nelem)
  q = reform(q,nz,nt,nelem)
  
  z = z/1000.

  !p.multi = [0,1,2]
  loadct, 39

  ;Calculate the column mass, which should be conserved.
  mass    = fltarr(nelem,nt)
  
  for ielem = 0, nelem-1 do begin
    for it = 0, nt-1 do begin
      mass[ielem,it] = total(q[*,it,ielem]*dz[*])
    endfor
  endfor
  

  for it = 0, nt-1 do begin
    plot, q[*,0,0], z[*], yrange=[1.e-2,15], xrange=[1.e-15,2.e-10], $
         title = 'time = '+string(time[it])+' seconds', $
         xtitle='Particle Concentration [g cm-3]', ytitle = 'Altitude [km]', thick=3, $
	 /xlog, /ylog
 
    ; Add a legend
    plots, [1.5e-10,2.5e-10], 0.1, thick=3, color=66, lin=0
    plots, [1.5e-10,2.5e-10], 0.2, thick=3, color=66, lin=1
    xyouts, 3.0e-10, 0.1, 'DryDep', color=66
    xyouts, 3.0e-10, 0.2, 'NoDryDep', color=66

    for ielem = 0, nelem-1 do begin
      oplot, q[*,it,ielem], z[*], lin=ielem, thick=3, color=66
    endfor

    ; Show the mass evolution.
    plot, mass[0,*], xtitle = 'Time Step', ytitle = 'Column Mass [g cm-2]', thick=6, $
        title = 'Total mass evolution'
    for ielem = 0, nelem-1 do begin
      oplot, mass[ielem,*], thick=6, lin=ielem
    endfor
    for ielem = 0, nelem-1 do begin
      oplot, mass[ielem,0:it], thick=6, color=66, lin=ielem
    endfor

    wait, 0.01
  endfor
  
  plot, q[*,0,0], z[*], yrange=[1.e-2,15], xrange=[1.e-15,2.e-10], $
       title = 'Falling History', $
       xtitle='Particle Concentration [g cm-3]', ytitle = 'Altitude [km]', thick=3, $
       /xlog,/ylog
  oplot, [1.e-15,1.e-9], [8,8], thick=1, lin=1
  xyouts, 2.5e-10, 2.0, 't = 0 sec'

  ; Add a legend
  plots, [1.5e-10,2.5e-10], 0.1, thick=3, lin=0
  plots, [1.5e-10,2.5e-10], 0.2, thick=3, lin=1
  xyouts, 3.e-10, 0.1, 'DryDep'
  xyouts, 3.e-10, 0.2, 'NoDryDep'


  it = 200
  for ielem = 0, nelem-1 do begin
    oplot, q[*,it,ielem], z[*], color=86, thick=3, line=ielem
  endfor
 ;oplot, [1.e-2,2.5e-10], [4,4], color=86, thick=1, lin=1
  xyouts, 2.5e-10, 1.0, 't = 200000 sec', color=86
  
  it = 400
  for ielem = 0, nelem-1 do begin
    oplot, q[*,it,ielem], z[*], color=126, thick=3, line=ielem
  endfor
  oplot, [0,2.5e-10], [0,0], color=126, thick=1, lin=1
  xyouts, 2.5e-10, 0.5, 't = 400000 sec', color=126
  
  plot, mass[0,*], xtitle = 'Time Step', ytitle = 'Column Mass [g cm-2]', thick=6, $
        title = 'Total mass evolution'
  for ielem = 0, nelem-1 do begin
    oplot, mass[ielem,*], thick=6, color=66, lin=ielem
  endfor
end
