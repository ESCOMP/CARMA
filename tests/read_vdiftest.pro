  openr, lun, 'carma_vdiftest.txt', /get_lun

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
  mmr_  = fltarr(nz)
  q_    = fltarr(nz)

  nt = 0
  while(not(eof(lun))) do begin
    readf, lun, t1
   
    for iz = 0, nz-1 do begin
      readf, lun, data
      
      mmr_[iz]  = data[1]
      q_[iz]    = data[2]
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

  mmr = reform(mmr,nz,nt)
  q = reform(q,nz,nt)
  
  z = z/1000.

  !p.multi = [0,1,2]
  loadct, 39

  ;Calculate the column mass, which should be conserved.
  mass    = fltarr(nt)
  
  for it = 0, nt-1 do begin
    mass[it] = total(q[*,it]*dz[*])
  endfor
  
  print
  
  for it = 0, nt-1 do begin
    plot, q[*,0], z[*], yrange=[80,104], xrange=[0,2.5e-10], $
         title = 'time = '+string(time[it])+' seconds', $
         xtitle='Particle Concentration [g cm-3]', ytitle = 'Altitude [km]', thick=3
    oplot, q[*,it], z[*], lin=2, thick=3, color=66

    ; Show the mass evolution.
    plot, mass, xtitle = 'Time Step', ytitle = 'Column Mass [g cm-2]', thick=6, $
        title = 'Total mass evolution'
    oplot, mass[0:it], thick=6, color=66, lin=0

    wait, 0.01
  endfor
  
  plot, q[*,0], z[*], yrange=[80,104], xrange=[0,2.5e-10], $
       title = 'Falling History', $
       xtitle='Particle Concentration [g cm-3]', ytitle = 'Altitude [km]', thick=3
  oplot, [0,2.5e-10], [8,8], thick=1, lin=1
  xyouts, 1.5e-10, 8.2, 't = 0 sec'

  it = 100
  oplot, q[*,it], z[*], color=66, thick=3
  oplot, [0,2.5e-10], [6,6], color=66, thick=1, lin=1
  xyouts, 1.5e-10, 6.2, 't = 100000 sec'
  
  it = 200
  oplot, q[*,it], z[*], color=86, thick=3
  oplot, [0,2.5e-10], [4,4], color=86, thick=1, lin=1
  xyouts, 1.5e-10, 4.2, 't = 200000 sec'
  
  it = 300
  oplot, q[*,it], z[*], color=106, thick=3
  oplot, [0,2.5e-10], [2,2], color=106, thick=1, lin=1
  xyouts, 1.5e-10, 2.2, 't = 300000 sec'
  
  it = 400
  oplot, q[*,it], z[*], color=126, thick=3
  oplot, [0,2.5e-10], [0,0], color=126, thick=1, lin=1
  xyouts, 1.5e-10, 0.2, 't = 400000 sec'
  
  plot, mass, xtitle = 'Time Step', ytitle = 'Column Mass [g cm-2]', thick=6, $
        title = 'Total mass evolution'
  oplot, mass, thick=6, color=66, lin=0
end
