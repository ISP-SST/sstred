; docformat = 'rst'

;+
; Fit pinhole positions to a rectangular grid. Robust with respect to
; missing (blocked) pinholes.
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl
; 
; 
; 
; 
; :Returns:
;
;    The grid parameters as a structure, { x0:x0, y0:y0, dx:dx, dy:dy, theta:theta, Nx:Nx, Ny:Ny}.
; 
; :Params:
; 
;    im : in, required, type=array
; 
;      A pinhole image.
;   
;   
; 
; 
; :History:
;  
;    2013-09-16 : MGL. First version.
; 
; 
;-
function red_fitgrid, im

  Nblock = 30                   ; Size of area containing a pinhole spot.

  ;; Get starting values for the fit

  ;; X direction
  xproj = total(im, 2)
  xproj = xproj-median(xproj)

  xsz = n_elements(xproj)
  xprojj = xproj
  mx = max(xproj)
  ix = 0
  while max(xproj) gt mx*.2 do begin
     lx = max(xproj, loc)
     indx = loc+indgen(Nblock)-Nblock/2 >0 <(xsz-1)
     peak = xproj[indx]
     cgplot, indx, peak, color = 'red'
     peakfit = MPFITPEAK(indx, peak, A, NTERMS=4)
     cgplot, indx, peakfit, color = 'blue', /over
     if ix eq 0 then xlocs = [A[1]] else xlocs = [xlocs, A[1]]
     xproj[indx] = 0            ; Zero this peak so we won't find it again
     ix += 1
  endwhile   
  cgplot, xprojj
  cgplot, color = 'red', /over, psym = 9, xlocs, xprojj[xlocs]

  xlocs = xlocs(sort(xlocs))
  x0 = xlocs[0]
  dx = median(deriv(xlocs(sort(xlocs))))
  Nx = n_elements(xlocs)

  ;; Y direction

  yproj = total(im, 1)
  yproj = yproj-median(yproj)

  ysz = n_elements(yproj)
  yprojj = yproj
  mx = max(yproj)
  iy = 0
  while max(yproj) gt mx*.2 do begin
     ly = max(yproj, loc)
     indx = loc+indgen(Nblock)-Nblock/2 >0 <(ysz-1)
     peak = yproj[indx]
     cgplot, indx, peak, color = 'red'
     peakfit = MPFITPEAK(indx, peak, A, NTERMS=4)
     cgplot, indx, peakfit, color = 'blue', /over
     if iy eq 0 then ylocs = [A[1]] else ylocs = [ylocs, A[1]]
     yproj[indx] = 0            ; Zero this peak so we won't find it again
     iy += 1
  endwhile   
  cgplot, yprojj
  cgplot, color = 'red', /over, psym = 9, ylocs, yprojj[ylocs]

  ylocs = ylocs(sort(ylocs))
  y0 = ylocs[0]
  dy = median(deriv(ylocs(sort(ylocs))))
  Ny = n_elements(ylocs)

  ;; Rotation angle
  ;theta = atan(dy,dx)-!pi/4.    ; Does not work due to different spacing in x and y.
  theta = 0.01                  ; Expected order of magnitude.


  ;; Now find the 2D coordinates
  imm = im
  ii = 0
  mx = max(imm)
  while max(im) gt mx*0.2 do begin
     ly = max(im, loc)
     xpos = loc mod xsz
     ypos = loc / xsz
     peak = im[xpos-Nblock/2:xpos+Nblock/2-1,ypos-Nblock/2:ypos+Nblock/2-1]

     yfit = mpfit2dpeak(peak, A )

     if ii eq 0 then begin
        xlocs = [xpos-Nblock/2+A[4]]
        ylocs = [ypos-Nblock/2+A[5]]
     endif else begin 
        xlocs = [xlocs, xpos-Nblock/2+A[4]]
        ylocs = [ylocs, ypos-Nblock/2+A[5]]
     endelse

     im[xpos-Nblock/2:xpos+Nblock/2-1,ypos-Nblock/2:ypos+Nblock/2-1] = 0.

     ii += 1
  endwhile   

  tighttv, imm, 0

  parinfo = replicate({ limited:[0,0], limits:[0.D,0] } $
                      , 5)

  parinfo.limited[0] = [1, 1, 1, 1, 0] ; Lower limited
  parinfo.limited[1] = [1, 1, 1, 1, 0] ; Upper limited
  parinfo[0].limits = [0d, dx]         ; x0 limits
  parinfo[1].limits = [0d, dy]         ; y0 limits
  parinfo[2].limits = dx*[0.9, 1.1]    ; dx limits
  parinfo[3].limits = dy*[0.9, 1.1]    ; dy limits

  Pinit = [x0, y0, dx, dy, theta]

  functargs = {X : xlocs, Y : ylocs, Nx : Nx, Ny : Ny}

  print, 'Do mpfit now'
  P = mpfit("red_fitgrid_deviates", Pinit, FUNCTARGS=functargs, PARINFO = parinfo)

  fitgrid = red_gridpoints(p[0], p[1], p[2], p[3], p[4], Nx, Ny)

  tighttv, imm, 1 
  
  ;; Display image
  cgimage, -imm, /axes, /keep

  grid = red_gridpoints(pinit[0], pinit[1], pinit[2], pinit[3], pinit[4], Nx, Ny)

  ;; Plot initial grid
  cgplot, /over, grid.x, grid.y, psym = 9, xrange = [0, xsz], /xstyle, yrange = [0, ysz], /ystyle $
          , color = 'red', symsize = 9

  ;; Plot approximate locations
  cgplot, /over, xlocs, ylocs, psym = 9, color = 'blue'
  
  ;; Plot fitted grid
  cgplot, /over, fitgrid.x, fitgrid.y, psym = 2, xrange = [0, xsz], /xstyle, yrange = [0, ysz], /ystyle $
          , color = 'green', symsize = 2

  print, 'x0,y0 : ', P[0], P[1]
  print, 'dx,dy : ', P[2], P[3]
  print, 'Angle [deg] : ', P[4]*180/!pi

  return, { x0:P[0], y0:P[1], dx:P[2], dy:P[3], theta:P[4], Nx:Nx, Ny:Ny}

end
