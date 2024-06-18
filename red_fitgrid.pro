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
; :Returns:
;
;    The grid parameters as a structure, { x0:x0, y0:y0, dx:dx, dy:dy,
;    theta:theta}, where (X0,y0) are the pixel coordinates of one
;    point in the grid, dx and dy are the grid spacings in X and Y
;    resp., and theta is a rotation angle in radians.
; 
; :Params:
; 
;    im : in, required, type=array
; 
;      A pinhole image.
;   
; :Keywords:  
; 
;    debug : in, optional, type=boolean
; 
;      Use when debugging.
; 
; :History:
;  
;    2013-09-16 : MGL. First version.
;  
;    2019-12-06 : MGL. New keywords dx_init and dy_init. Remove
;                 pinholes near array limits. Improve the initial
;                 fitting in X and Y.
; 
;    2024-06-18 : MGL. Re-implemented
; 
;-
function red_fitgrid, im, dx_init = dx_init, dy_init = dy_init, debug = debug

  dims = size(im, /dim)
  x = indgen(dims[0])#replicate(1, dims[1])
  y = indgen(dims[1])##replicate(1, dims[0])
  
  ;; Segment the pinholes
  seg = red_separate_mask(im GE max(im)/5)
  Nholes = max(seg)
  if Nholes lt 10 then begin
    return, { x0:0d, y0:0d, dx:0d, dy:0d, theta:0d}
  endif

  xy = fltarr(2, Nholes)
  for ihole = 0, Nholes-1 do begin
    ;; If this is not good enough, we could run 2d peak fitting on the
    ;; individual pinhole spots.  
    indx = where(seg eq ihole+1)
;    xy[0, ihole] = mean(x[indx])
;    xy[1, ihole] = mean(y[indx])
;    print, xy[*, ihole]
    ;;centroids
    xy[0, ihole] = total(x[indx]*im[indx])/total(im[indx])
    xy[1, ihole] = total(y[indx]*im[indx])/total(im[indx])
;    print, xy[*, ihole]
  endfor                        ; ihole

  ;; Remove spots too close to the border.
  margin = 40
  indx = where(xy[0,*] ge margin  and $
               xy[0,*] lt dims[0]-margin and $
               xy[1,*] ge margin and $
               xy[1,*] lt dims[1]-margin $
               , Nwhere)
  if Nwhere eq 0 then stop
  xy = xy[*, indx]

  ;; Display image
  cgimage, im, /axes, /keep

  ;; Plot approximate locations
  cgplot, /over, xy[0, *], xy[1, *], psym = 9, color = 'cyan', symsize = 4

  ;; Find coordinates of a central spot
  tmp = min(sqrt((xy[0,*]-mean(xy[0,*]))^2 + (xy[1,*]-mean(xy[1,*]))^2),ipos)
  x0 = xy[0,ipos]
  y0 = xy[1,ipos]

  cgplot, /over, [x0], [y0], psym = 2, color = 'orange', symsize=2

  sim = smooth((seg gt 0)*im,11)

  ;; Two passes to get initial dx,dy,theta
  rim = im*(sim gt 0)
  xproj = total(rim, 2)
  yproj = total(rim, 1)
  xpos = where(red_differential(xproj gt max(xproj)/10.) gt 10)
  ypos = where(red_differential(yproj gt max(yproj)/10.) gt 10)

  dx = median(red_differential(xpos))
  dy = median(red_differential(ypos))

  ;; Rotation angle
;  indx = where(abs(xy[0,*] - x0) lt dx/2)
;  cgplot, /over, xy[0,indx], xy[1,indx], psym = 4, symsize=2, color = 'pink'
;  theta = atan(max(xy[1,indx])-min(xy[1,indx]),max(xy[0,indx])-min(xy[0,indx])) - !pi/2
  indx =  where(abs(xy[1,*] - y0) lt dy/2)
  cgplot, /over, xy[0,indx], xy[1,indx], psym = 4, symsize=2, color = 'pink'
                                ;theta = atan(max(xy[1,indx])-min(xy[1,indx]),max(xy[0,indx])-min(xy[0,indx]))
  theta = atan(xy[1,indx[-1]]-xy[1,indx[0]],xy[0,indx[-1]]-xy[0,indx[0]])
  if abs(theta - !pi) lt abs(theta) then theta = theta - !pi
  if abs(theta + !pi) lt abs(theta) then theta = theta + !pi
  
  if 0 then begin
    rim = rot(im*(sim gt 0),theta*180/!pi)
    xproj = total(rim, 2)
    yproj = total(rim, 1)
    xpos = where(red_differential(xproj gt max(xproj)/10.) gt 10)
    ypos = where(red_differential(yproj gt max(yproj)/10.) gt 10)

    dx = median(red_differential(xpos))
    dy = median(red_differential(ypos))

    ;; Rotation angle
    indx =  where(abs(xy[1,*] - y0) lt dy/2)
    cgplot, /over, xy[0,indx], xy[1,indx], psym = 12, symsize=2, color = 'pink'
    theta = atan(max(xy[1,indx])-min(xy[1,indx]),max(xy[0,indx])-min(xy[0,indx])) 
  endif
  
  Nx = n_elements(xpos) + 2
  Ny = n_elements(ypos) + 2

  
  
  ;;
;  x0 = min(xy[0,*])
;  y0 = min(xy[1,*])

  
;  ;; X direction
;  xproj = total(sim, 2)
;  yproj = total(sim, 1)
;  ;;xproj = smooth(xproj-median(xproj), dx_init/8)
;
;  xpos = where(red_differential(xproj gt max(xproj)/10.) gt 10)
;  ypos = where(red_differential(yproj gt max(yproj)/10.) gt 10)
; dx = median(red_differential(xpos)) 
; dy = median(red_differential(ypos)) 
;

  parinfo = replicate({ limited:[0,0], limits:[0.D,0] } $
                      , 5)

  parinfo.limited[0] = [1, 1, 1, 1, 1]    ; Lower limited
  parinfo.limited[1] = [1, 1, 1, 1, 1]    ; Upper limited
                                ;parinfo[0].limits = [0d, dx_init*1.5]    ; x0 limits
;;parinfo[1].limits = [0d, dy_init*1.5]    ; y0 limits
  parinfo[0].limits = x0+dx*[-1d, 1d]/2   ; x0 limits
  parinfo[1].limits = y0+dy*[-1d, 1d]/2   ; y0 limits
  parinfo[2].limits = dx*[0.9, 1.1]       ; dx limits
  parinfo[3].limits = dy*[0.9, 1.1]       ; dy limits
  parinfo[4].limits = 5.*[-1, 1] *!pi/180. ; angle limits +/- 5 deg

  Pinit = [x0, y0, dx, dy, theta]
  grid = red_gridpoints(pinit[0], pinit[1], pinit[2], pinit[3], pinit[4], Nx, Ny)

;; Plot initial grid
  cgplot, /over, grid.x, grid.y, psym = 9, color = 'green', symsize = 6

;  stop

  functargs = {X : reform(xy[0, *]), Y : reform(xy[1, *])}

                                ;help, red_fitgrid_deviates(Pinit, x = reform(xy[0, *]), y = reform(xy[1, *]))
  
  print, 'Do mpfit now'
  P = mpfit("red_fitgrid_deviates", Pinit, FUNCTARGS=functargs, PARINFO = parinfo, status = status, errm = errm, quiet = ~keyword_set(debug))

  if errm ne '' then begin
    print, errm
    return, { x0:0d, y0:0d, dx:0d, dy:0d, theta:0d}
    stop
  endif
  
  fitgrid = red_gridpoints(p[0], p[1], p[2], p[3], p[4], Nx, Ny)

;  red_show, imm, w = 1 


;; Plot approximate locations
;  cgplot, /over, xlocs, ylocs, psym = 9, color = 'cyan', symsize = 5

;; Plot fitted grid
  cgplot, /over, fitgrid.x, fitgrid.y, psym = 9, color = 'yellow', symsize = 3

;  gridspacing_in_arcsec = 5.12  ; From the CRISPRED paper
  print, 'Initial:'
  print, 'x0,y0 : ', Pinit[0], Pinit[1]
  print, 'dx,dy : ', Pinit[2], Pinit[3]
  print, 'Angle [deg] : ', Pinit[4]*180/!pi

  print, 'Fitted:'
  print, 'x0,y0 : ', P[0], P[1]
  print, 'dx,dy : ', P[2], P[3]
  print, 'Angle [deg] : ', P[4]*180/!pi

  if keyword_set(debug) then stop

  p = float(p)
  return, { x0:P[0], y0:P[1], dx:P[2], dy:P[3], theta:P[4]}

end

cd, '/scratch/mats/2016.09.19/CRISP-aftersummer'


cams = 'cam' + ['XX', 'XIX', 'XXV']

for icam = 0, 2 do begin

  dark = file_search('darks/'+cams[icam]+'.dark.fits')

  pinhs = file_search('pinhs/'+cams[icam]+'_6302_6302_+0.pinh.fits', count = Nfiles)

  im = red_readdata(pinhs[0]) - red_readdata(dark[0])
  

  p = red_fitgrid(im)

  stop
  
endfor


end
