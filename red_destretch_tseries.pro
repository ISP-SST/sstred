; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :author:
; 
; 
; 
; 
; :returns:
; 
; 
; :Params:
; 
;    cub : 
;   
;   
;   
;    platescale : 
;   
;   
;   
;    grids : 
;   
;   
;   
;    clips : 
;   
;   
;   
;    tstep : 
;   
;   
;   
; 
; :Keywords:
; 
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
; 
;-
function red_destretch_tseries, cub, platescale, grids, clips, tstep
                                ;
                                ; T. Bergers (LMSAL) routine, adapted to work with cubes instead of 
                                ; loading images from disk.
                                ;
                                ;
  dim = size(cub, /dim)
  nfiles = dim[2]
  refim = cub[*,*,0]
;
  maxg = MAX(grids)
  maxstr = 2 * platescale       ;2" maximum allowed displacement for any one grid cell.
  sz =  SIZE(refim)
  IF sz[1] GT sz[2] THEN BEGIN
     maxgx =  maxg
     maxgy =  ROUND(FLOAT(sz[2])/sz[1]*maxgx)
  END ELSE BEGIN
     maxgy = maxg
     maxgx = ROUND(FLOAT(sz[1])/sz[2]*maxgy)
  END
  
  delta = FLTARR(dim[2], 2, maxgx, maxgy)

  for i=1, dim[2]-1 do begin
     im = cub[*,*,i]
     dq = dsgridnest(refim, im, grids, clips)
     badg = WHERE(dq gt maxstr, num)
     if num gt 0 then begin
        dq[badg] = 0.0
     end 
     delta[i,*,*,*] = dq
     refim = temporary(im)
  end

  idx = WHERE(~FINITE(delta), num)
  IF num GT 0 THEN delta(idx) = 0.0

;Detrend and unsharp mask the displacements:
  delta = destretch_gridprep(delta, tstep)

  RETURN, delta
END

function red_rotation, img, angle, sdx, sdy, cubic = cubic
                                ;
                                ;
  if n_elements(sdx) eq 0 then sdx = 0.0
  if n_elements(sdy) eq 0 then sdy = 0.0
                                ;
  dim = size(img, /dim)
                                ;
                                ; get the index of each matrix element
                                ;
  xgrid = findgen(dim[0]) # (fltarr(dim[1]) + 1.0)
  ygrid = (fltarr(dim[0]) + 1.0) # findgen(dim[1])
  
                                ;
                                ; get rotation of indexes
                                ; 
  xsi = dim[0] * 0.5
  ysi = dim[1] * 0.5
  dx = cos(angle) * (xgrid - xsi) - sin(angle) * (ygrid - ysi) + xsi - sdx
  dy = sin(angle) * (xgrid - xsi) + cos(angle) * (ygrid - ysi) + ysi - sdy
  
                                ;
                                ; Interpolation onto new grid
                                ;
  
  return, interpolate(img, dx, dy, missi = median(img), cubic = cubic)
end
