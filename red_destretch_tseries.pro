; docformat = 'rst'

;+
; Destretch a time sequence of images in order to remove rubber sheet
; stretching caused by differential seeing.
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :author:
; 
;    Tom Bergers (LMSAL) routine, adapted to work with cubes instead
;    of loading images from disk.
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
  
  dim = size(cub, /dim)
  nfiles = dim[2]
  refim = cub[*,*,0]
;
  maxg = max(grids)
  maxstr = 2 * platescale       ;2" maximum allowed displacement for any one grid cell.
  sz =  SIZE(refim)
  IF sz[1] GT sz[2] THEN BEGIN
     maxgx =  maxg
     maxgy =  round(float(sz[2])/sz[1]*maxgx)
  END ELSE BEGIN
     maxgy = maxg
     maxgx = round(float(sz[1])/sz[2]*maxgy)
  END
  
  delta = fltarr(dim[2], 2, maxgx, maxgy)

  for i=1, dim[2]-1 do begin
     im = cub[*,*,i]
     dq = dsgridnest(refim, im, grids, clips)
     badg = where(dq gt maxstr, num)
     if num gt 0 then begin
        dq[badg] = 0.0
     end 
     delta[i,*,*,*] = dq
     refim = temporary(im)
  end

  idx = where(~finite(delta), num)
  IF num GT 0 THEN delta(idx) = 0.0

  ;;Detrend and unsharp mask the displacements:
  delta = destretch_gridprep(delta, tstep)

  return, delta
end

