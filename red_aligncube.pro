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
;    np : 
;   
;   
;   
; 
; :Keywords:
; 
;    xbd  : 
;   
;   
;   
;    ybd  : 
;   
;   
;   
;    cubic  : 
;   
;   
;   
;    aligncube  : 
;   
;   
;   
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2014-01-14 : PS Use red_box_cursor for interactive ROI selection
;                   default size set to 256
;-
function red_aligncube, cub, np, xbd = xbd, ybd = ybd, cubic = cubic, aligncube = aligncube
  if n_elements(xbd) eq 0 then xbd = 256
  if n_elements(ybd) eq 0 then ybd = 256

  inam = 'red_aligncube : '
  dim = size(cub, /dim)
  
  if n_elements(np) eq 0 then np = 4
  np = np < dim[2]
  
   ;; split the series into subcubes
  maxj = red_get_iterbounds(dim[2], np)
  shifts = dblarr(2, dim[2])
  k = 0L
  noref = 1
                                ;
  for i = 0L, n_elements(maxj)-1 do begin
       ;;; current sub-cube
     cube = dblarr(dim[0], dim[1], maxj[i])
     rms = fltarr(maxj[i])
     for j = 0, maxj[i] - 1 do begin
        cube[*,*,j] = cub[*,*,k]
        rms[j] = red_getrms(cub[*,*,k])
        K += 1L
     endfor
        ; Position of the image with highest RMS
     dum = max(rms, idx)
                                 ; Align
     if(noref) then begin
        refoffs=[0.0, 0.0]
        window,0, xs = dim[0], ys = dim[1], title = 'Select subfield:  LMB moves,  RMB accepts'
        tvscl, histo_opt(cube[*,*,idx])
        red_box_cursor, x0, y0, xbd, ybd, /FIXED
        sec = [x0, x0+xbd-1, y0, y0+ybd-1]
        window,0, xs = xbd, ys = ybd, title = 'Reference for subset'
        tvscl, histo_opt(cube[sec[0]:sec[1],sec[2]:sec[3], idx])
        tempoff = red_getcoords(cube[sec[0]:sec[1],sec[2]:sec[3], *], idx)
        shifts[0, 0] = tempoff
        last = maxj[i]
        oldref = cube[sec[0]:sec[1],sec[2]:sec[3], idx] 
        noref = 0
     endif Else begin
        nref = cube[sec[0]:sec[1], sec[2]:sec[3], idx]
        tvscl, nref
        tempoff = red_getcoords(cube[sec[0]:sec[1],sec[2]:sec[3], *], idx)
        refoffs += red_shc(oldref, nref, /filt, /int) ;Alignment between old reference and new reference
        FOR ii = 0, maxj[i]-1 DO tempoff[*, ii] += refoffs
        shifts[0, last] = tempoff
        last += maxj[i]
        oldref = temporary(nref)
     endelse
  endfor
  
      ;Substract the average
  shifts[0,*]-= median(shifts[0,*])
  shifts[1,*]-= median(shifts[1,*])
  
  if keyword_set(aligncube) then for ii = 0L, dim[2]-1 DO $
    cub[*, *, ii] = red_shift_im(cub[*, *, ii], shifts[0, ii], shifts[1, ii], cubic = cubic)

  return, shifts
end
