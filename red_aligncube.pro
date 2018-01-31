; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
; 
; 
; 
; :Returns:
;
;    An array with the shifts needed to align the images in the cube.
; 
; 
; :Params:
; 
;    cub : in, type="array(x,y,t)"
;   
;       A data cube with the sequence of images to align.
;   
;    np : in, optional, type=integer, default=4
;     
;       Length of subcubes to use for alignment.
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
;    aligncube  : in, optional, type=boolean
;   
;      Set this to apply the shifts to the image cube.
;   
; 
; 
; :History:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2014-01-14 : PS. Use red_box_cursor for interactive ROI selection.
;                Default size set to 256.
;
;   2016-10-24 : MGL. Bypass interactive selection of ROI by providing
;                xc and yc.
;
;
;-
function red_aligncube, cub, np $
                        , xbd = xbd, ybd = ybd $
                        , cubic = cubic $
                        , aligncube = aligncube $
                        , xc = xc, yc = yc, centered = centered
  
  if n_elements(xbd) eq 0 then xbd = 256
  if n_elements(ybd) eq 0 then ybd = 256

  inam = 'red_aligncube : '
  dim = size(cub, /dim)

  if keyword_set(centered) then begin
     xc = dim[0]/2
     yc = dim[1]/2
  endif
  if n_elements(xc) gt 0 and n_elements(yc) gt 0 then begin
    sec = [ (xc+[-1, 1]*xbd/2) >0 <dim[0], (yc+[-1, 1]*ybd/2) >0 <dim[1]]
  endif

  if n_elements(np) eq 0 then np = 4
  np = np < dim[2]
  
  ;; split the series into subcubes
  maxj = red_get_iterbounds(dim[2], np)
  shifts = dblarr(2, dim[2])
  k = 0L
  noref = 1
                                ;
  for i = 0L, n_elements(maxj)-1 do begin
    
    red_progressbar, i, n_elements(maxj), 'Calculate image shifts'

    ;; Current sub-cube
    cube = dblarr(dim[0], dim[1], maxj[i])
    rms = fltarr(maxj[i])
    for j = 0, maxj[i] - 1 do begin
      cube[*,*,j] = cub[*,*,k]
      rms[j] = red_getrms(cub[*,*,k])
      K += 1L
    endfor
    ;; Position of the image with highest RMS
    dum = max(rms, idx)

    ;; Align
    if(noref) then begin
      refoffs=[0.0, 0.0]
      if n_elements(sec) eq 0 then begin
        window,0, xs = dim[0], ys = dim[1] $
               , title = 'Select subfield:  LMB moves,  RMB accepts'
        tvscl, red_histo_opt(cube[*,*,idx])
        red_box_cursor, x0, y0, xbd, ybd, /FIXED
        sec = [x0, x0+xbd-1, y0, y0+ybd-1]
      endif
      window,0, xs = xbd, ys = ybd, title = 'Reference for subset'
      tvscl, red_histo_opt(cube[sec[0]:sec[1],sec[2]:sec[3], idx])
      tempoff = red_getcoords(cube[sec[0]:sec[1],sec[2]:sec[3], *], idx)
      shifts[0, 0] = tempoff
      last = maxj[i]
      oldref = cube[sec[0]:sec[1],sec[2]:sec[3], idx] 
      noref = 0
    endif else begin
      nref = cube[sec[0]:sec[1], sec[2]:sec[3], idx]
      tvscl, nref
      tempoff = red_getcoords(cube[sec[0]:sec[1],sec[2]:sec[3], *], idx)
      ;; Alignment between old reference and new reference
      refoffs += red_shc(oldref, nref, /filt, /int) 
      for ii = 0, maxj[i]-1 do tempoff[*, ii] += refoffs
      shifts[0, last] = tempoff
      last += maxj[i]
      oldref = temporary(nref)
    endelse                     ; i
  endfor
  
  ;; Subtract the average
  shifts[0,*] -= median(shifts[0,*])
  shifts[1,*] -= median(shifts[1,*])
  
  if keyword_set(aligncube) then for ii = 0L, dim[2]-1 do begin

    red_progressbar, ii, dim[2], 'Apply the image shifts'

    cub[*, *, ii] = red_shift_im(cub[*, *, ii], shifts[0, ii], shifts[1, ii] $
                                 , cubic = cubic)

  endfor                        ; ii

  return, shifts
  
end
