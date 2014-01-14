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
; 
;-
function red_aligncube, cub, np, xbd = xbd, ybd = ybd, cubic = cubic, aligncube = aligncube
  if n_elements(xbd) eq 0 then xbd = 255
  if n_elements(ybd) eq 0 then ybd = 255

  inam = 'red_aligncube : '
  dim = size(cub, /dim)
  if n_elements(np) eq 0 then np = 4


  np = np < dim[2]
  
                                ;
                                ; start cycling
                                ;
  maxj = red_get_iterbounds(dim[2], np)
  shifts = dblarr(2, dim[2])
  k = 0L
  noref = 1
                                ;
  for i = 0L, n_elements(maxj)-1 do begin
     cube = dblarr(dim[0], dim[1], maxj[i])
     namt = strarr(maxj[i])
     rms = fltarr(maxj[i])
     for j = 0, maxj[i] - 1 do begin
        cube[*,*,j] = cub[*,*,k]
        rms[j] = red_getrms(cub[*,*,k])
        K += 1L
     endfor
                                ;
                                ; Position of the image with highest RMS
                                ;
     dum = max(rms, idx)
                                ;
                                ; Align
                                ;
     if(noref) then begin
        refoffs=[0.0, 0.0]
        window,0, xs = dim[0], ys = dim[1], title = 'Click where you want to place the reference subfield'
        tvscl, histo_opt(cube[*,*,idx])
        cursor, x, y, /device, /down 
        sec = [(x - xbd/2)>0, (x + xbd/2)<(dim[0]-1), (y - ybd/2)>0, (y + ybd/2)<(dim[1]-1)]
        nx = sec[1] - sec[0] + 1
        ny = sec[3] - sec[2] + 1 
        window,0, xs = nx, ys = ny, title = 'Reference for subset'
        tvscl, histo_opt(cube[sec[0]:sec[1],sec[2]:sec[3], idx])
        tempoff = red_getcoords(cube[sec[0]:sec[1],sec[2]:sec[3], *], idx)

        dim2 = size(tempoff, /dimension)
        shifts[*,0:dim2[1]-1] = tempoff
        last = dim2[1]
        oldref = cube[sec[0]:sec[1],sec[2]:sec[3], idx] 
        noref = 0
     endif Else begin
        nref = cube[sec[0]:sec[1],sec[2]:sec[3],idx]
        tempoff = red_getcoords(cube[sec[0]:sec[1],sec[2]:sec[3], *], idx)
        refoffs+=red_shc(oldref,nref, /filt, /int) ;Alignment between old reference and new reference
        dim2=size(tempoff,/dimension)
        if n_elements(dim2) eq 1 then shifts[*,last]=refoffs $ ;Case when it is only one image in the subset
        else begin
           shifts[*,last:last+dim2[1]-1] = tempoff + $
                                           transpose([[replicate(refoffs[0],dim2[1])],[replicate(refoffs[1],dim2[1])]])
           last+=dim2[1]
        endelse
        oldref = temporary(nref)
     endelse
  endfor
  
                                ;
                                ;Substracts the average
                                ;
  shifts[0,*]-= median(shifts[0,*])
  shifts[1,*]-= median(shifts[1,*])
  
                                ;
                                ; Want to align images in the cube?
                                ;
  if keyword_set(aligncube) then begin
     for ii = 0L, dim[2]-1 do cub[*,*,ii] = red_shift_im(cub[*,*,ii], shifts[0,ii], shifts[1,ii], cubic = cubic)
  endif
  return, shifts
end
