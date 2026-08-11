; docformat = 'rst'

;+
; Find the linear shift between two 2D images, evaluating only a
; region defined by a mask.
;
; Call much like red_shc but always interpolate so that keyword is not
; defined. Uses squared difference function rather than FFT cross
; correlation and uses 2D quadratic interpolation.
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, Institute for Solar Physics
; 
; 
; :Returns:
; 
;     The shifts [dx,dy] of im relative to imref.
; 
; :Params:
; 
;     imref_in : in, type=array
; 
;       The reference image.
; 
;     im_in : in, type=array
; 
;       The other image.
; 
;     mask_in : in, type=array
; 
;       The mask.
; 
; 
; :Keywords:
; 
;   
;   
;   
; 
; 
; :History:
; 
;   2021-06-30 : MGL. First version.
; 
;   2024-07-11 : MGL. A version that actually works.
; 
;-
function red_shc_mask, imref_in, im_in, mask_in $ 
                       , poly = pfit $
                       , range = range $
                       , status = status $
                       , hipass = hipass
  
  dims = size(imref_in, /dim)
  
  if n_elements(range) eq 0 then range = min(dims)/10
  if n_elements(mask_in) eq 0 then mask_in = replicate(1, dims[0], dims[1])

  ;; Close holes in the mask
  mask = red_centerpic(mask_in gt 0, xs = dims[0]+20, ys = dims[1]+20)
  structure = replicate(1, 5, 5) 
  mask = red_centerpic(morph_close(mask, structure), xs = dims[0], ys = dims[1])
  
  ;; Erode mask by max_shift pixels
  mask = morph_distance(mask,neigh=3) gt range
  
  imref = imref_in    
  im    = im_in    
  
  ;; Is a removal of polynomial fits requested?
  if n_elements(pfit) gt 0 then begin
    imref = imref - sfit(imref, pfit)
    im    = im    - sfit(im,    pfit)
  endif

  ;; Highpass filter?
  if keyword_set(hipass) then begin
    
    im    = im    - smooth(im,25)      
    imref = imref - smooth(imref,25)
    
  endif
  
  ;; Calculate abs diff squared of im and imref, shifted by
  ;; -max_shift,...,0,...max_shift pixels. Find minimum. Interpolate
  ;; to get subpixel shifts.
  cc = fltarr(2*range+1, 2*range+1)
  for ix = -range, range do $  
     for iy = -range, range do $  
        cc[ix+range, iy+range] = total(abs(mask*(imref-shift(im, ix, iy))^2))
  
  tmp = min(cc, indx)
  pos = array_indices(cc, indx) ; The integer position of the minimum

  if max(pos) eq 2*range or min(pos) eq 0 then begin
    status = 1
    return, [0, 0]
  endif 
  
  ;; Interpolate
  a2 = ( cc[pos[0]+1, pos[1]+0] - cc[pos[0]-1, pos[1]+0] ) / 2.
  a3 = ( cc[pos[0]+1, pos[1]+0] - 2*cc[pos[0], pos[1]] + cc[pos[0]-1, pos[1]+0] ) / 2.
  a4 = ( cc[pos[0]+0, pos[1]+1] - cc[pos[0]+0, pos[1]-1] ) / 2.
  a5 = ( cc[pos[0]+0, pos[1]+1] - 2*cc[pos[0], pos[1]] + cc[pos[0]+0, pos[1]-1] ) / 2.
  a6 = ( cc[pos[0]+1, pos[1]+1] - cc[pos[0]-1, pos[1]+1] - cc[pos[0]+1, pos[1]-1] + cc[pos[0]-1, pos[1]-1] ) / 4.
  
  dx = pos[0] + (2*a2*a5 - a4*a6)/(a6^2 - 4*a3*a5 )
  dy = pos[1] + (2*a3*a4 - a2*a6)/(a6^2 - 4*a3*a5 )

  status = 0
  return, [dx, dy] - range
  
  
;; Set status keyword to ... for success and ... for failure (when
;; the minimum is on the edge).
  
end

sz = 200

imref = red_centerpic(red_readdata('/scratch/mats/2016.09.19/CHROMIS-jan19/momfbd/09:28:36/3950/cfg/results/camXXVIII_2016-09-19T09:28:36_00000_3950.momfbd'), sz = sz)


mask = makewindow(sz, z = 30) and ~round(rot(makewindow(sz, z = 40), 20))



dx =  2.3
dy = -1.2

im = shift_sub(imref, -dx, -dy)

shifts = red_shc(imref, im, /int, /filt, range = 5)
shifts2 = red_shc_mask(imref, im, mask, range = 5, poly = 1)

print, 'Real shifts', dx, dy
print, 'red_shc: ', shifts
print, 'shc_mask: ', shifts2

end
