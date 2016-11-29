; docformat = 'rst'

;+
; Determine the offsets (image shifts) of an image with respect to a
; reference image.
;
; FFT method is used to maximize the cross-correlation.
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;     J. Chae, May 1999
; 
; 
; :Returns:
; 
;      A two-element array of the offset values defined by OFFSET =
;      (i, j) - (l, m) where (i, j) is the object image coordinates of
;      a feature and (l, m), its reference image coordinates.
; 
; :Params:
; 
;      
;      cor :  out, optional
;
;         The maximum correlation coefficient.
;
;      image : in
;
;         The object image.
;
;      reference : in
;
;         The reference image.
;
; 
; :Keywords:
; 
;      window : in, optional, type=boolean
;   
;         Set this to multiply the images with a Hanning-like window
;         function. 
;
;      fitplane : in, optional, type=boolean
;    
;         Set this to subtract a fitted plane (default is to subtract
;         the median intensity).
;
; 
; :History:
;
;      1999-05-?? : J. Chae. First version (as alignoffset.pro).
;
;      2016-11-22 : MGL. New keywords window and fitplane.
;
;     2016-11-29 : MGL. Renamed into the red_ namespace.
;
; 
;-
function red_alignoffset, image, reference, cor, window = window, fitplane = fitplane

  si = size(image)
  sr = size(reference)
  
  if not(si[1] eq sr[1] and si[2] eq sr[2]) then begin
     print, 'Incompatbile Images : getoffset'
     return, [0,0.]
  endif
  
  if keyword_set(fitplane) then begin
     if sr[1] ne sr[2] then stop ; Only square windows for now
     xp = x_coord(si[1]/2, 1)
     yp = transpose(xp)
     ;; Subtract fitted planes, use planefit from the robust directory
     ;; in astrolib.
     c = planefit(xp, yp, reference, 0., rplane)
     reference1 = reference - rplane
     c = planefit(xp, yp, image, 0., iplane)
     image1 = image - iplane
  endif else begin
     ;; Subtract median values 
     reference1 = reference-median(reference)
     image1     = image-median(image)
  endelse

  if keyword_set(window) then begin
     if sr[1] ne sr[2] then stop ; Only square windows for now
     w = makewindow(sr[1]*[1, 0, 1/8.])
     image1 *= w
     reference1 *= w
  endif

  ;; Maximize cross-correlation over the indeteger pixels
  cor = float(fft(fft(image1, 1)*fft(reference1, -1), -1))
  
  tmp = max(cor, s)
  x0 = s(0) mod si(1)  & x0 = x0 - si(1)*(x0 gt si(1)/2)
  
  y0 = s(0) / si(1)  & y0 = y0 - si(2)*(y0 gt si(2)/2)
  
  ;; Maximize the cross-correlation over the subpixels
  cc = (shift(cor, -x0+1, -y0+1))(0:2, 0:2)
  x1 = (cc(0,1)-cc(2,1))/(cc(2,1)+cc(0,1)-2.*cc(1,1))*.5
  y1 = (cc(1,0)-cc(1,2))/(cc(1,2)+cc(1,0)-2.*cc(1,1))*.5
  x = x0+x1
  y = y0+y1
  if n_params() ge 3 then begin
     image1= shift_sub(image1, -x, -y)
     i=findgen(si(1))#replicate(1., si(2))+x
     j=replicate(1., si(1))#findgen(si(2))+y
     w=i ge 0 and i le si(1)-1 and j ge 0 and j le si(2)-1.
     cor = total(image1*reference1*w)/sqrt(total(image1^2*w)*total(reference1^2*w))
  endif
  return, [x, y]
end

