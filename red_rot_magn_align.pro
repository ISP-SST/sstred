; docformat = 'rst'

;+
; Find rotation, magnification, and alignment parameters that align
; one image to another.
;-

;+
; Apply a parameter array to an image.
; 
; Helper function to red_rot_magn_align.
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
;    The resulting image.
; 
; :Params:
; 
;   im_in : in
; 
;     The image.
; 
;   p : in, type="dblarr(4)"
; 
;     The parameter array.
; 
; :History:
; 
;   2020-04-14 : MGL. First version.
; 
;-
function red_rot_magn_align_apply_p, im_in, p

  im = rot(im_in, p[0], p[1], cubic=-0.5) ; Rotate and magnify
  im = red_shift_sub(im, p[2], p[3])      ; Shift
  return, im
end

;+
; Calculate model vs data deviation.
; 
; Helper function to red_rot_magn_align.
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
;   The deviations.
; 
; :Params:
; 
;   p : in, type="dblarr(4)"
; 
;     The parameter array.
; 
; :Keywords:
; 
;   im_ref : in, type=image
;
;      The reference image.   
;   
;   im_other : in, type=image
; 
;      The image to be aligned to the reference image.
; 
; :History:
; 
;   2020-04-14 : MGL. First version.
; 
;-
function red_rot_magn_align_deviation, p, im_ref = im_ref, im_other = im_other

  im_ret = red_rot_magn_align_apply_p(im_other, p)
  diff = im_ret - im_ref

  ;; Don't evaluate outermost rows and columns
  margin = 20
  dims = long(size(diff, /dim))
  diff = red_centerpic(diff, xs = dims[0]-2*margin, ys = dims[1]-2*margin)
  dims = long(size(diff, /dim))

  return, reform(diff, dims[0]*dims[1]) ; return 1D array

end

;+
; Find rotation, magnification, and alignment parameters that align
; one image to another.
;
; The images have to be fairly close to alignment on input. Just a few
; degrees of rotation, a few percent of relative magnification, and a
; shift of a few pixels. The images also have to have structures that
; are similar in the two images and contain enough information to
; determine a good match.
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
;   Parameters that make the image match the reference image as well
;   as possible: [angle (deg), relative magnification, x shift, y
;   shift]
; 
; :Params:
; 
;   im_in : in, type=image
; 
;      Image to be aligned to the reference image.
; 
;   im_ref_in : in, type=image
; 
;      The reference image.
; 
; :Keywords:
; 
;   display : in, optional, type=boolean
;   
;     Display some image data.
; 
; 
; :History:
; 
;   2020-04-14 : MGL. First version.
; 
;-
function red_rot_magn_align, im_ref_in, im_in, display = display

;  dims = size(im_in, /dim)

  ;; Normalize intensities and remove bias
  im_ref = im_ref_in/median(im_ref_in) - 1
  im     = im_in/median(im_in)         - 1

  ;; Window function
;  sz = dims[0]
;  w = makewindow(sz,softmargin=sz/16.) ; window
  w = 1.
  
  ;; Fourier transform
  fim_ref = fft(im_ref * w)
  fim     = fft(im     * w)

  ;; Fourier spectrum
  spec_ref = abs(fim_ref)
;  spec     = abs(fim)

  ;; Back to image space w/ common Fourier spectrum
  im = float(fft(spec_ref * exp(complex(0, 1)*atan(fim, /phase)), /inv))
  
  if keyword_set(display) then begin
    red_show, im,     w = 0
    red_show, im_ref, w = 1
  endif
  
  ;; Settng up MPFIT 
  Nparam = 4
  parinfo = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0.D], mpside:0}, Nparam)

  parinfo[0:1].limited[0:1] = 1
  parinfo[0].limits[0:1] = [-5, 5]    ; Rotation limits (degrees)
  parinfo[1].limits[0:1] = [0.90,1.1] ; Magnification limits
  parinfo[*].mpside = 2               ; Symmetric derivatives
  parinfo.value = [0d, 1d, 0d, 0d]    ; Initial: zero angle, unit magninfication, zero shifts

  functargs = { im_other:im, im_ref:im_ref }

  p = mpfit('red_rot_magn_align_deviation' $
            , functargs=functargs $
            , parinfo=parinfo $
            , maxiter=500 $
           )
  
  if keyword_set(display) then begin
    ;; Compare reference HMI image to SST image with optimum
    ;; parameters applied.
    opt = red_rot_magn_align_apply_p(im, p)
    red_show, opt, w = 0, /reuse
    blink, [0, 1]
  endif

  return, p
  
end

cd, '/scratch/mats/2016.09.19/CRISP-aftersummer'

restore, 'testdata.sav'
p = red_rot_magn_align(red_pic_at_coord(subim_hmi_before,1030,820,512,512) $
                       , red_pic_at_coord(im_sst,1030,820,512,512), /displ)



end
