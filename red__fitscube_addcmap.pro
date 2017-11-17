; docformat = 'rst'

;+
; Add a cavity map cube to a fitscube file using the mechanisms in the
; WCS distortions draft (Calabretta et al. 2004).
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
; 
; :Params:
; 
;    filename : in, type=string
;
;       The name of the file in which to add the cavity maps.
;
;    cmaps : in, type="fltarr(Nx,Ny,Nscans)"
; 
;       The 3D cube with cavity maps, each adapted to one of the scans
;       in the fitscube file. Unit is nm.
; 
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
;    2017-10-31 : MGL. First version.
; 
; 
; 
; 
;-
pro red::fitscube_addcmap, filename, cmaps

  ;; Keywords numbered with j=3 because WAVE is coordinate 3
  j = '3'

  ;; Modify the main header --------------------------------------------------------------------

  hdr = headfits(filename)

  ;; If the main header doesn't have the EXTEND keyword, add it now.
  red_fitsaddkeyword, hdr, 'EXTEND', !true, 'The file has extension(s).'
  
  red_fitsaddkeyword, hdr, 'CPDIS'+j, 'Lookup', 'WAVE distortions in lookup table', after = 'CDELT'+j
  red_fitsaddkeyword, hdr, 'CPERR'+j, max(abs(cmaps)), 'Max distortion', after = 'CDELT'+j

  ;; Record-valued keyword DPj NAXIS 1-3 because the cavity map cube
  ;; has 3 axes (dimensions), they have to be added together: 
  names = 'DP'+j+' ' + ['NAXES', 'NAXIS.1', 'NAXIS.2', 'NAXIS.3']
  ;; There are 3 axes in the cavity map cube, they correspond to
  ;; axes 1, 2, and 5 in the main HDU:
  values = [3, 1, 2, 5]       
  comments = ['3 axes in the lookup table', 'Spatial X', 'Spatial Y', 'Scan number']
  red_fitsaddkeyword, hdr, names, values, comments, anchor = 'CPDIS'+j

  ;; Construct a header for the image extension with the lookup table. ---------------------------
  
  mkhdr, chdr, cmaps, /image
  anchor = 'DATE'
  red_fitsaddkeyword, anchor = anchor, chdr, 'EXTNAME', 'WCSDVARR', 'WCS distortion array'
  red_fitsaddkeyword, anchor = anchor, chdr, 'EXTVER', 1, 'Distortion array version number'
  red_fitsaddkeyword, anchor = anchor, chdr, 'PCOUNT', 0, 'Special data area of size zero' ; ?? 
  red_fitsaddkeyword, anchor = anchor, chdr, 'GCOUNT', 1, 'One data group'                 ; ??
  red_fitsaddkeyword, anchor = anchor, chdr, 'CRPIX1', 1, 'Distortion array reference pixel' 
  red_fitsaddkeyword, anchor = anchor, chdr, 'CDELT1', 1, 'Grid step size in 1st coordinate' 
  red_fitsaddkeyword, anchor = anchor, chdr, 'CRVAL1', 1, 'Image array pixel coordinate' 
  red_fitsaddkeyword, anchor = anchor, chdr, 'CRPIX2', 1, 'Distortion array reference pixel' 
  red_fitsaddkeyword, anchor = anchor, chdr, 'CDELT2', 1, 'Grid step size in 2nd coordinate' 
  red_fitsaddkeyword, anchor = anchor, chdr, 'CRVAL2', 1, 'Image array pixel coordinate' 
  red_fitsaddkeyword, anchor = anchor, chdr, 'CRPIX3', 1, 'Distortion array reference pixel' 
  red_fitsaddkeyword, anchor = anchor, chdr, 'CDELT3', 1, 'Grid step size in 3rd coordinate'
  red_fitsaddkeyword, anchor = anchor, chdr, 'CRVAL3', 1, 'Image array pixel coordinate'

  red_fitsaddkeyword, anchor = anchor, chdr, ['', 'HISTORY'], ['', 'These wavelength coordinate distortions were generated from the cavity map, shifted, rotated, and cropped as the science data.']

  ;; Write the image extension and the updated header. -----------------------------------------
  
  mgl_fxhmodify, filename, new_header = hdr
  writefits, filename, cmaps, chdr, heap, /append

end
