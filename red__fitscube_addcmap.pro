; docformat = 'rst'

;+
; Add a cavity map cube to a fitscube file using an extension of the
; mechanisms in the WCS distortions draft (Calabretta et al. 2004).
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
; :History:
; 
;    2017-10-31 : MGL. First version.
; 
;    2017-03-22: MGL. Implement extension of (Calabretta et al. 2004). 
; 
;-
pro red::fitscube_addcmap, filename, cmaps

  ;; Keywords numbered with j=3 because WAVE is coordinate 3
  j = '3'

  ;; Modify the main header --------------------------------------------------------------------

  hdr = headfits(filename)

  ;; If the main header doesn't have the EXTEND keyword, add it now.
  red_fitsaddkeyword, hdr, 'EXTEND', !true, 'The file has extension(s).'

  anchor = 'CDELT'+j
  red_fitsaddkeyword, anchor = anchor, hdr, 'CWERR'+j, max(abs(cmaps)), '[nm] Max total distortion'
  red_fitsaddkeyword, anchor = anchor, hdr, 'CWDIS'+j, 'Lookup', 'WAVE distortions in lookup table'

  ;; Make and write the record-valued DWj keyword, Start with EXTVER,
  ;; which should match the corresponding keyword in the extension
  ;; header.
  names = 'EXTVER'
  values = 1                    ; Default but specify anyway
  comments = 'Extension version number'
  ;; Record-valued keyword DPj NAXIS 1-3 because the cavity map cube
  ;; has 3 axes (dimensions): 1, 2, and 5 in the main HDU.
  red_append, names, ['NAXES', 'AXIS.1', 'AXIS.2', 'AXIS.3']
  red_append, values, [3, 1, 2, 5]
  red_append, comments, ['3 axes in the lookup table' $
                         , 'Spatial X', 'Spatial Y', 'Scan number']
  ;; In the extended distortions notation, the distortions are
  ;; associated to "stage 1" and should be applied at "stage 6".
  red_append, names, ['ASSOCIATE', 'APPLY']
  red_append, values, [1, 6]
  red_append, comments, ['Association stage (pixel coordinates)' $
                         , 'Application stage (world coordinates)']
  ;; Error and distortion mechanism
  red_append, names, ['CWERR', 'CWDIS.LOOKUP']
  red_append, values, [max(abs(cmaps)), 1]
  red_append, comments, ['[nm] Max distortion (this correction step)' $
                         , 'Distortions in lookup table']
  ;; Add the DWj keyword name
  names = 'DW'+j+' ' + names
  ;; Write to the header
  red_fitsaddkeyword, anchor = anchor, hdr, names, values, comments

  ;; Construct a header for the image extension with the lookup table. ---------------------------
  
  mkhdr, chdr, cmaps, /image
  anchor = 'DATE'
  red_fitsaddkeyword, anchor = anchor, chdr, 'EXTNAME', 'WCSDVARR', 'WCS distortion array'
  red_fitsaddkeyword, anchor = anchor, chdr, 'EXTVER', 1, 'Distortion array version number'
  red_fitsaddkeyword, anchor = anchor, chdr, 'PCOUNT', 0, 'Special data area of size zero' 
  red_fitsaddkeyword, anchor = anchor, chdr, 'GCOUNT', 1, 'One data group'
  red_fitsaddkeyword, anchor = anchor, chdr, 'CRPIX1', 1, 'Distortion array reference pixel' 
  red_fitsaddkeyword, anchor = anchor, chdr, 'CDELT1', 1, 'Grid step size in 1st coordinate' 
  red_fitsaddkeyword, anchor = anchor, chdr, 'CRVAL1', 1, 'Image array pixel coordinate' 
  red_fitsaddkeyword, anchor = anchor, chdr, 'CRPIX2', 1, 'Distortion array reference pixel' 
  red_fitsaddkeyword, anchor = anchor, chdr, 'CDELT2', 1, 'Grid step size in 2nd coordinate' 
  red_fitsaddkeyword, anchor = anchor, chdr, 'CRVAL2', 1, 'Image array pixel coordinate' 
  red_fitsaddkeyword, anchor = anchor, chdr, 'CRPIX3', 1, 'Distortion array reference pixel' 
  red_fitsaddkeyword, anchor = anchor, chdr, 'CDELT3', 1, 'Grid step size in 3rd coordinate'
  red_fitsaddkeyword, anchor = anchor, chdr, 'CRVAL3', 1, 'Image array pixel coordinate'

  red_fitsaddkeyword, anchor = anchor, chdr $
                      , ['', 'HISTORY'] $
                      , ['', 'These wavelength coordinate distortions were generated from ' $
                         + 'the cavity map, shifted, rotated, and cropped as the science data.']

  ;; Write the image extension and the updated header. -----------------------------------------
  modfits, filename, 0, hdr

;  writefits, filename, cmaps, chdr, heap, /append
  writefits, filename, cmaps, chdr, /append

end
