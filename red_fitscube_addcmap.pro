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
; :Keywords:
; 
;    cmap_number : in, optional, type=integer, default=1
; 
;       Used to distinguish between cavity maps for multiple (CHROMIS)
;       prefilters. The value of the EXTVER keyword in the extensions
;       header, as well as DW3.EXTVER in the main header.
; 
;    prefilter : in, optional, type=integer
; 
;       The prefilter tag.
; 
;    indx : in, optional, type=array, default="all tunings"
; 
;       The tuning pixel coordinate indices for which this cavity map
;       applies. 
; 
; 
; :History:
; 
;    2017-10-31 : MGL. First version.
; 
;    2017-03-22: MGL. Implement extension of (Calabretta et al. 2004). 
; 
;    2018-04-06 : MGL. The Calabretta mechanism is not enough, instead
;                 use new extension of this mechanism with keywords
;                 stored both with the Calabretta record-valued
;                 mechansim and the HIERARCH convention. 
; 
;    2019-09-25 : MGL. Implement writing multiple wavelength
;                 distortions (cavity maps).
; 
;    2019-09-30 : MGL. Make it a regular subroutine, not a class
;                 method. 
; 
;-
pro red_fitscube_addcmap, filename, cmaps $
                          , cmap_number = cmap_number $
                          , prefilter = prefilter $
                          , indx = indx

  inam = red_subprogram(/low, calling = inam1)

  if n_elements(cmap_number) eq 0 then cmap_number = 1
  if n_elements(prefilter) eq 0 then prefilter = ''
  
  ;; Keywords numbered with j=3 because WAVE is coordinate 3
  j = '3'

  ;; Modify the main header --------------------------------------------------------------------

  hdr = headfits(filename)

  ;; Check that dimensions match
  naxis = fxpar(hdr, 'NAXIS*')
  cdims = size(cmaps, /dim)
  if n_elements(cdims) eq 2 then red_append, cdims, 1
  if naxis[0] ne cdims[0] || naxis[1] ne cdims[1] || naxis[4] ne cdims[4] then begin
    ;;  if max(naxis ne cdims) eq 1 then begin
    print, inam + ' : Dimensions do not match'
    print, '   Cube: ', naxis
    print, '   Cmap: ', cdims
    stop
  endif
  
  ;; If the main header doesn't have the EXTEND keyword, add it now.
  red_fitsaddkeyword, hdr, 'EXTEND', !true, 'The file has extension(s).'

  anchor = 'CDELT'+j
  red_fitsaddkeyword, anchor = anchor, hdr, 'CWERR'+j, max(abs(cmaps)), '[nm] Max total distortion'
  red_fitsaddkeyword, anchor = anchor, hdr, 'CWDIS'+j, 'Lookup', 'WAVE distortions in lookup table'

  if n_elements(indx) ne 0 then begin
    ;; Calculate SCALE and OFFSET for the tuning coordinates. Should
    ;; map data cube tuning indices to the range ].5,1.5[, so that
    ;; rounding gives the index 1, the index of the only cmap in this
    ;; WCSDVARR extension. Note that integer tuning indices refer to
    ;; the centers of pixels, so the range of the tuning indices is
    ;; ]indx[0]-0.5,indx[-1]+0.5[. This will include the exact border
    ;; point between tuning pixels twice when there are multiple cmaps
    ;; but we'll have to live with that.
    eps = 1e-3
    scale = (1. - 2.*eps) / ((indx[-1]+.5) - (indx[0]-0.5))
    offset = ( (indx[-1]+.5)*(.5+eps) - (1.5-eps)*(indx[0]-0.5) ) / (1. - 2.*eps)
    print, 'Adding cmap for tun_indx='+red_collapserange(indx)
  endif
  
  ;; Make and write the record-valued DWj keyword. Avoid dots in the
  ;; names, please!

  ;; The HIERARCH representation of the records is an array of lists.
  ;; The lists consist of: The keyword names, the field names, the
  ;; value, the comment. The only list element that can be omitted is
  ;; the comment.
  undefine, hierarch_fields
  if prefilter eq '' then begin
    red_append, hierarch_fields, list('NAME'      , 'Cavity error' ,                'Type of correction'          )
  endif else begin
    red_append, hierarch_fields, list('NAME'      , 'Cavity error for '+prefilter , 'Type of correction'          )
  endelse
  red_append, hierarch_fields, list('EXTVER'      , cmap_number    , 'Extension version number'                   )
  red_append, hierarch_fields, list('NAXES'       , 5              , 'Number of axes in the extension'            )
  red_append, hierarch_fields, list('AXIS1'       , 1              , 'Spatial X'                                  )
  red_append, hierarch_fields, list('AXIS2'       , 2              , 'Spatial Y'                                  )
  red_append, hierarch_fields, list('AXIS3'       , 3              , 'Tuning'                                     )
  red_append, hierarch_fields, list('AXIS4'       , 4              , 'Stokes'                                     )
  red_append, hierarch_fields, list('AXIS5'       , 5              , 'Scan number'                                )
  if n_elements(indx) ne 0 then begin
    red_append, hierarch_fields, list('OFFSET3'   , offset         , 'Tuning coordinates offset'                  )
    red_append, hierarch_fields, list('SCALE3'    , scale          , 'Tuning coordinates scale'                   )
  endif
  red_append, hierarch_fields, list('CWERR'       , max(abs(cmaps)), '[nm] Max distortion (this correction step)' )
  red_append, hierarch_fields, list('CWDIS LOOKUP', 1              , 'Distortions in lookup table'                )
  red_append, hierarch_fields, list('ASSOCIATE'   , 1              , 'Association stage (pixel coordinates)'      )
  red_append, hierarch_fields, list('APPLY'       , 6              , 'Application stage (world coordinates)'      )
  ;; APPLY should be the last keyword

  ;; Translate the hierarch_fields to the kind of record-valued
  ;; keywords defined in the distortions paper.
  for ifield = 1, n_elements(hierarch_fields)-1 do begin ; Skip NAME, strings not allowed
    hfield = red_strreplace((hierarch_fields[ifield])[0],' ','.',n=10)
    ;; Need a dot in some of the record-valued keywords.
    if strmid(hfield, 0, 4) eq 'AXIS'   then hfield = 'AXIS.'+strmid(hfield, 4)
    if strmid(hfield, 0, 4) eq 'SCALE'  then hfield = 'SCALE.'+strmid(hfield, 5)
    if strmid(hfield, 0, 4) eq 'OFFSET' then hfield = 'OFFSET.'+strmid(hfield, 6)
    red_append, names,    hfield
    red_append, values,   (hierarch_fields[ifield])[1]
    red_append, comments, (hierarch_fields[ifield])[2]
  endfor                        ; ifield

  ;; Add the DWj keyword name
  names = 'DW'+j+' ' + names

  oldanchor = anchor

  ;; Make sure to add multiple cmaps in order: the first one will
  ;; remove existing DW3 keywords.
  print, 'nodelete = ', cmap_number NE 1
  
  ;; Write the extended keywords to the header
  red_fitsaddkeyword, anchor = anchor, hdr, names, values, comments, nodelete = cmap_number NE 1
  
  ;; Write the HIERARCH DW3 keywords to the header, delete earlier
  ;; instances if this is the first cmap.
  red_fitsaddkeyword_hierarch, anchor = oldanchor, hdr, 'DW'+j, hierarch_fields, nodelete = cmap_number NE 1
  
  ;; Construct a header for the image extension with the lookup table. ---------------------------
  
  mkhdr, chdr, cmaps, /image
  anchor = 'DATE'
  red_fitsaddkeyword, anchor = anchor, chdr, 'EXTNAME', 'WCSDVARR', 'WCS distortion array'
  red_fitsaddkeyword, anchor = anchor, chdr, 'EXTVER', cmap_number, 'Distortion array version number'
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
  red_fitsaddkeyword, anchor = anchor, chdr, 'CRPIX4', 1, 'Distortion array reference pixel' 
  red_fitsaddkeyword, anchor = anchor, chdr, 'CDELT4', 1, 'Grid step size in 4th coordinate'
  red_fitsaddkeyword, anchor = anchor, chdr, 'CRVAL5', 1, 'Image array pixel coordinate'
  red_fitsaddkeyword, anchor = anchor, chdr, 'CRPIX5', 1, 'Distortion array reference pixel' 
  red_fitsaddkeyword, anchor = anchor, chdr, 'CDELT5', 1, 'Grid step size in 5th coordinate'
  red_fitsaddkeyword, anchor = anchor, chdr, 'CRVAL5', 1, 'Image array pixel coordinate'

  if prefilter eq '' then begin
    history = 'These wavelength coordinate distortions were generated from' $
              + ' the cavity map, shifted, rotated, and cropped as the science data.'
  endif else begin
    history = 'These wavelength coordinate distortions were generated from' $
              + ' the cavity map for the ' + prefilter $
              + ' prefilter, shifted, rotated, and cropped as the science data.'
  endelse 
  
  red_fitsaddkeyword, anchor = anchor, chdr $
                      , ['', 'HISTORY'] $
                      , ['', history  ]

  ;; Write the image extension and the updated header. -----------------------------------------
  ;; modfits, filename, 0, hdr
  red_fitscube_newheader, filename, hdr

  ;;  writefits, filename, cmaps, chdr, heap, /append
  writefits, filename, cmaps, chdr, /append

  
end
