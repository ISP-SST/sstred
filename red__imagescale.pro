; docformat = 'rst'

;+
; Calculate the image scale based on pinhole grid pitch measurements. 
;
; The CRISP image scale depends on the prefilter for some years. If no
; measurements exist, just use the default from the config file.
;
; Try to get a value, using these sources in order: 1.
; self.image_scales (from raw pinholes, written in config file), 2.
; measurements from file (by the pinhole alignment step), 3. measure
; from summed pinholes, 4. self.image_scale. The first 3 are
; individual measurments for different prefilters. The last one is
; common to all prefilters.
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
;   The image plate scale in arcsec per pixel.
; 
; :Params:
; 
;   pref : in, type=string
; 
;     The prefilter for which to calculate the image scale.
; 
; 
; :Keywords:
; 
;   use_config : in, optional, type=boolean
;   
;     Use the image scales from the config file, self.image_scale or
;     self.image_scales, and not those measured for this day.
;
;   use_measurment : in, optional, type=boolean
;   
;     Use measurments based on this day's pinhole data.
; 
;   no_gridfile : in, optional, type=boolean
;   
;     Do not use the pre-measured gridhole pitch values.
; 
;   verbose : in, optional, type=boolean
; 
;     Write messages about what we are doing.
; 
; :History:
; 
;   2024-06-17 : MGL. First version. 
; 
;-
function red::imagescale, pref $
                          , no_gridfile = no_gridfile $
                          , use_config = use_config $
                          , use_measurement = use_measurement $
                          , verbose = verbose

  inam = red_subprogram(/low, calling = inam1)

  if ~keyword_set(use_measurement) && n_elements(self.image_scales) ne 0 then begin
    if ~(n_elements(self.image_scales) eq 1 && self.image_scales eq !null) then begin
      if keyword_set(verbose) then print, inam + ' : Image scales from config.txt : ' $
                                          + strtrim(self.image_scales[pref], 2) + '"/pix'
      return, self.image_scales[pref]
    endif
  endif

  if keyword_set(use_config) then goto, use_config
  
  grid_pitch_arcsec = self.pinhole_spacing ; ["]

  if keyword_set(no_gridfile) then begin
    if keyword_set(verbose) then begin
      print
      print, inam + ' : Prefilter: ' + pref
      print, 'Do not use the gridfile.'
    endif
    goto, dont_use_gridfile
  endif

  gridfile = 'calib/refgrid.txt'
  if file_test(gridfile) then begin
    
    ;; Get the grid pitch in pixels from the gridfile, written during
    ;; pinhole calibration.
    readcol, gridfile, fpi_states, x0, y0, grid_pitch_pixels, angle, COMMENT='#', /NAN $
             , FORMAT = 'A,F,F,F,F', /SILENT, COUNT = Nlines
    
    ;; Select relevant lines, either the same prefilter as pref or the
    ;; matching WB prefilter. I.e., if pref eq '3950', accept all NB
    ;; prefilters matching that WB filter.
    
    prefilters = strmid(fpi_states, 0, 4)
    mtch = self -> match_prefilters(prefilters, pref)
    indx = where(mtch, Nmatch)

    case Nmatch of

      0 : begin
        if keyword_set(verbose) then begin
          print, inam + ' : Prefilter: ' + pref
          print, inam + ' : No matching prefilter in ' + gridfile
        endif
        goto, dont_use_gridfile
      end

      1 : image_scale = grid_pitch_arcsec / grid_pitch_pixels[indx] ; ["/pix]

      else : image_scale = grid_pitch_arcsec / median(grid_pitch_pixels[indx])
      
    endcase
    
    if keyword_set(verbose) then begin
      print
      print, inam + ' : Prefilter: ' + pref
      print, inam + ' : Image scale from config.txt : ' + strtrim(self.image_scale, 2) + '"/pix'
      print, inam + ' : Using the image scale from pinholes : ' + strtrim(image_scale, 2) + '"/pix'
      print
    endif
    
    return, image_scale[0]

  endif else begin
    if keyword_set(verbose) then begin
      print
      print, inam + ' : File missing: '+gridfile
      print, inam + ' : Prefilter: ' + pref
    endif
  endelse

  dont_use_gridfile:
  
  ;; Do fitgrid on the appropriate pinhole image
  image_scale = red_imagescale_from_pinholes(pref, self.isodate)

  if image_scale ne 0 then begin
    if keyword_set(verbose) then begin
      print, inam + ' : Image scale from config.txt : ' + strtrim(self.image_scale, 2) + '"/pix'
      print, inam + ' : Using the image scale from pinhole image : ' + strtrim(image_scale, 2) + '"/pix'
      print
    endif
    return, image_scale
  endif

  use_config:
  if ~keyword_set(use_measurement) && self.image_scale ne 0.0 then begin
    if keyword_set(verbose) then begin
      print, inam + ' : Using the image scale from the config file: ' + strtrim(self.image_scale, 2) + '"/pix'
      print
    endif
    return, self.image_scale
  endif

  ;; Failure
  if keyword_set(verbose) then begin
    print
    print, inam + ' : No value returned.'
  endif
  return, 0
  
end

verbose = 1

if 0 then begin
  
  cd, '/scratch/mats/2024-04-21/CHROMIS'

  a = chromisred("config.txt",/no, /dev)

  pref = '3950'
  image_scale = a -> imagescale(pref, verbose = verbose)
  print
  image_scale = a -> imagescale(pref, /no_gridfile, verbose = verbose)
  print
  image_scale = a -> imagescale('2590', verbose = verbose)

endif else begin

  cd, '/scratch/mats/2016.09.19/CRISP-imscale/'

  a = crispred("config.txt", /no, /dev)

  print
  image_scale = a -> imagescale('6302', verbose = verbose)

  print
  image_scale = a -> imagescale('6563', verbose = verbose)

  print
  image_scale = a -> imagescale('8542', verbose = verbose)
  
endelse

  
end
