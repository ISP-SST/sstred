; docformat = 'rst'

;+
; Calculate the image scale based on pinhole grid pitch measurements. 
;
; The CRISP image scale depends on the prefilter for some years. If no
; measurements exist, just use the default from the config file.
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
;   no_gridfile : in, optional, type=boolean
;   
;     Do not use the pre-measured gridhole pitch values.
; 
; 
; :History:
; 
;   2024-06-17 : MGL. First version. 
; 
;-
function red::imagescale, pref, no_gridfile = no_gridfile

  inam = red_subprogram(/low, calling = inam1)

  grid_pitch_arcsec = self.pinhole_spacing ; ["]

  if keyword_set(no_gridfile) then begin
    print
    print, inam + ' : Prefilter: ' + pref
    print, 'Do not use the gridfile.'
    goto, dont_use_gridfile
  endif

  gridfile = 'calib/refgrid.txt'
  if file_test(gridfile) then begin
    
;    if self.isodate lt '2023-01-01' then begin
;      grid_pitch_arcsec = 5.12  ; ["] The old pinhole array, installed in 2004 (?)
;    endif else begin
;      grid_pitch_arcsec = 4.13  ; ["] The new pinhole array, installed before the 2023 season
;    endelse
    
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
        print, inam + ' : Prefilter: ' + pref
        print, inam + ' : No matching prefilter in ' + gridfile
        goto, dont_use_gridfile
      end

      1 : image_scale = grid_pitch_arcsec / grid_pitch_pixels[indx] ; ["/pix]

      else : image_scale = grid_pitch_arcsec / median(grid_pitch_pixels[indx])
      
    endcase
    
    print
    print, inam + ' : Prefilter: ' + pref
    print, inam + ' : Image scale from config.txt : ' + strtrim(self.image_scale, 2) + '"/pix'
    print, inam + ' : Using the image scale from pinholes : ' + strtrim(image_scale, 2) + '"/pix'
    print
    
    return, image_scale

  endif else begin
    print
    print, inam + ' : File missing: '+gridfile
    print, inam + ' : Prefilter: ' + pref
  endelse

  dont_use_gridfile:
  
  ;; Maybe offer to do fitgrid on the appropriate pinhole image?
  image_scale = red_imagescale_from_pinholes(pref, self.isodate)

  if image_scale ne 0 then begin
    print, inam + ' : Image scale from config.txt : ' + strtrim(self.image_scale, 2) + '"/pix'
    print, inam + ' : Using the image scale from pinhole image : ' + strtrim(image_scale, 2) + '"/pix'
    print
    return, image_scale
  endif

  
;  
;  pfiles = file_search('pinhs/*', count = Nfiles)
;  if Nfiles gt 0 then begin
;    self -> extractstates, pfiles, pstates
;
;    mtch = pstates.prefilter eq pref
;
;
;    indx = where(pstates.is_wb and pstates.camera eq (*self.cameras)[self.refcam] and mtch, Nmatch)
;    if Nmatch eq 0 then begin
;      print, inam + ' : No matching pinhole images.'
;      goto, dont_use_data
;    endif
;  
;    print
;    print, inam + ' : Using pinhole image in ' + pfiles[indx[0]]
;
;    pim = red_readdata(pfiles[indx[0]])
;    fitinfo = red_fitgrid(pim)
;    grid_pitch_pixels = (fitinfo.dx+fitinfo.dy)/2.
;    image_scale = grid_pitch_arcsec / grid_pitch_pixels ; ["/pix]
;    
;    print, inam + ' : Image scale from config.txt : ' + strtrim(self.image_scale, 2) + '"/pix'
;    print, inam + ' : Using the image scale from pinhole image : ' + strtrim(image_scale, 2) + '"/pix'
;    print
;
;    return, image_scale
;
;  endif

  print, inam + ' : Using the image scale from the config file: ' + strtrim(self.image_scale, 2) + '"/pix'
  print
  
  return, self.image_scale

end

if 1 then begin
  
  cd, '/scratch/mats/2024-04-21/CHROMIS'

  a = chromisred("config.txt",/no, /dev)

  pref = '3950'
  image_scale = a -> imagescale(pref)
  print
  image_scale = a -> imagescale(pref, /no_gridfile)
  print
  image_scale = a -> imagescale('2590')

endif else begin

  cd, '/scratch/mats/2018.10.03/CRISP/'

  a = crispred("config.txt",/no, /dev)

  print
  image_scale = a -> imagescale('5896')

  print
  image_scale = a -> imagescale('6173')

  print
  image_scale = a -> imagescale('6302')

  print
  image_scale = a -> imagescale('6563')

  print
  image_scale = a -> imagescale('8542')

  
endelse

  
end
