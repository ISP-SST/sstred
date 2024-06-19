; docformat = 'rst'

;+
; Calculate the image scale by use of pinhole grid pitch measurements.
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
;   The image scale in arcsec/pixel. (0 if failed.)
; 
; :Params:
; 
;    date : in, type=string
; 
;      The date in ISO format YYYY-MM-DD.
; 
;    pref : in, type=string
; 
;      The prefilter.
; 
; :Keywords:
; 
;   fitinfo : out, optional, type=struct
;
;     The output from red_fitgrid.
; 
; :History:
; 
;   2024-06-19 : MGL. First version.
; 
;-
function red_imagescale_from_pinholes, pref, isodate, fitinfo = fitinfo

  inam = red_subprogram(/low, calling = inam1)
  
  grid_pitch_arcsec = red_pinhole_pitch_arcsec(isodate)

  pfiles = file_search('pinhs/*', count = Nfiles)
  if Nfiles eq 0 then return, 0

  cameras = strtrim(red_fitsgetkeyword_multifile(pfiles, 'CAMERA'), 2)
  prefilters = strtrim(red_fitsgetkeyword_multifile(pfiles, 'FILTER1'), 2)

  is_wb = ((strsplit(cameras, '-', /extract)).toarray())[*,1] eq 'W'

  indx = where(is_wb and prefilters eq pref, Nmatch)
  if Nmatch eq 0 then begin
    print, inam + ' : No matching pinhole images.'
    return, 0
  endif
    
  print, inam + ' : Using pinhole image in ' + pfiles[indx[0]]

  pim = red_readdata(pfiles[indx[0]])
  fitinfo = red_fitgrid(pim)
  grid_pitch_pixels = (fitinfo.dx+fitinfo.dy)/2.
  image_scale = grid_pitch_arcsec / grid_pitch_pixels ; ["/pix]
  
  return, image_scale
  
end

if 1 then begin
  
  cd, '/scratch/mats/2024-04-21/CHROMIS'
  isodate = '2024-04-21'
  
;  a = chromisred("config.txt",/no, /dev)

  pref = '3950'
  image_scale = red_imagescale_from_pinholes(pref, isodate)
  print, image_scale
  print
  image_scale = red_imagescale_from_pinholes('2590', isodate)
  print, image_scale

endif else begin

  cd, '/scratch/mats/2018.10.03/CRISP/'
  isodate = '2018-10-03'

;  a = crispred("config.txt",/no, /dev)

  print
  image_scale = red_imagescale_from_pinholes('5896', isodate)
  print, image_scale
  
  print
  image_scale = red_imagescale_from_pinholes('6173', isodate)
  print, image_scale

  print
  image_scale = red_imagescale_from_pinholes('6302', isodate)
  print, image_scale

  print
  image_scale = red_imagescale_from_pinholes('6563', isodate)
  print, image_scale

  print
  image_scale = red_imagescale_from_pinholes('8542', isodate)
  print, image_scale

  
endelse


end
