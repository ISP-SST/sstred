; docformat = 'rst'

;+
; Generate good default tiles and clips stretch parameters.
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
;   dims : in, type=array
; 
;     Size of the FOV in pixels, [Nx,Ny].
; 
;   image_scale : in, type=float
; 
;     The image scale in arcsec/pixel
; 
; :Keywords:
; 
;   clips : out, optional, type=array
; 
;      
;
;   mode : in, optional, type=string, default='old_crisp'
; 
;
;
;   temporal : in, optional, type=boolean
; 
;
;
;   tiles : out, optional, type=array
; 
;
; :History:
; 
;      2026-03-06 : MGL. First version.
; 
;-
pro red_tiles_and_clips, dims, image_scale $
                         , clips = clips $
                         , mode = mode $
                         , tiles = tiles $
                         , temporal = temporal $
                         , verbose = verbose
  
  ;;  For now, we'll start from old defaults from CRISP/Sarnoff and
  ;;  scale them to the FOV of the instrument/cameras at hand. We may
  ;;  want to change this later.
  if n_elements(mode) eq 0 then mode = 'old_crisp'

  case mode of

    'old_crisp' : begin

      ;;
      ;;  CRISP in 2019:
      ;;
      ;;  Sarnoff = [1024,1024]
      ;;  arcsecperpix =  0.0592
      ;;
      ;;  --> FOV = 1024 * 0.0592" = 60.6"
      ;;
      ;;  The units of clips as well as dim/tiles is pixels.
      ;;
      ;;  "tiles" sets the number of cross-correlation subfields along the
      ;;  (longer) dimension of the FOV. So the size in arcsec of each
      ;;  subfield is approx. tiles / arcsecperpix.
      ;;
      ;;  "clips" limits the length of the increment of the stretch vectors in
      ;;  each iteration.

      if keyword_set(temporal) then begin
        ;; This is for time-series destretching, where the images are not
        ;; supposed to be identical. From old method red__polish_tseries.
        tiles_crisp = [6, 8, 14, 24]
        clips_crisp = [12, 4, 2, 1]
      endif else begin
        ;; This is for destretching of images that are supposed to be
        ;; identical. From red__measure_data_destretch.
        tiles_crisp = [5, 12, 32, 48, 64] 
        clips_crisp = [32, 16,  8,  4,  2]
      endelse
      sz_crisp = 1024.          ; Sarnoff detector size
      imscale_crisp = 0.0592 

      ;; Decouple from CRISP
      t = (sz_crisp * imscale_crisp) / tiles_crisp ; [arcsec]
      c = clips_crisp * imscale_crisp              ; [arcsec]

      if keyword_set(verbose) then begin
        print
        print, 'Old CRISP parameters:'
        print, 'Subfield spacing in arcsec : ', t
        print, 'Clips in arcsec : ', c
        print, 'Clips/spacing : ', c/t
        print
      endif
      
      ;; Calculate for current instrument
      clips = round(c / image_scale)
      tiles = round(dims[0] * image_scale / t) ; Use X dimension? Or max?
      
    end

    'final_size' : begin
 
      ;; One idea is to decide what the final tile spacing should be,
      ;; like 1", which means the final number of tiles should be
      ;; equal to the size of the FOV in arcsec. Then recursively make
      ;; every tile size half of the one that comes after, until we
      ;; either get to, say, two tiles. Or decide that the first tile
      ;; size should not be larger than some size, like 10". And the
      ;; clips could just be something like 1/5 of the tile sizes.

      dims_arcsec = dims * image_scale
      
      if keyword_set(temporal) then begin
        ;; This is for time-series destretching, where the images are not
        ;; supposed to be identical.
        final_spacing = 2.5     ; [arcsec]
      endif else begin
        final_spacing = 1.      ; [arcsec]
      endelse
      red_append, t, final_spacing
      for ii = 1, 10 do begin
        if t[-1] gt 10 then break
        red_append, t, t[-1]*1.25^ii
      endfor
      t = reverse(t)

      ;; Calculate for current instrument
      tiles = round(dims[0] * image_scale / t) ; Use X dimension? Or max?
      clips = round(t / 2.)

      c = clips * image_scale   ; [arcsec]

      if keyword_set(verbose) then begin
        print, 'Subfield spacing in arcsec : ', t
        print, 'Clips in arcsec : ', c
        print, 'Clips/spacing : ', c/t
      endif
      
      
    
      
    end
    
  endcase

  

end

instrument = 'CRISP2'
mode = 'final_size'
mode = 'old_crisp'

case instrument of 
  'CRISP' : begin
    image_scale = 0.0592
    dims = 1024 * [1L, 1L]
  end
  'CRISP2' : begin
    image_scale = 0.050
    dims = 2560 * [1L, 1L]
  end
  'CHROMIS' : begin
    image_scale = 0.035
    dims = 2256 * [1L, 1L]
  end
endcase


red_tiles_and_clips, dims, image_scale, /verbose $
                     , clips = clips $
                     , tiles = tiles $
                     , mode = mode

print
print, instrument + ' Non-temporal'
print, 'clips : ['+strjoin(strtrim(clips,2),', ')+']'
print, 'tiles : ['+strjoin(strtrim(tiles,2),', ')+']'
print

red_tiles_and_clips, dims, image_scale, /verbose $
                     , clips = clips $
                     , tiles = tiles $
                     , /temporal  $
                     , mode = mode

print
print, instrument + ' Temporal' 
print, 'clips : ['+strjoin(strtrim(clips,2),', ')+']'
print, 'tiles : ['+strjoin(strtrim(tiles,2),', ')+']'
print

end
