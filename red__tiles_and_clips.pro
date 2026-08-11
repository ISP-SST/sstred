; docformat = 'rst'

;+
; Wrapper around red_stretch_tiles.
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
; :Params:
; 
;   dims : in, type=array
; 
;     Size of the FOV in pixels, [Nx,Ny].
; 
;   pref : in, type=string
; 
;     The prefilter for which to calculate the image scale.
; 
; :Keywords:
; 
;     _ref_extra : in, out, optional
;   
;        Keywords are passed on to red_stretch_tiles.pro.
; 
; 
; :History:
; 
;     2026-03-06 : MGL. First version.
; 
;-
pro red::tiles_and_clips, dims, pref, _ref_extra = extra

  ;; Call red_stretch_tiles with the relevant image_scale
  image_scale = self -> imagescale, pref
  red_tiles_and_clips, dims, image_scale, _strict_extra = extra 
  
end

