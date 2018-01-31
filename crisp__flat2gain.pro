; docformat = 'rst'
;+
; Make (inverse) gain table from flat field (or sum thereof). 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;    Tomas Hillberg, 2016
; 
; 
; :Returns:
; 
; 
; :Params:
; 
;    flat : in, type="2D array"
;   
;      A flat field image.
;   
; 
; :Keywords:
;
;    preserve : in, optional, type=boolean
;   
;      If set, don't zero borders between Sarnoff taps.
;
;
;    (Other keywords will be forwarded to red_flat2gain)
;
;
; :History:
; 
; 
;-
function crisp::flat2gain, flat, preserve = preserve, _REF_EXTRA = ex

    ; call generic function
    g = red_flat2gain( flat, _EXTRA = ex )   ; will silently ignore unrecognized keywords
    
    ; crisp/sarnoff specific: zero the boundaries between CCD segments, so they can be filled
    ; later with fillpix.
    if(~keyword_set(preserve)) then for ii = 1,7 do g[ii*128,*] = 0.0                   

    return, g
    
end
