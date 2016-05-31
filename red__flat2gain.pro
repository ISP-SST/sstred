; docformat = 'rst'
;+
;    Make (inverse) gain table from flat field (or sum thereof).
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :author:
; 
;    Tomas Hillberg, 2016
; 
; 
; :returns:
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
;    All (accepted) keywords will be forwarded to red_flat2gain
;
;
; :history:
; 
; 
;-
function red::flat2gain, flat, _REF_EXTRA = ex
    
    return, red_flat2gain( flat, _EXTRA = ex )   ; will silently ignore unrecognized keywords
    
end
