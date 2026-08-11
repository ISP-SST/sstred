; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :author:
; 
; 
; 
; 
; :returns:
; 
; 
; :Params:
; 
;    scan : 
;   
;   
;   
; 
; :Keywords:
; 
;    hscan  : 
;   
;   
;   
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
; 
;-
function red_decode_scan, scan, hscan = hscan
  hscan = strmid(scan,0,1)
  return, '0'+strmid(scan,1)
end
