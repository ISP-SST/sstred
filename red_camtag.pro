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
; :returns:
; 
; 
; :Params:
; 
;    file : 
;   
;   
; 
; :Keywords:
; 
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
; 
;-
function red_camtag, file
  return, (strsplit(file_basename(file),'.',/extract))[0]
end
