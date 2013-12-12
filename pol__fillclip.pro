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
;    x0 : 
;   
;   
;   
;    x1 : 
;   
;   
;   
;    y0 : 
;   
;   
;   
;    y1 : 
;   
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
pro pol::fillclip, x0, x1, y0, y1
  self.x0 = x0
  self.x1 = x1
  self.y0 = y0
  self.y1 = y1
  return
end
