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
;   im_in : 
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
function red_fillnan,im_in
                                ;
  im = im_in
  pos = where(~finite(im_in), count, complement = pos1)
  if(count gt 0 AND (count lt n_elements(im_in))) then im[pos] = median(im[pos1])
  if(count eq n_elements(im_in)) then im[*] = 0.0
  
                                ;
  return, im
end
