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
;   pro pol::unloadimages : 
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
pro pol::unloadimages
   
  nt = n_elements(self.tfiles)
   
  ;; Free the memory that is pointed at!
   
  ptr_free, self.timg
  ptr_free, self.rimg
   
  return
end
