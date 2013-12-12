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
;    dir : 
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
;   2013-12-10 : PS  adapt for multiple flat_dir
; 
;-
pro red::setflatdir, dir
*self.flat_dir = dir
FOR i = 0, n_elements(*self.flat_dir)-1 DO $
  IF strmid((*self.flat_dir)[i], strlen((*self.flat_dir)[i])-1) NE '/' THEN $
    (*self.flat_dir)[i] += '/'
END

