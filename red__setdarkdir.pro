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
;     Mats Löfdahl, ISP
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
;   2014-03-21 : MGL. New method, based on red::setflatdir.
; 
;-
pro red::setdarkdir, dir
  *self.dark_dir = dir
  FOR i = 0, n_elements(*self.dark_dir)-1 DO $
     IF strmid((*self.dark_dir)[i], strlen((*self.dark_dir)[i])-1) NE '/' THEN $
        (*self.dark_dir)[i] += '/'
END

