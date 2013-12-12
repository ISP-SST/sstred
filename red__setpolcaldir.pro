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
; :Returns:
; 
; 
; :Params:
; 
;    dir : in, type=string
;   
;      The directory you want to use as polcal_dir.
;   
; 
; :Keywords:
; 
; 
; 
; :History:
; 
;   2013-12-04 : MGL. Modelled on red::setflatdir.pro
; 
; 
;-
pro red::setpolcaldir, dir
  self.polcal_dir = dir+'/'
end
