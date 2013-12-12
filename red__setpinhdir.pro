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
;    dir : in, type=string
;   
;      The directory you want to use as pinh_dir.
;   
; 
; :Keywords:
; 
; 
; 
; :history:
; 
;   2013-07-09 : MGL. Modelled on red::setflatdir.pro
; 
; 
;-
pro red::setpinhdir, dir
  self.pinh_dir = dir+'/'
end
