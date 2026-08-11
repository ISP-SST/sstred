; docformat = 'rst'

;+
; Class SLITJAW.
;
; :Author:
; 
;    Mats Löfdahl, ISP
;
; :Params:
; 
;   filename : in, optional, type=string, default="config.txt"
;
;      The name of the config file.
; 
; :Keywords:
; 
;    develop : in, optional, type=boolean
; 
;       Run in developer mode.
;
; :History:
;
;   2017-03-13 : MGL. First version.
;
;-
function slitjawred, filename, develop = develop

  return, obj_new('slitjaw', filename, develop = develop)
  
end
