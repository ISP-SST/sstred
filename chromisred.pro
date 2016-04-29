; docformat = 'rst'

;+
; Class CHROMIS
;
; :Author:
; 
;    Tomas Hillberg
;
; :Params:
; 
;   filename : in, optional, type=string, default="config.txt"
;
; :History:
;
;   2016-04-29 : THI. Split class RED into a base-class (instrument
;                independent parts) and derived classes (CRISP/CHROMIS). 
;
;-
function chromisred, filename

  return, obj_new('chromis', filename)
  
end
