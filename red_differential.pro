; docformat = 'rst'

;+
; Calculate the differential of the input array.
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, Institute for Solar Physics
; 
; 
; :Returns:
;
;      x - [0, x]
; 
; :Params:
; 
;      x : in, type=array 
; 
;        The input array.
; 
; 
; 
; :History:
; 
;      2017-06-08 : MGL. First version.
; 
; 
; 
; 
;-
function red_differential, x

  return, x - [0, x]
  
end
