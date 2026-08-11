; docformat = 'rst'

;+
; 
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
;    TRUE for odd integers, FALSE for even integers.
;
; :Params:
; 
;    i : in, type="integer or integer array"
;
;      Input data, the oddity of which is needed.
; 
; 
; :History:
; 
;-
function red_odd, i

  return, abs(i mod 2)

end

