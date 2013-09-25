; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
; 
; 
; 
; :Returns:
; 
; 
; :Params:
; 
; 
;    ny : 
;
;
;    maxs : 
; 
; 
; :History:
; 
;    2013-09-12 : MGL. Brought into red_ namespace.
; 
; 
; 
;-
function red_findthedim, ny, maxs
  res = 0L
  for ii = 1L, maxs do begin
     if((ny / ii) * ii EQ ny) then res = ii
  endfor
  print, 'findthedim: using ny1 = '+string(res, format='(I5)')
  return, res
end
