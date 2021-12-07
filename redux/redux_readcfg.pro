; docformat = 'rst'

;+
; Read a cfg file.
; 
; :Categories:
;
;    MFBD methods.
; 
; 
; :Author:
; 
;    Mats Löfdahl, Institute for Solar Physics
; 
; 
; :Returns:
;
;    The cfg file contents as an ordered hash.
; 
; :Params:
; 
;    filename : in, type=string
; 
;      A cfg file name.
; 
; :History:
; 
;     2020-03-09 : MGL. First version.
; 
;-
function redux_readcfg, filename

  return, redux_cfg2hash(filename)

end
