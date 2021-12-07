; docformat = 'rst'

;+
; Write a cfg file.
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
; 
; :Params:
; 
;    filename : in, optional, type=string
;
;       Where to store the output.
; 
;    hcfg : in, type=orderedhash
; 
;       The cfg information to be converted.
; 
; :History:
; 
;    2020-03-09 : MGL. First version.
; 
;-
pro redux_writecfg, filename, hcfg

  redux_hash2cfg, hcfg, filename = filename
  
end
