; docformat = 'rst'

;+
; Delete a keyword from a FITS header.
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
; 
; :Params:
; 
; 
; 
; 
; 
; 
; :Keywords:
; 
;   
;   
;   
; 
; 
; :History:
; 
;   2017-09-08 : MGL. First version.
; 
; 
; 
; 
;-
pro red_fitsdelkeyword, hdr, name

  ;; There is no fxdelpar so we have to use scdelpar. However,
  ;; sxdelpar is not aware of the OGIP LONGSTRN convention, so if the
  ;; keyword is continued on the next line, the CONTINUE line is not
  ;; removed. So we first make sure the keyword is short, then we
  ;; remove it.

  fxaddpar, hdr, name, 'del', 'del'
  sxdelpar, hdr, name
  
end
