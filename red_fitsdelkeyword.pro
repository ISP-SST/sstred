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
;   hierarch_only : in, optional, type=boolean
;   
;     Remove only the HIERARCH version of the keyword.
; 
; 
; :History:
; 
;   2017-09-08 : MGL. First version.
; 
;   2018-04-06 : MGL. Delete also HIERARCH keywords.
; 
;   2018-04-06 : MGL. New keyword hierarch_only.
; 
;-
pro red_fitsdelkeyword, hdr, name, hierarch_only = hierarch_only

  ;; There is no fxdelpar so we have to use sxdelpar. However,
  ;; sxdelpar is not aware of the OGIP LONGSTRN convention, so if the
  ;; keyword is continued on the next line, the CONTINUE line is not
  ;; removed. So we first make sure the keyword is short, then we
  ;; remove it.

  if ~keyword_set(hierarch_only) then begin
    fxaddpar, hdr, name, 'del', 'del'
    sxdelpar, hdr, name
  endif
  
  ;; The above commands take care of ordinary keywords, as well as
  ;; record-valued keywords. But not HIERARCH keywords! So we rewrite
  ;; such lines as "normal" (but protected) keywords and then  remove
  ;; them. 

  hindx = where(strmatch(hdr, 'HIERARCH ' + name + ' *'), Nmatch)
  for imatch = 0, Nmatch-1 do $
     hdr[hindx[imatch]] $
     = "+DEL+    = 'DEL'                                                                "
  sxdelpar, hdr, "+DEL+"
  
end

mkhdr, h, 0

;; Add a long string keyword:
fxaddpar, h, 'TEST', 'A very long string that will need to be continued on the next line using the OGIP LONGSTRN convention.'
print
print, h

;; Remove the keyword with sxdelpar, notice that the CONTINUE line is
;; still there:
sxdelpar, h, 'TEST'
print
print, h



mkhdr, h, 0

;; Try again, but now delete with red_fitsdelkeyword
fxaddpar, h, 'TEST', 'A very long string that will need to be continued on the next line using the OGIP LONGSTRN convention.'
red_fitsdelkeyword, h, 'TEST'
print
print, h




end
