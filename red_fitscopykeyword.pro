; docformat = 'rst'

;+
; Copy a keyword from one FITS header to another. 
;
; Do not use for record-valued keywords or SOLARNET variable-keywords.
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
; :Params:
;
;     hdr : in, type=strarr
; 
;        The target FITS header. 
; 
;     name : in, type=strarr
;     
;        The name of the parameter/keyword to be copied. See
;        documentation for fxaddpar.
;
;     oldhdr : in, type=strarr
; 
;        The source FITS header. 
; 
; 
; :Keywords:
; 
;     Any keywords are used when calling red_fitsaddkeyword.
; 
; :History:
; 
;    2017-12-06 : MGL. First version.
; 
;-
pro red_fitscopykeyword, hdr, name, oldhdr, _ref_extra = extra

  value = red_fitsgetkeyword(oldhdr, name, comment = comment, count = count)
  
  if count eq 1 then red_fitsaddkeyword, hdr, name, value, comment, _strict_extra = extra

end

