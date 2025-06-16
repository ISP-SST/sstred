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
;     oldhdr1 : in, type=strarr
; 
;        The source FITS header. 
; 
;     oldhdr2 : in, optional,type=strarr
; 
;        A second source FITS header. 
; 
; 
; :Keywords:
; 
;     _ref_extra : in, optional
;
;        Any keywords are used when calling red_fitsaddkeyword.
; 
; :History:
; 
;    2017-12-06 : MGL. First version.
; 
;-
pro red_fitscopykeyword, hdr, name, oldhdr1, oldhdr2, _ref_extra = extra

  value1 = red_fitsgetkeyword(oldhdr1, name, comment = comment1, count = count1)
  if n_elements(oldhdr2) ne 0 then begin
    value2 = red_fitsgetkeyword(oldhdr2, name, comment = comment2, count = count2)
    if count2 eq 1 then $
       red_fitsaddkeyword, hdr, name,  strtrim(value1, 2) $
                           + ',' + strtrim(value2, 2), comment1 $
                           , _strict_extra = extra
  endif else begin
    if count1 eq 1 then red_fitsaddkeyword, hdr, name, value1, comment1 $
                                            , _strict_extra = extra
  endelse
  

end

