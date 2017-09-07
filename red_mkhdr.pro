; docformat = 'rst'

;+
; Wrapper, make FITS header they way we want them.
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, ISP
; 
; 
; :Returns:
; 
; 
; :Params:
; 
;    header : out
;
;      See documentation for mkhdr.pro.
;
;    im : in
;
;      See documentation for mkhdr.pro.
;
;    naxisx : in
;
;      See documentation for mkhdr.pro.
; 
; :Keywords:
; 
;    image : in
;
;      See documentation for mkhdr.pro.
;   
;    extend : in
;
;      See documentation for mkhdr.pro.
; 
; 
; :History:
;
;   2017-03-17 : MGL. First version. 
;
;   2017-03-22 : MGL. Add TIMESYS. 
;
;   2017-05-30 : MGL. Use red_fitsaddpar. 
;
;   2017-09-06 : MGL. Use new features in red_fitsaddpar.
;
;   2017-09-07 : MGL. Changed red_fitsaddpar --> red_addfitskeyword. 
; 
;-
pro red_mkhdr, header, im, naxisx, _ref_extra = extra 

  if n_elements(naxisx) gt 0 then begin
    mkhdr, header, im, naxisx, _strict_extra = extra 
  endif else begin
    mkhdr, header, im, _strict_extra = extra 
  endelse

  red_fitsaddkeyword, header, 'TIMESYS', 'UTC', before = 'DATE'

  ;; Add the OGIP warning near the end of the header so it doesn't
  ;; appear near the first added long keyword. Do it "by hand" rather
  ;; than calling FXADDPAR_CONTWARN because 1) that routine may put
  ;; the comment lines in reverse order and 2) it might not be
  ;; compiled as it is defined in fxaddpar.pro as a helper routine.
  
  red_fitsaddkeyword, header, 'LONGSTRN', 'OGIP 1.0', 'The OGIP long string convention may be used.' $
                      , before = 'COMMENT'
  red_fitsaddkeyword, header, '', '', before = 'LONGSTRN' 
  red_fitsaddkeyword, header, ['COMMENT', '', ''], ["This FITS file may contain long string keyword values that are continued over multiple keywords.  This convention uses the '&' character at the end of a string which is then continued on subsequent keywords whose name = 'CONTINUE'.", '', ''] $
              , after = 'LONGSTRN'
  
end

red_mkhdr,header,0 
print, header[where(strtrim(header, 2) ne '')], format = '(a0)'

end
