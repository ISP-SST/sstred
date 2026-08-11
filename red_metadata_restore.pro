; docformat = 'rst'

;+
; Restore previously stored metadata, optionally adding it to an
; existing FITS header.
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
; :Params:
; 
;    head : in, out, type=strarr
;
;      The FITS header.
; 
; :Keywords:
;   
;   anchor : in, optional, type=string
;
;     See documentation of red_fitsaddkeyword.
; 
;   fname : in, optional, type=string, default='info/metadata.fits'
;
;     Name of file from which to read the metadata.
;
; :History:
; 
;    2017-03-13 : MGL. First version.
; 
;    2017-11-14 : MGL. Check that the file exists.
; 
;    2018-05-25 : MGL. New keyword anchor.
; 
; 
;-
pro red_metadata_restore, head, fname = fname, anchor = anchor

  if n_elements(fname) eq 0 then fname = 'info/metadata.fits'
  if n_elements(anchor) eq 0 then anchor = 'DATE'
  
  if file_test(fname) then begin
    
    metahead = headfits(fname)

    if n_elements(head) eq 0 then begin
      head = metahead
      return
    endif
    
    ;; Add metadata info to the provided header
    Nlines = n_elements(metahead)
    for iline = 0, Nlines-1 do begin

      ;; Skip lines that don't have keywords with values.
      if strmid(metahead[iline], 8, 1) ne '=' then continue
      
      keyword = (strsplit(metahead[iline], ' =', /extract))[0]

      ;; Skip some administrative keywords
      if keyword eq 'SIMPLE' $
         or keyword eq 'BITPIX' $
         or keyword eq 'NAXIS' $
         or keyword eq 'EXTEND' $
         or keyword eq 'DATE' $
         or keyword eq 'LONGSTRN' $
      then continue

      red_fitscopykeyword, anchor = anchor, head, keyword, metahead
      
    endfor                      ; iline
    
  endif  

end
