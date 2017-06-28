; docformat = 'rst'

;+
; Store metadata in a fits header-only file.
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
;    metadata : in, type=structarr
;
;      An array of structs containing metadata in the form {keyword,
;      value, comment}. 
;   
; :Keywords:
;   
;    fname : in, optional, type=string, default='info/metadata.fits'
;   
;      The name of the file in which to store the metadata. 
;      
;    force : in, optional, type=boolean
;
;      Set this to overwrite existing metadata with the specified
;      keyword. 
; 
; :History:
; 
;    2017-03-06 : MGL. First version.
; 
;    2017-06-28 : MGL. Use red_fitsaddpar.
;
;-
pro red_metadata_store, metadata $
                        , fname = fname $
                        , force = force

  Nkeywords = n_elements(metadata)

  if n_elements(fname) eq 0 then fname = 'info/metadata.fits'
  file_mkdir, file_dirname(fname)

  if file_test(fname) then begin
    ;; Read existing metadata fits header
    head = headfits(fname)
  endif else begin
    ;; Create an empty header
    mkhdr, head, 0
  endelse

  header_changed = 0

  anchor = 'DATE'
  ;; Loop through the keywords
  for ikeyword = 0, Nkeywords-1 do begin
    
    ;; Check if the keyword is already there
    aa = fxpar(head, metadata[ikeyword].keyword, count = N)
    
    ;; Add it if it's not or if force flag is set
    if keyword_set(force) or N eq 0 then begin
      ;; Add a leading space to comment if needed
      comment = metadata[ikeyword].comment
      if strlen(comment) gt 0 then begin
        if strmid(comment, 0, 1) ne ' ' then comment = ' '+comment
      endif
      ;; Add the keyword to the header
      red_fitsaddpar, anchor = anchor, head $
                , metadata[ikeyword].keyword $
                , metadata[ikeyword].value $
                , comment
      header_changed = 1
    endif

  endfor                        ; ikeyword

  ;; Write the updated header to disk if it has changed.
  if header_changed then fxwrite, fname, head 
  
end
