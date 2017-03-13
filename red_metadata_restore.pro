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
;   fname : in, optional, type=string, default='info/metadata.fits'
;   
;   
; 
; 
; :History:
; 
;    2017-03-13 : MGL. First version.
; 
; 
;-
pro red_metadata_restore, head, fname = fname

  if n_elements(fname) eq 0 then fname = 'info/metadata.fits' 

  metahead = headfits(fname)

  if n_elements(head) eq 0 then begin
    head = metahead
    return
  endif
  
  ;; Add metadata info to the provided header
  Nlines = n_elements(metahead)
  for iline = 0, Nlines-1 do begin
    
    keyword = (strsplit(metahead[iline], ' =', /extract))[0]
    if keyword eq 'COMMENT' then continue

    oldvalue = sxpar(head, keyword, comment = comment, count = count)
    if count eq 0 then begin
      metavalue = sxpar(metahead, keyword, comment = metacomment)
      sxaddpar, head, keyword, metavalue, metacomment
    endif
    
  endfor                        ; iline
  
  

end
