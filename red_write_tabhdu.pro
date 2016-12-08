; docformat = 'rst'

;+
; Helper procedure procedure for red_writedata. Write a bintable.
;
; Some funniness involved with uint data. FXB* adds the TSCALE,TZERO
; keywords for uints, but FXBWRITE doesn't know what to do with the
; uint data (it doesn't test for type code 12) so you have to scale it
; down yourself.
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;    Lewis Fox, ISP
; 
; 
; :Returns:
; 
; 
; :Params:
; 
;    tabhdu : in
; 
;   
;   
;    filename : in, type=string
; 
;
; 
; :History:
;
;   2016-09-16 : JLF. First version.
; 
;   2016-12-07 : MGL. Split into a separate file.
; 
; 
;-
pro red_write_tabhdu, tabhdu, filename

  extnames = tag_names(tabhdu)

  for i = 0, n_elements(extnames)-1 do begin

    ;; Create the extension table in the file.
    fxbcreate,lun,filename,tabhdu.(i).bdr,ext_no
    
    ;; Get column data keyword names
    keys = tag_names(tabhdu.(i))

    ;; COMMENT and BDR aren't data keywords, everything else is.
    idx = where(keys ne 'COMMENT' and keys ne 'BDR')
    
    ;; Write the data column by column, everything is in row 1
    for j = 0, n_elements(idx)-1 do begin
      data = tabhdu.(i).(idx[j]).val
      ;; fxbwrite can't handle uint data!?!
      if size(data,/type) eq 12 then data = fix(data-32768)
      fxbwrite,lun,data,j+1,1   ; FITS tables are unit indexed
    endfor
    ;; Done with this extension

    fxbfinish,lun
    
  endfor

end

