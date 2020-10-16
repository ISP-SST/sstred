; docformat = 'rst'

;+
; Copy a binary extension from one FITS file to another.
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
;    file_orig, in, type=string
; 
;       The file from which to copy the extension.
; 
;    file_target, in, type=string
; 
;       The file to which the extension is to be copied.
; 
;    extension, in, type="string or integer"
; 
;       The name or number of the extension. Note that when file_orig
;       has multiple extensions with the same name, specifying the
;       number gives more control.
; 
; 
; :History:
; 
; 
;-
pro red_fits_copybinext, file_orig, file_target, extension
  
  ;; Open the binary extension in the original file and transfer
  ;; header and data to an extension in the target file.

  ;; Header
  fxbopen,   olun, file_orig, extension, bdr
  fxbcreate, tlun, file_target, bdr ; Create the extension in the file
  
  ;; Find the columns in the extension, transfer them to the new file
  fxbfind, olun, 'TTYPE', column_numbers, column_names, Ncolumns
  for icolumn = 0, Ncolumns-1 do begin
    fxbread,  olun, column_data, column_numbers[icolumn], 1
    fxbwrite, tlun, column_data, column_numbers[icolumn], 1
  endfor                        ; icolumn

  ;; Close the binary extensions
  fxbfinish, tlun             
  fxbclose,  olun

end
