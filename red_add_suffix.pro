; docformat = 'rst'

;+
; Parses the filename and adds a specified suffix before the
; extension.
;
; This function parses the filename and separates it into the name itself and
; the extension.  If there is no extension, only the name is taken.  It adds a
; specified suffix inbetween or nothing if the suffix is not given.
; 
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Andrii V. Sukhorukov, Institute for Solar Physics
; 
; 
; :Returns:
; 
; 
; :Params:
; 
;   filename : in, type=string
; 
;      The file name to which the suffix is to be added.
; 
; 
; :Keywords:
; 
;   suffix : in, optional, type=string
;   
;      The suffix.
; 
; 
; :History:
; 
;    2017-08-23 : AVS. First version.
; 
; 
; 
; 
;-
function red_add_suffix, filename, suffix = suffix

  if ( n_elements( suffix ) eq 0 ) the return
  if suffix eq '' then return
  
  filenamelen = strlen( filename )
  if ( filenamelen eq 0 ) then begin
    message, 'filename must not be empty.'
    retall
  endif
  dotpos = strpos( filename, '.', /reverse_search )
  if ( dotpos ne -1 ) then begin
    name      = strmid( filename, 0,      dotpos               )
    extension = strmid( filename, dotpos, filenamelen - dotpos )
  endif else begin
    name      = filename
    extension = ''
  endelse
  
  return, name + suffix + extension
  
end
