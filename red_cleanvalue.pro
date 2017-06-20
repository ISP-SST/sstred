; docformat = 'rst'

;+
; Clean a keyval "value" string, strip blanks and quotes from
; beginning and end, detect array representation.
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
;     The cleaned value.
; 
; :Params:
; 
;     x : in, type=string
; 
;       The string to be stripped.
; 
; 
; :Keywords:
; 
;   
;   
;   
; 
; 
; :History:
; 
;    2017-06-20 : MGL. First version.
; 
; 
; 
; 
;-
function red_cleanvalue, x, status = status

  ;; Trim whitespace from begining and end of string
  xx = strtrim(x, 2)

  ;; Is this actually an array?
  if strmid(xx,0,1) eq '[' and strmid(xx,strlen(xx)-1,1) eq ']' then begin
    dum = execute('tmp='+xx)
    status = 0
    return, tmp
  endif
  
  ;; Split at
  ;; quotes. Single and double quotes are both OK.
  xx = strsplit(xx, "'" + '"', /extr)

  ;; If the only quotes are at the beginning and end of the trimmed
  ;; string (or if there are no quotes), then the splitting will
  ;; result in an array of length 1.
  if n_elements(xx) eq 1 then begin
    ;; Success
    status = 0
    return, xx[0]
  end

  ;; Failed, result undefined. Return empty string.
  status = -1
  return, ''

end
