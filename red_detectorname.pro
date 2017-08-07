; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author: 
;
; 
; 
; :Returns:
; 
; 
; :Params:
; 
;    file : 
;   
;   
; 
; :Keywords:
; 
;    head : in, type=string array
;
;	An input header to check for DETECTOR, instead of calling
;	red_readhead.
; 
; :History:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2016-08-31 : MGL. Rename red_getcamtag to red_detectorname. Get
;                the detector name from the header if possible.
;
;   2016-09-01 : JLF. Added head=head keyword to bypass an infinite loop;
;                red_detectorname and red_readhead were calling each other.
; 
;-
function red_detectorname, file, head=head

  ;; if we've been passed a header already don't try to use red_readhead.
  ;; This bypasses an infinite loop!

  if n_elements(head) ne 0 then $
    status = 0 $
  else head = red_readhead(file, status = status)
  
  ;; First try to get the detector from the headers
  
  if status eq 0 then begin
     detector = fxpar(head, 'DETECTOR')
     if detector ne '' then return, strtrim(detector, 2)
  endif

  ;; Backup: try with the file name
  detector = (strsplit(file_basename(file),'.',/extract))[0]
  return, detector


end
