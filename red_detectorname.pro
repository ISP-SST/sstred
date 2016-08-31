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
; 
; 
; :History:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2016-08-31 : MGL. Rename red_getcamtag to red_detectorname. Get
;                the detector name from the header if possible.
; 
;-
function red_detectorname, file

  ;; First try to get the detector from the headers
  head = red_readhead(file, status = status)

  if status eq 0 then begin
     detector = sxpar(head, 'DETECTOR')
     if detector ne '' then return, detector
  endif

  ;; Backup: try with the file name
  detector = (strsplit(file_basename(file),'.',/extract))[0]
  return, detector


end
