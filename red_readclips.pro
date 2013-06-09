; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :author:
; 
; 
; 
; 
; :returns:
; 
; 
; :Params:
; 
;    file : 
;   
;   
;   
; 
; :Keywords:
; 
;    tc  : 
;   
;   
;   
;    rc  : 
;   
;   
;   
;    wc  : 
;   
;   
;   
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
; 
;-
pro red_readclips, file, tc = tc, rc = rc, wc = wc
  inam = 'red_readclips : '
  if(~file_test(file)) then begin
     print, inam + 'ERROR -> file not found '+file
     STOP
  endif
  openr, lun, file, /get_lun
  tc = ' '
  rc = ' '
  wc = ' '
  readf, lun, wc
  readf, lun, tc
  readf, lun, rc
  free_lun, lun 
  return
end
