; docformat = 'rst'

;+
; Sort file names so frame numbers are ascending.
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
;    The sorted array of file names.
; 
; 
; :Params:
; 
;    files : in, type=strarr
;   
;      An array of file names.
;   
; 
; :Keywords:
; 
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
; 
;-
function red_sortfiles, files
  nt = n_elements(files)
  num = lonarr(nt)
  for ii = 0L, nt -1 do begin
     tmp = strsplit(files[ii],'.',/extract)
     num[ii] = long(tmp[n_elements(tmp) - 1])
  endfor
  pos = sort(num)
  res = files[pos]
  return, res
end
