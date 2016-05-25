; docformat = 'rst'

;+
; Sort file names so frame numbers are ascending.
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
; 
; :Returns:
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
; :History:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2016-05-25 : MGL. Make it handle files with a .fits extension. 
; 
; 
;-
function red_sortfiles, files
  nt = n_elements(files)
  num = lonarr(nt)
  for ii = 0L, nt -1 do begin
     tmp = strsplit(file_basename(files[ii],'.fits'),'.',/extract)
     num[ii] = long(tmp[n_elements(tmp) - 1])
  endfor
  pos = sort(num)
  res = files[pos]
  return, res
end
