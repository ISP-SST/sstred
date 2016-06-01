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
;   2016-06-01 : MGL. Split at underscores as well as dots. Split all
;                strings in one go.
; 
; 
;-
function red_sortfiles, files

  nt = n_elements(files)
  num = lonarr(nt)

  tmp = strsplit(file_basename(files,'.fits'),'._',/extract)
  for ii = 0L, nt-1 do num[ii] = (tmp[ii])[-1]

  return, files[sort(num)]

end
