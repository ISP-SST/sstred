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
;    root : 
;   
;   
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
function red_select_sub, root
  inam = 'red_select_sub : '
  dir = file_search(root + '/*', /test_dir, count = ndir)
  idx = 0L
  if(ndir gt 1) then begin
     print, inam + 'found '+stri(ndir)+' sub-folders:'
     for jj = 0, ndir - 1 do print, stri(jj,ni='(I2)')+' -> '+dir[jj]
     read,idx,prompt = inam+'Please select state ID: '
     dir = dir[idx]
  endif
  print, inam + 'Selected folder -> '+dir

  return, dir
end
