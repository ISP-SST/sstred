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
;    file_list : 
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
function red_inum, file_list
  res = strarr(n_elements(file_list))
                                ;
  for i = 0L, n_elements(res) - 1 do begin
     tmp = strsplit(file_list[i],'.',/extract)
     res[i] = tmp[n_elements(tmp)-1]
  endfor
                                ;
  return, res
end
