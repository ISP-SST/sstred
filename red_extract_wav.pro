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
;    state : 
;   
;   
;   
; 
; :Keywords:
; 
;    lc  : 
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
function red_extract_wav, state, lc = lc
  nt = n_elements(state)
  res = dblarr(nt)
  lc = strarr(nt)

  for ii = 0L, nt - 1 do begin
     tmp = strsplit(state[ii],'.',/extract)
     lc[ii] = tmp[2]

     tmp = tmp[1]
     if(strpos(tmp,'_') eq -1) then res[ii] = double(0.0d0) else begin
        tmp = strsplit(tmp[1], '_', /extract)
        res[ii] = double(tmp[0]) + double(tmp[1]) * 0.001d0
     endelse

  endfor
  return, res
end
