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
;    cub : 
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
function red_get_cross, cub
                                ;
  inam = 'get_red_cross : '
  stkname = ['I', 'Q', 'U', 'V']
                                ;
  crt = fltarr(4)
  crt[0] = median(cub[*,*,0])
                                ;
  for ii = 1, 3 do begin
     tmp = red_fillnan(reform(cub[*,*,ii]) / cub[*,*,0])
     me = median(tmp)
     mask = abs(tmp-me) LE (2.0 * stdev(tmp))
     idx = where(mask, count)
                                ;
     if(count ge 1) then begin
        crt[ii] = mean(tmp[idx]) 
        if(~finite(crt[ii])) then crt[ii] = median(tmp[idx])
     endif else begin
        crt[ii] = mean(tmp)
        if(~finite(crt[ii])) then crt[ii] = median(tmp[idx])
        print, inam + 'WARNING, something could be wrong with I -> '+stkname[ii]+' = ' + crt[ii]
     endelse
  endfor
  return, crt
end
