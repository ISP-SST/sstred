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
;    st : 
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
function red_flagchange, st, remove = remove
  nt = n_elements(st)
  star = bytarr(nt)
  if(n_elements(remove) eq 0) then return, temporary(star)

                                ;
  os = st[0]
  for ii = 0L, nt - 1 do begin
     if(st[ii] ne os) then begin
        for jj = 0, remove-1 do star[ii+jj] = 1B
        os = st[ii]
     endif
  endfor
                                ;
  return, star
end
