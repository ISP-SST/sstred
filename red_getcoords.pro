; docformat = 'rst'

;+
; Calculates the shifts relative to the new reference
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
;    var : 
;   
;   
;   
;    pos : 
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
;   2013-07-24 : MGL. Use red_show rather than show. Changed dimension
;                check from dim1[0] to dim1[2] when deciding whether
;                to call align.
;
; 
;-
function red_getcoords, var, pos
                                
  dim1=size(var,/dimension)
  ref=reform(var[*,*, pos])
  red_show,ref,/nowin
  if(n_elements(dim1) eq 2) then dim1 = [dim1, 1L] 
  dum=indgen(dim1[2])
  
  pos1=where(dum ne total(pos))


  res=dblarr(dim1)
  offs=dblarr(2,dim1[2])
  if dim1[2] eq 1 then begin
     offs[*,0]=[0.d0,0.d0]
  endif else begin
     for l=0,n_elements(pos1)-1 do begin
        t=reform(var[*,*,pos1[l]])
        offs[*,pos1[l]]=align(ref,t)
     endfor
  endelse
  return,offs
end
