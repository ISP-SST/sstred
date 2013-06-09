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
;    mm : in, type="fltarr(4,4,n,m)"
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
function red_invert_mmatrix, mm
  inam = 'red_invert_mmatrix : '
  dim = size(mm, /dim)
  imm = fltarr(16,dim[1], dim[2])
                                ;
  bb = string(13B)
  ntot = 100. / (dim[2] - 1.0)
                                ;
  for jj = 0L, dim[2] - 1 do begin
     for ii = 0L, dim[1] - 1 do begin
        imm[*,ii,jj] = reform(invert(reform(mm[0:15,ii,jj], [4,4])),16)
     endfor
     print, bb, inam + 'inverting modulation matrix over the FOV -> ', $
            ntot * jj,'%', format = '(A, A, F5.1, A, $)'
  endfor
  print, ' '
                                ;
  return, imm
end
