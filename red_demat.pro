; docformat = 'rst'

;+
; 
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
; 
; :Params:
; 
;    g : 
;   
;   
;   
;    flat : 
;   
;   
;   
;    imm : 
;   
;   
;   
; :History:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2018-11-13 : MGL. Don't change the input.
; 
; 
;-
function red_demat, g, flat, imm

  lim = 0.001

  idx = where((g[*,*,0] gt lim) and $
              (g[*,*,1] gt lim) and $
              (g[*,*,2] gt lim) and $
              (g[*,*,3] gt lim), count)

  if count ge 1 then me = mean(flat[idx]) else me = mean(flat)
  
  gf0 = red_fillnan(me/red_fillpix(red_mask_ccd_tabs(g[*,*,0] * flat)))
  gf1 = red_fillnan(me/red_fillpix(red_mask_ccd_tabs(g[*,*,1] * flat)))
  gf2 = red_fillnan(me/red_fillpix(red_mask_ccd_tabs(g[*,*,2] * flat)))
  gf3 = red_fillnan(me/red_fillpix(red_mask_ccd_tabs(g[*,*,3] * flat)))

  imm_out = imm
  
  for ii = 0, 3 do begin
    ;; Transmitted
    imm_out[ii*4+0,*,*] *= gf0 
    imm_out[ii*4+1,*,*] *= gf1
    imm_out[ii*4+2,*,*] *= gf2
    imm_out[ii*4+3,*,*] *= gf3
  endfor
  
  return, imm_out

end
