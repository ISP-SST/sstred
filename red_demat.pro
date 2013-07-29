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
function red_demat, g, flat, imm
  lim = 0.001
  idx = where((g[*,*,0] gt lim) AND $
              (g[*,*,1] gt lim) AND $
              (g[*,*,2] gt lim) AND $
              (g[*,*,3] gt lim), count)
  if(count ge 1) then me = mean(flat[idx]) else me = mean(flat)
  


  gf0 = red_fillnan(me/red_fillpix(red_mask_ccd_tabs(g[*,*,0] * flat)))
  gf1 = red_fillnan(me/red_fillpix(red_mask_ccd_tabs(g[*,*,1] * flat)))
  gf2 = red_fillnan(me/red_fillpix(red_mask_ccd_tabs(g[*,*,2] * flat)))
  gf3 = red_fillnan(me/red_fillpix(red_mask_ccd_tabs(g[*,*,3] * flat)))

  for ii = 0, 3 do begin
                                ; Transmitted
     imm[ii*4+0,*,*] *= gf0 
     imm[ii*4+1,*,*] *= gf1
     imm[ii*4+2,*,*] *= gf2
     imm[ii*4+3,*,*] *= gf3
  endfor
  
  return, imm
end
