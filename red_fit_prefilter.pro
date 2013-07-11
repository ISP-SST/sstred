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
;    pp : 
;   
;   
;   
; 
; :Keywords:
; 
;    xs  : 
;   
;   
;   
;    ys  : 
;   
;   
;   
;    dpr  : 
;   
;   
;   
;    mm  : 
;   
;   
;   
;    pref  : 
;   
;   
;   
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2013-07-11 : MGL. Use red_intepf, not intepf.
; 
; 
; 
;-
function red_fit_prefilter, pp, xs = xs, ys = ys, dpr = dpr, mm = mm, pref = pref
  device, decompose=0
  
  p0 = pp[0]
  p1 = pp[1]
  p2 = pp[2]
  p3 = pp[3]
  p4 = pp[4]
  p5 = pp[5]
  p6 = pp[6]
  p7 = pp[7]
                                ;
;  wav = mm.wav * p7
  xs *= p7
  ys1 = red_intepf(xs, ys, mm.wav - p4, /linear)
  pref =  1.d0 / (1.d0+((2.d0*(mm.wav - p1) / p2)^2.)^p3) * (1.d0 + (mm.wav * p5) + (mm.wav^3. * p6))
                                ;
  res = p0 * pref * ys1 
                                ;
  loadct,1,/silent
  plot, mm.wav, res, psym = -4, yrange=[0,2] * median(mm.yl1), /ystyle
  oplot, mm.wav, mm.yl1, /line, psym=-1
  oplot, mm.wav, pref*median(mm.yl1), color = 175, thick = 1, psym= -1
  legend, ['Model', 'Observed', 'Prefilter'], line = [0, 1, 0], psym = [-4, -1, -1], color = [255, 255, 175], $
          /bottom, /right, charsize=2
  loadct,0,/silent
                                ;
  return, res - mm.yl1
end
