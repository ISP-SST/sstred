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
; :History:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2013-07-11 : MGL. Use red_intepf, not intepf.
; 
;   2014-01-07 : PS. red_legend -> al_legend (from astron). Color
;                handling.
;
;   2016-09-26 : MGL. Improve plotting: cgplot, more colors, works for
;                a single wavelength point.
;
;-
function red_fit_prefilter, pp, xs = xs, ys = ys, dpr = dpr, mm = mm, pref = pref, w=w

  p0 = pp[0]
  p1 = pp[1]
  p2 = pp[2]
  p3 = pp[3]
  p4 = pp[4]
  p5 = pp[5]
  p6 = pp[6]
  p7 = pp[7]

;  wav = mm.wav * p7
  xs *= p7
  ys1 = red_intepf(xs, ys, mm.wav - p4, /linear)
  pref =  1.d0 / (1.d0+((2.d0*(mm.wav - p1) / p2)^2.)^p3) * (1.d0 + (mm.wav * p5) + (mm.wav^3 * p6))

  res = p0 * pref * ys1 

  cgplot, mm.wav, res, psym = 4 $
          , xrange = [min(mm.wav)-1, max(mm.wav)+1], yrange=[0,2] * median(mm.yl1) $
          , /ystyle, symsize = 2, thick = 2
  cgplot, /over, mm.wav, mm.yl1, line = 0, psym=-1, color = 'red', symsize = 2, thick = 2
  cgplot, /over, mm.wav, pref*median(mm.yl1), color = 'blue', thick = 2, line = 2

  ;al_legend, ['Model', 'Observed', 'Prefilter'], line = [0, 1, 0], psym = [-4, -1, -1], $
  ;           color = [255, 255, 175], /bottom, /right, charsize=2
  
  return, (res - mm.yl1)*w

end
