; docformat = 'rst'

;+
; Find SIMX and SIMY grid coordinates for pinhole array image. Assume
; the grid is fairly well aligned to the image X and Y directions.
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :author:
; 
;    From Pit's setup_ph.pro
; 
; 
; :returns:
; 
; 
; :Params:
; 
;    pinholeimage : 
;   
;   
;   
;    simx : 
;   
;   
;   
;    simy : 
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
pro red_findpinholegrid, pinholeimage, simx, simy

  ;; Each pinhole gets a unique number
  mask = red_separate_mask(pinholeimage gt .05*max(pinholeimage))
  ;; # pinholes found
  nph = max(mask)

  ;; Compute PH positions
  cc = fltarr(2, nph)
  FOR i=0, nph-1 DO cc(*, i) = red_com(mask EQ i+1)
  cx = reform(cc(0, *))
  cy = reform(cc(1, *))

  ;; Sort values. PHs should be aligned hor/vert, so this will give a
  ;; clear step shape
  cx = cx(sort(cx))
  cy = cy(sort(cy))

  ;; Locate the steps and average the values of each step
  dcx = cx(1:*)-cx
  scx = [-1, where(dcx GT avg(dcx), nx), nph-1]
  simx = intarr(nx+1)
  FOR i=1, nx+1 DO simx(i-1) = round(avg(cx(scx(i-1)+1:scx(i))))

  ;; Now the same for y
  dcy = cy(1:*)-cy
  scy = [-1, where(dcy GT avg(dcy), ny), nph-1]
  simy = intarr(ny+1)
  FOR i=1, ny+1 DO simy(i-1) = round(avg(cy(scy(i-1)+1:scy(i))))

end 
