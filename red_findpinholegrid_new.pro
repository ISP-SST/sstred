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
;    pinholeimage : in, type=array
;   
;   
;   
;    x : out
;   
;   
;   
;    y : out
;   
;   
;   
; 
; :Keywords:
; 
;    thres : in, optional, type=float
; 
;       Threshold.
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
;
;   2013-09-10 : MGL. Return coordinates for each pinhole rather than
;                for a regular grid. Remove pinholes that are closer
;                to an edge than half the grid spacing.
;
;   2014-01-23 : MGL. New keyword thres.
; 
; 
;-
pro red_findpinholegrid_new, pinholeimage, x, y, thres = thres

  if n_elements(thres) eq 0 then thres = 0.05

  ;; Each pinhole gets a unique number
  mask = red_separate_mask(pinholeimage gt thres*max(pinholeimage))
  ;; # pinholes found
  nph = max(mask)
stop
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

  ;; Grid spacing:
  dx = median(deriv(simx))
  dy = median(deriv(simy))
  
  ;; Pick only pinholes that are at least half a grid spacing away
  ;; from the array border.
  dim = size(pinholeimage, /dim)
  indx = where(cc(0, *) gt dx/2 and cc(1, *) gt dy/2 $
               and  (dim[0] - cc(0, *)) gt dx/2 $
               and  (dim[1] - cc(1, *)) gt dy/2)
  
  ;; The returned coordinates:
  x = cc[0, indx]
  y = cc[1, indx]

end 
