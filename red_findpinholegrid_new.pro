; docformat = 'rst'

;+
; Find X and Y coordinates for pinholes in an image. 
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
;      The image with pinholes.
;   
;    x : out, type=fltarr
;   
;      The X coordinates of the found pinholes.
;   
;    y : out, type=fltarr
;   
;      The Y coordinates of the found pinholes.
;   
; 
; :Keywords:
; 
;    thres : in, optional, type=float, default=0.1
; 
;       Threshold for identifying a strong enough pinhole.
; 
;    dx : out, optional, type=float
;   
;       The average grid spacing in X. 
; 
;    dy : out, optional, type=float
;   
;       The average grid spacing in Y. 
; 
;    Npinh : out, optional, type=integer
; 
;       The number of pinholes found.
;
;    edgemargin : in, optional, type=integer
; 
;       Margin for momfbd to use for swapping tilts for image shifts.
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
;
;   2013-09-10 : MGL. Return coordinates for each pinhole rather than
;                for a regular grid. Remove pinholes that are closer
;                to an edge than half the grid spacing.
;
;   2014-01-23 : MGL. New keywords thres, dx, dy, Npinh.
;                Documentation. 
;
;   2014-01-27 : MGL. New keyword: margin.
;
;   2015-05-07 : THI. renamed keyword margin to edgemargin to avoid name-clash in the piholecalib routine.
; 
; 
;-
pro red_findpinholegrid_new, pinholeimage, x, y $
                             , dx = dx $
                             , dy = dy $
                             , Npinh = Npinh $
                             , thres = thres $
                             , edgemargin = edgemargin

  if n_elements(thres) eq 0 then thres = 0.1
  if n_elements(edgemargin) eq 0 then edgemargin = 0

  ;; Each pinhole gets a unique ROI number
  mask = red_separate_mask(pinholeimage gt thres*max(pinholeimage))
  ;; # pinholes found
  nph = max(mask)

  ;; Compute PH positions
  cc = fltarr(2, nph)
  FOR i=0, nph-1 DO cc(*, i) = red_com(mask EQ i+1)


  ;; Calculation of grid spacing -----
  ;; Should be rewritten to take rotational misorientation into account.

  cx = reform(cc(0, *))
  cy = reform(cc(1, *))

  ;; Sort values. PHs should be aligned hor/vert, so this will give a
  ;; clear step shape
  cx = cx(sort(cx))
  cy = cy(sort(cy))

  ;; Locate the steps and average the values of each step
  dcx = cx(1:*)-cx
  scx = [-1, where(dcx GT mean(dcx), nx), nph-1]
  simx = intarr(nx+1)
  FOR i=1, nx+1 DO simx(i-1) = round(mean(cx(scx(i-1)+1:scx(i))))

  ;; Now the same for y
  dcy = cy(1:*)-cy
  scy = [-1, where(dcy GT mean(dcy), ny), nph-1]
  simy = intarr(ny+1)
  FOR i=1, ny+1 DO simy(i-1) = round(mean(cy(scy(i-1)+1:scy(i))))

  ;; Grid spacing:
  dx = median(deriv(simx))
  dy = median(deriv(simy))
  
  dd = (dx+dy)/2.

  ;; Selection of pinholes to return ------

  ;; Pick only pinholes that are at least half a grid spacing (plus
  ;; edgemargin) away from the array border.
  dim = size(pinholeimage, /dim)
  indx = where(cc(0, *) gt (dd/2+edgemargin) and cc(1, *) gt (dd/2+edgemargin) $
               and  (dim[0] - cc(0, *)) gt (dd/2+edgemargin) $
               and  (dim[1] - cc(1, *)) gt (dd/2+edgemargin) )
  
  ;; The returned coordinates:
  x = reform(cc[0, indx])
  y = reform(cc[1, indx])

  Npinh = n_elements(x)

end 
