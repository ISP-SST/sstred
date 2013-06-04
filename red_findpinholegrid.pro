pro red_findpinholegrid, pinholeimage, simx, simy
;; Find SIMX and SIMY grid coordinates for pinhole array image. Assume
;; the grid is fairly well aligned to the image X and Y directions.

;; From Pit's setup_ph.pro

  ;;; each pinhole gets a unique number
  mask = red_separate_mask(pinholeimage gt .05*max(pinholeimage))
  ;;; # pinholes found
  nph = max(mask)
  ;;; compute PH positions
  cc = fltarr(2, nph)
  FOR i=0, nph-1 DO cc(*, i) = red_com(mask EQ i+1)
  cx = reform(cc(0, *))
  cy = reform(cc(1, *))
  ;;; sort values - PHs should be aligned hor/vert, so this will give a clear
  ;;;               step shape
  cx = cx(sort(cx))
  cy = cy(sort(cy))
  ;;; locate the steps and average the values of each step
  dcx = cx(1:*)-cx
  scx = [-1, where(dcx GT avg(dcx), nx), nph-1]
  simx = intarr(nx+1)
  FOR i=1, nx+1 DO simx(i-1) = round(avg(cx(scx(i-1)+1:scx(i))))
                                ;simx = simx(where((simx GT 32) AND (simx LT sx-32)))
  ;;; now the same for y
  dcy = cy(1:*)-cy
  scy = [-1, where(dcy GT avg(dcy), ny), nph-1]
  simy = intarr(ny+1)
  FOR i=1, ny+1 DO simy(i-1) = round(avg(cy(scy(i-1)+1:scy(i))))
                                ;simy = simy(where((simy GT 32) AND (simy LT sy-32)))

end 
