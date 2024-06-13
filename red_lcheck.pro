FUNCTION red_lcheck, l_pos, l_ind, gridstep, LIMIT=lim, RATIO=rat

IF ~keyword_set(rat) THEN rat = [7, 3]
IF ~keyword_set(lim) THEN lim = 0.04

  ;;; find the corner
dr = fltarr(3)
dr[0] = sqrt(total((l_pos[*,0]-l_pos[*,1])^2))
dr[1] = sqrt(total((l_pos[*,0]-l_pos[*,2])^2))
dr[2] = sqrt(total((l_pos[*,1]-l_pos[*,2])^2))
corner = ([2, 1, 0])[where(dr EQ max(dr))]
  ;;; long and short edges
CASE corner OF
    0: BEGIN
        IF dr[0] LT dr[1] THEN BEGIN
            edge_s = 1
            edge_l = 2
        ENDIF ELSE BEGIN
            edge_s = 2
            edge_l = 1
        ENDELSE
    END
    1: BEGIN
        IF dr[0] LT dr[2] THEN BEGIN
            edge_s = 0
            edge_l = 2
        ENDIF ELSE BEGIN
            edge_s = 2
            edge_l = 0
        ENDELSE
    END
    2: BEGIN
        IF dr[1] LT dr[2] THEN BEGIN
            edge_s = 0
            edge_l = 1
        ENDIF ELSE BEGIN
            edge_s = 1
            edge_l = 0
        ENDELSE
    END
ENDCASE
  ;;; re-sort?
ix = [corner, edge_l, edge_s]
IF total(ix NE [0, 1, 2]) NE 0 THEN BEGIN
    l_ind = l_ind[ix]
    l_pos = l_pos[*, ix]
    dr[0] = sqrt(total((l_pos[*,0]-l_pos[*,1])^2))
    dr[1] = sqrt(total((l_pos[*, 0]-l_pos[*, 2])^2))
    dr[2] = sqrt(total((l_pos[*, 1]-l_pos[*, 2])^2))
ENDIF
  ;;; check that it is a right angle triangle with correct side ratio
sl = dr/[rat,sqrt(total(rat^2))]
gridstep = mean(sl)
IF total(abs(sl-gridstep)) GT lim*gridstep THEN BEGIN
      ;;; not our L triangle.
    return, 0
ENDIF

return, 1

END
