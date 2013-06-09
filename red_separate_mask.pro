; docformat = 'rst'

;+
;
; :author:
;
;    Pit Sutterlin
;
;-
FUNCTION red_separate_mask, m, AREA_LIMIT=alim, NO_EDGE_REMOVE=ner

  IF n_elements(alim) eq 0 THEN alim = 1

  ;;; increase mask size by 1 for easier rim handling
  s = size(m)
  sx = s[1]+2
  sy = s[2]+2
  ;;; this are the 8 neighbour points
  sr = [sx, -1, 1, -sx, sx-1, sx+1, -sx-1, -sx+1]

  mask = bytarr(sx, sy)
  mask[1, 1] = m GT 0

  IF NOT keyword_set(ner) THEN BEGIN
     ;;; First delete all areas in contact with the edges
     ix = 0
     tmp = where(mask[1:sx-2, 1] EQ 1, count)+sx+1
     IF count GT 0 THEN ix = [ix, tmp]
     tmp = sx*(sy-2)+where(mask[1:sx-2, sy-2] EQ 1, count)+1
     IF count GT 0 THEN ix = [ix, tmp]
     tmp = sx*(where(mask[1, 2:sy-3] EQ 1, count)+2)+1
     IF count GT 0 THEN ix = [ix, tmp]
     tmp = sx*(where(mask[sx-2, 2:sy-3] EQ 1, count)+3)-2
     IF count GT 0 THEN ix = [ix, tmp]
     
     mask[ix] = 0
     
     WHILE n_elements(ix) GT 1 DO BEGIN
        ix1 = 0
        FOR i=1, n_elements(ix)-1 DO BEGIN
           FOR j=0, 7 DO BEGIN
              k = ix[i]+sr[j]
              IF mask[k] EQ 1 THEN BEGIN
                 mask[k] = 0
                 ix1 = [ix1, k]
              ENDIF
           ENDFOR
        ENDFOR
        ix = ix1
     ENDWHILE
  ENDIF

  newmask = intarr(sx, sy)
  n_found = 0
  
  ix = where(mask EQ 1, count)
  WHILE count GT 0 DO BEGIN
     x0 = ix[0] MOD sx
     y0 = ix[0]/sx
     mask[x0, y0] = 2
     ix = [0, sx*y0+x0]
     ixall = 0
     
     REPEAT BEGIN
        ixall = [ixall, ix[1:*]]
        ix1 = 0
        FOR i=1l, n_elements(ix)-1 DO BEGIN
           FOR j=0, 7 DO BEGIN
              k = ix[i]+sr[j]
              IF mask[k] EQ 1 THEN BEGIN
                 mask[k] = 2
                 ix1 = [ix1, k]
              ENDIF
           ENDFOR
        ENDFOR
        Ix = ix1
     ENDREP UNTIL n_elements(ix) EQ 1
     ixall = ixall[1:*]
     
     IF n_elements(ixall) GT alim THEN BEGIN
        n_found = n_found+1
        newmask[ixall] = n_found
     ENDIF
     
     mask[ixall] = 0
     ix = where(mask EQ 1, count)
     
  ENDWHILE

  return, newmask[1:sx-2, 1:sy-2]

END
