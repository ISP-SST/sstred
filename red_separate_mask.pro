FUNCTION red_separate_mask, m, AREA_LIMIT=alim, NO_EDGE_REMOVE=ner, UNLOAD=unl
;+
; NAME:
;       Separate_Mask
; PURPOSE:
;       Find connected areas in a binary mask
; CALLING SEQUENCE:
;       Newmask = Separate_Mask(Mask [, AREA_LIMIT=alim] [, /NO_EDGE_REMOVE])
; INPUTS:
;       Mask:    (byte) 2D threshold mask
;       
; KEYWORDS:
;       AREA_LIMIT:      (int) only consider areas lager than this value
;
;       NO_EDGE_REMOVE:  (bool) do not remove areas connected to the
;                        outer boundary of the mask
;
; OUTPUTS:
;       Newmask: (int) integer array of same size as Mask, entries of
;                connected areas are set to the same value (first
;                area=1, second=2 etc.)
;
; PROCEDURE:
;       
; MODIFICATION HISTORY:
;       17-Jul-2002  P.Suetterlin, SIU
;       15-Jan-2004  Add default for alim
;       15-Oct-2020  Convert to C external
;-

IF NOT keyword_set(alim) THEN alim = 1l
IF NOT keyword_set(ner) THEN ner = 0l

libfile = red_libfile('creduc.so')
IF libfile EQ '' THEN BEGIN
    message, 'Could not find library file sep_mask.so in your DLM_PATH', /info
    return, -1
ENDIF

  ;;; increase mask size by 1 for easier rim handling
s = size(m, /dim)
sx = s(0)+2
sy = s(1)+2
  ;;; this are the 8 neighbour points
sr = [sx, -1, 1, -sx, sx-1, sx+1, -sx-1, -sx+1]

newmask = lonarr(sx, sy)
mask = bytarr(sx, sy)
mask(1, 1) = m GT 0

dummy = call_external(libfile, 'sep_mask', mask, newmask, sx, sy, $
                      long(alim), long(ner), UNLOAD=unl)

return, newmask(1:sx-2, 1:sy-2)

END
