FUNCTION red_markpinholes, ref_m, ref_pos, ratio, title=tit, keep=keep

siz = size(ref_m, /dim)
np = (size(ref_pos,/dim))[1]
device, get_screen_size=s_scr
s_scr -= [20, 50]
FOR scale=1., 8. DO $
  IF min(scale*s_scr/siz) GT 1 THEN BREAK
window, /free, xsi=siz[0]/scale, ysi=siz[1]/scale, title=tit
tvscl, rebin(ref_m, siz[0]/scale, siz[1]/scale)
ref_lind = intarr(3)
FOR i=0, 2 DO BEGIN
    cursor, a, b, /dev
    tmp = total(abs(ref_pos-scale*[a, b]#replicate(1, np)), 1)
    ref_lind[i] = where(tmp EQ min(tmp))
    plots, ref_pos[*, ref_lind[i]]/scale, /dev, psy=6, symsi=2
    wait, 1
ENDFOR
IF NOT keyword_set(keep) THEN wdelete
return, ref_lind
END
