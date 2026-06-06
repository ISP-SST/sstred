FUNCTION red_markpinholes, ref_m, ref_pos, ratio, title=tit, keep=keep

siz = size(ref_m, /dim)
np = (size(ref_pos,/dim))[1]
device, get_screen_size=s_scr
s_scr -= [20, 50]
FOR scale=1., 8. DO $
  IF min(scale*s_scr/siz) GT 1 THEN BREAK
window, /free, xsi=siz[0]/scale, ysi=siz[1]/scale, title=tit
s_img = rebin(ref_m, siz[0]/scale, siz[1]/scale)
tvscl, s_img
s_win = !d.window
window, /free, xsi=400, ysi=400, title='Zoom'
z_win = !d.window
ref_lind = intarr(3)
FOR i=0, 2 DO BEGIN
    REPEAT BEGIN
        wset, s_win
        cursor, a, b, /dev, /change
        IF !mouse.button EQ 1 THEN BREAK
        
        z_img = ref_m[(a*scale-100) > 0 : (a*scale+99) < (siz[0]-1), $
                      (b*scale-100) > 0 : (b*scale+99) < (siz[1]-1)]
        z_img = rebin(z_img, 2*size(z_img, /dim))
        wset, z_win
        tvscl, z_img
    ENDREP UNTIL 0
    tmp = total(abs(ref_pos-scale*[a, b]#replicate(1, np)), 1)
    ref_lind[i] = where(tmp EQ min(tmp))
    plots, ref_pos[*, ref_lind[i]]/scale, /dev, psy=6, symsi=2
    wait, 1
ENDFOR
wdelete, z_win
IF NOT keyword_set(keep) THEN wdelete, s_win
return, ref_lind
END
