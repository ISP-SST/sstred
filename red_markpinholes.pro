FUNCTION red_markpinholes, ref_m, ref_pos, ratio, title=tit, keep=keep, camera=cam

l_ori = -1
IF keyword_set(cam) THEN BEGIN
    CASE cam OF
        'Crisp2-W': l_ori = 7
        'Crisp2-T': l_ori = 2
        'Crisp2-R': l_ori = 5
        'Crisp2-D': l_ori = 0
        'Chromis-W': l_ori = 0
        'Chromis-D': l_ori = 7
        'Chromis-N': l_ori = 5
        'Chromis-T': l_ori = 2
        'Chromis-R': l_ori = 5
        ELSE: l_ori = -1
    ENDCASE
    l_im = bytarr(10, 6)
    l_im[1, 1] = 1 & l_im[1, 4] = 1 & l_im[8, 1] = 1
    l_im = rebin(l_im, 200, 120, /sample)
    l_im = rotate(l_im, l_ori)
ENDIF

siz = size(ref_m, /dim)
z_siz = 400
z_4 = z_siz/4
np = (size(ref_pos, /dim))[1]
device, get_screen_size=s_scr
s_scr -= [20, 50]
FOR scale=1., 8. DO $
  IF min(scale*s_scr/siz) GT 1 THEN BREAK
window, /free, xsi=siz[0]/scale, ysi=siz[1]/scale, title=tit
s_img = rebin(ref_m, siz[0]/scale, siz[1]/scale)
tvscl, s_img
s_win = !d.window
IF l_ori GE 0 THEN BEGIN
    window, /free, xsi=z_siz+400, ysi=z_siz, title='Zoom'
    tvscl, l_im, z_siz+100, 140
    xyouts, z_siz+200, 20, /dev, alig=.5, charsiz=2, "Expected orientation for "+cam
ENDIF ELSE window, /free, xsi=z_siz, ysi=z_siz, title='Zoom'
z_win = !d.window
ref_lind = intarr(3)
FOR i=0, 2 DO BEGIN
    REPEAT BEGIN
        wset, s_win
        cursor, a, b, /dev, /change
        IF !mouse.button EQ 1 THEN BREAK
        
        z_img = ref_m[(a*scale-z_4) > 0 : (a*scale+z_4-1) < (siz[0]-1), $
                      (b*scale-z_4) > 0 : (b*scale+z_4-1) < (siz[1]-1)]
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
