function red_conv_cmap, cmap, img
  inam = 'red_conv_cmap : '
  chan = 0
  nx = img.clip[chan,0,1] - img.clip[chan,0,0] + 1
  ny = img.clip[chan,1,1] - img.clip[chan,1,0] + 1
                                ;
  print, inam + 'clipping cavity map to '+red_stri(nx)+' x '+red_stri(ny)
                                ;
  dim = size(cmap, /dimension)
  mnx = dim[0] - 1
  mny = dim[1] - 1
                                ;print, img.clip[chan,0,0],img.clip[chan,0,1],img.clip[chan,1,0],img.clip[chan,1,1], format='(I5, I5, I5, I5)'
                                ;
                                ; Copy structure to accomodate the chopped cmap
                                ;
  ccmap = {patch:img.patch, clip:img.clip}
  dim = size(img.patch, /dimension)
  for y = 0, dim[1]-1 do begin
     for x = 0, dim[0]-1 do begin
        pnx = img.patch[x,y].xh - img.patch[x,y].xl + 1
        pny = img.patch[x,y].yh - img.patch[x,y].yl + 1
        pmm = dblarr(pnx,pny)
        clip = red_getclips(img, x, y)
                                ; print, clip
                                ; print,  img.patch[x,y].xl, img.patch[x,y].xh, img.patch[x,y].yl, img.patch[x,y].yh
        tpmm = cmap[(clip[0]>0):(clip[1]<mnx),(clip[2]>0):(clip[3]<mny)]
        pmm[*] = median(tpmm)
        pmm[(-clip[0])>0,(-clip[2])>0] = red_fillnan(temporary(tpmm))
        ccmap.patch[x,y].img = red_convolve(temporary(pmm), $
                                            reform(img.patch[x,y].psf[*,*,0]))
        
     endfor
  endfor
  return, temporary(ccmap)
end
