function red_mask_ccd_tabs, img
  
  dim = size(img, /dim)
  x = dindgen(5)-2.d0

  idx0 = [-2,1,1,2]

  for xx = 1,7 do begin     
     idx1 = idx0 + 128*xx
     for yy = 0, dim[1]-1 do begin
        img[128*xx,yy] = cbezier3(x,img[idx1,yy],[0.0d0])
     endfor
  endfor

  return, temporary(img)
end
