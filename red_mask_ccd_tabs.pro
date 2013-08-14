function red_mask_ccd_tabs, img
  
  dim = size(img, /dim)
  idx0 = [-1,1]
  
  for xx = 1,7 do begin     
     idx1 = idx0 + 128*xx
     img[128*xx,*] = total(img[idx1,*],1) * 0.5
  endfor

  return, temporary(img)
end
