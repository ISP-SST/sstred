function pol_mozaic, img, clip = clip
                                ;
  return, mozaic(img.patch.img, $
                 img.patch[*,0].xl, $
                 img.patch[*,0].xh, $
                 img.patch[0,*].yl, $
                 img.patch[0,*].yh, clip= clip)
                                ;
end
