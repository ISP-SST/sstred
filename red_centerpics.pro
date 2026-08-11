function red_centerpics, ims, _ref_extra = extra
  
  dims = size(ims, /dim)
  Nframes = dims[2]

  cim =  red_centerpic(ims[*, *, 0], _extra = extra)
  cdim = size(cim, /dim)
  
  cims = fltarr([cdim, Nframes])
  cims[*, *, 0] = cim
  
  for iframe = 1, Nframes-1 do cims[*, *, iframe] = red_centerpic(ims[*, *, iframe], _extra = extra)

  return, cims

end
