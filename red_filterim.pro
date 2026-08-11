;; Help subprogram for red::make_pol_crispex
function red_filterim, im, filter
  me = median(im)
  dim = size(im, /dim)
  ft = fft((im - me) * red_taper2(dim[0], dim[1], 1./20.), /double)

  ft *= dcomplex(filter, filter)
  
  return, real_part(fft(ft, /inverse)) + me
end
