function red_prefilter, pp, wav
  res = fltarr(n_elemets(wav))
                                ;
  p0 = pp[0]
  p1 = pp[1]
  p2 = pp[2]
  p3 = pp[3]
;  p4 = pp[4]
  p5 = pp[5]
  p6 = pp[6]
                                ; p7 = pp[7]
                                ;
  res = p0 / (1.d0+((2.d0*(wav - p1) / p2)^2.)^p3) * (1.d0 + (wav * p5) + (wav^3. * p6))

  return, res
end
