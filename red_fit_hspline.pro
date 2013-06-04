function red_fit_hspline, pp, iwav = iwav, fl = fl, wl = wl
  res = float(intepf(iwav, pp, wl, /linear)) ; linear extrapolation outside bounds
  return, (res - fl)
end
