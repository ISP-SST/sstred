function chromis_prefilter, par, iwav, pref
    iwav1 = iwav - pref
    res = par[0] / (1.d0+((2.d0*(iwav1 - par[2]) / par[3])^2)^par[4]) * (1.0d0 + par[5]*iwav1 + par[6]*iwav1^3)
  return, res
end
