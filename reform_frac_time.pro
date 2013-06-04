function reform_frac_time, hh, mi, ss
  return, hh * 3600.d + mi * 60.d0 + ss
end
