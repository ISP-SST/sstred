function encode_lc, lcin, scan
  return, 'lc' + red_stri(long(strmid(lcin, 2)) + 4L * scan)
end
