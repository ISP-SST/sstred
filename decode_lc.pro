function decode_lc, lcin
  return, 'lc'+ red_stri(long(strmid(lcin, 2)) mod 4)
end
