function red_encode_scan, hscan, scan
  return, string(long(hscan) * 10000L + long(scan), FORMAT='(I05)')
end
