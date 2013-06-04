function decode_scan, scan, hscan = hscan
  hscan = strmid(scan,0,1)
  return, '0'+strmid(scan,1)
end
