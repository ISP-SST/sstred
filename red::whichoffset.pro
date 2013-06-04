pro red::whichoffset, state, xoff = xoff, yoff = yoff
                                ; method name
  inam = 'red::whichoffset : '
                                ;
                                ; get all states from calib
                                ;
  pref = (strsplit(state,'.',/extract))[0]
  files = file_search(self.out_dir+'/calib/*'+pref+'*.xoffs', count = count)
  states = file_basename(files, '.xoffs')
                                ;
  for ii = 0L, count - 1 do begin
     tmp = strsplit(states[ii], '.',/extract)
     states[ii] = tmp[1]+'.'+tmp[2]+'.'+tmp[3]
  endfor
  states = states[uniq(states, sort(states))]
                                ;
  if count eq 0 then begin
     print, inam + 'ERROR -> no offset files found in '+self.out_dir+'/calib'
     print, inam + '      -> you must calibrate pinholes for at least 1 state!'
     stop
  endif
                                ;
                                ; select case
                                ;
  pos = where(states eq state, ct)
                                ;
                                ; 
  if(ct eq 1) then begin
     xoff = state+'.xoffs'
     yoff = state+'.yoffs'
  ENDIF ELSE BEGIN
     ws = red_extract_wav(state, lc = lc)
     wss = red_extract_wav(states, lc = lcs)
     p = where((wss - ws) EQ min(wss-ws), ct1)
     IF ct1 GT 1 THEN BEGIN
        pos = (where(lcs[p] EQ lc))[0]
        xoff = states[p[pos]]+'.xoffs'
        yoff = states[p[pos]]+'.yoffs'
     ENDIF ELSE BEGIN
        xoff = states[p]+'.xoffs'
        yoff = states[p]+'.yoffs'
     endelse
                                ;
     
                                ;
  endelse
                                ;
  
  return
end
