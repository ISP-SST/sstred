pro pol::assign_states, state, tfiles, rfiles, pref, fdir,camt = camt, camr = camr, camwb = camwb
  self.tfiles[*] = tfiles
  self.rfiles[*] = rfiles
  self.state = state
  self.pref = pref
  self.camt = camt
  self.camr = camr
  self.camwb = camwb
                                ;
                                ; Wb images?
                                ;
  for ii = 0L, n_elements(self.tfiles) - 1 do begin
     dir = file_dirname(self.tfiles[ii])
     tmp = dir+'/'+camwb+'.'+strjoin((strsplit(file_basename(self.tfiles[ii]), '.',/extract))[1:*], '.')
     
     if(file_test(tmp)) then begin 
        self.wbfiles[ii] = tmp
        self.destretch = 2B
     endif
  endfor
                                ;
                                ; Total WB image
                                ;
                                ;state = '0'+strmid(state,1,80)
  tmp = dir+'/'+camwb+'.'+strjoin((strsplit(state,'.',/extract))[0:1], '.')+'.momfbd'
  if(file_test(tmp)) then begin
     self.wb = tmp
     
     if(self.destretch eq 2) then self.destretch = 1B 
  endif else self.destretch = 0B
                                ;

                                ;
                                ; Flats
                                ;
  self.utflat = fdir + '/flats/' + camt + '.' + strjoin((strsplit(file_basename(self.tfiles[0]), '.',/extract))[2:3], '.') + '.unpol.flat'
  self.urflat = fdir + '/flats/' + camr + '.' + strjoin((strsplit(file_basename(self.tfiles[0]), '.',/extract))[2:3], '.') + '.unpol.flat'
;  print, utflat, file_test(utflat), format='(A,I2)'
;  print, urflat, file_test(urflat), format='(A,I2)'

  for ii = 0L, n_elements(self.tfiles) - 1 do begin
     self.ftfiles[ii] = fdir + '/gaintables/' + camt + '.' + strjoin((strsplit(file_basename(self.tfiles[ii]), '.',/extract))[2:4], '.')+'.gain'
     self.frfiles[ii] = fdir + '/gaintables/' + camr + '.' + strjoin((strsplit(file_basename(self.tfiles[ii]), '.',/extract))[2:4], '.')+'.gain'
                                ;
  endfor
                                ;
  return
end
