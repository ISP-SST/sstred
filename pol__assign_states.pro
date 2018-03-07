; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :author:
; 
; 
; 
; 
; :returns:
; 
; 
; :Params:
; 
;    state : 
;   
;   
;   
;    tfiles : 
;   
;   
;   
;    rfiles : 
;   
;   
;   
;    pref : 
;   
;   
;   
;    fdir : 
;   
;   
;   
; 
; :Keywords:
; 
;   camt  : 
;   
;   
;   
;    camr  : 
;   
;   
;   
;    camwb  : 
;   
;   
;   
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2018-03-07 : THI:  Re-writing the CRISP polarimetry to properly
;                apply clips/offsets to the demodulation matrices.
;                Bugfix: exclude scan-number when constructing gain-name.
; 
;-
pro pol::assign_states, state, tfiles, rfiles, pref, fdir,camt = camt, camr = camr, camwb = camwb, newflats = newflats, ftype = ftype
  self.tfiles[*] = tfiles
  self.rfiles[*] = rfiles
  self.state = state
  self.pref = pref
  self.camt = camt
  self.camr = camr
  self.camwb = camwb
  self.scan = (strsplit(state,'.',/extract))[0]
  self.ftype = ftype
  self.tclip = [0,0,0,0]
  self.rclip = [0,0,0,0]
  
  inam = 'pol::assign_states : '

  ; load align_clips
  clipfile = fdir+'/calib/align_clips.'+pref+'.sav'
  if file_test(clipfile) then begin
    restore, clipfile
    wbclip = cl[*,0]
    self.tclip = cl[*,1]
    self.rclip = cl[*,2]
  endif
  
  ;; Wb images?

  for ii = 0L, n_elements(self.tfiles) - 1 do begin
     dir = file_dirname(self.tfiles[ii])
     tmp = dir+'/'+camwb+'.'+strjoin((strsplit(file_basename(self.tfiles[ii]), '.',/extract))[1:*], '.')
     
     if(file_test(tmp)) then begin 
        self.wbfiles[ii] = tmp
        self.destretch = 2B
     endif
  endfor

  ;; Total WB image
   
  ;;state = '0'+strmid(state,1,80)
  tmp = dir+'/'+camwb+'.'+strjoin((strsplit(state,'.',/extract))[0:1], '.')+'.'+ftype
  if(file_test(tmp)) then begin
     self.wb = tmp
     
     if(self.destretch eq 2) then self.destretch = 1B 
  endif else self.destretch = 0B

  ;; Flats

  self.utflat = fdir + '/flats/' + camt + '.' + strjoin((strsplit(file_basename(self.tfiles[0]), '.',/extract))[2:3], '.') + '.unpol.flat'
  self.urflat = fdir + '/flats/' + camr + '.' + strjoin((strsplit(file_basename(self.tfiles[0]), '.',/extract))[2:3], '.') + '.unpol.flat'
;  print, utflat, file_test(utflat), format='(A,I2)'
;  print, urflat, file_test(urflat), format='(A,I2)'

  dum = 1
  if( keyword_set(newflats) ) then dum=2
  for ii = 0L, n_elements(self.tfiles) - 1 do begin
     self.ftfiles[ii] = fdir + '/gaintables/' + camt + '.' + strjoin((strsplit(file_basename(self.tfiles[ii]), '.',/extract))[dum:4], '.')+'.gain'
     self.frfiles[ii] = fdir + '/gaintables/' + camr + '.' + strjoin((strsplit(file_basename(self.tfiles[ii]), '.',/extract))[dum:4], '.')+'.gain'
        
  endfor

  this_tuning = long((strsplit(state,'._',/extract))[3])
  xotfiles = file_search( fdir+'/calib/'+camt+'.*'+pref+'*.xoffs', count=nxot )
  if nxot gt 0 then begin
    tunings = lonarr(nxot)
    for i=0,nxot-1 do begin
      tunings[i] = long((strsplit(file_basename(xotfiles[i]), '._', /extract))[3])
    end
    tunings = abs(tunings-this_tuning)
    idx = where(tunings eq min(tunings))
    if min(idx) ge 0 then self.xotfile = xotfiles[idx[0]]
  endif
  
  if file_test(self.xotfile) then begin
    self.yotfile = red_strreplace( self.xotfile, '.xoffs', '.yoffs' )
    self.xorfile = red_strreplace( self.xotfile, camt, camr )
    self.yorfile = red_strreplace( self.xorfile, '.xoffs', '.yoffs' )
  endif else self.xotfile = ' '

  
end
