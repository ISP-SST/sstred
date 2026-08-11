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
; 
; :Keywords:
; 
;    xoff  : 
;   
;   
;   
;    yoff  : 
;   
;   
;   
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2014-01-03 : PS  fix computation of nearest state
;-
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
     wss = abs(wss-ws[0])
     p = where(wss EQ min(wss), ct1)
     IF ct1 GT 1 THEN BEGIN
         IF ct1 EQ 2 THEN BEGIN
               ;;; This is two wavelengths in same distance.  Just pick the first
             xoff = states[p[0]]+'.xoffs'
             yoff = states[p[0]]+'.yoffs'
         ENDIF ELSE BEGIN 
               ;;; one wavelength with 4 lc states
             pos = (where(lcs[p] EQ lc[0]))[0]
             xoff = states[p[pos]]+'.xoffs'
             yoff = states[p[pos]]+'.yoffs'
         ENDELSE
     ENDIF ELSE BEGIN
        xoff = states[p]+'.xoffs'
        yoff = states[p]+'.yoffs'
     endelse
  endelse
  return
end
