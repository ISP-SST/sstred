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
;    f : 
;   
;   
;   
; 
; :Keywords:
; 
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
; 
;-
function red_get_stkstates, f
  nt = n_elements(f)
  pref = strarr(nt)
  wav = strarr(nt)
  iwav = intarr(nt)
  scan = strarr(nt)
  lscan = lonarr(nt)
  ord = lindgen(nt)
  state = strarr(nt)
                                ;
                                ; Extract states for each file
                                ;
  for ii = 0L, nt -1 do begin
     tmp = strsplit(file_basename(f[ii]), '.', /extract)
     scan[ii] = tmp[1]
     lscan[ii] = long(tmp[1])
     pref[ii] = tmp[2]
     wav[ii] = tmp[3]
     iwav[ii] = fix((strsplit(tmp[3],'_',/extract))[1])
     state[ii] = strjoin(tmp[1:3],'.')
  endfor
                                ;
                                ; Uniq scan numbers and wavelengths
                                ;
  uscan = scan[uniq(lscan, sort(lscan))]
  uwav = wav[uniq(iwav, sort(iwav))]
  uiwav = iwav[uniq(iwav, sort(iwav))]
  nscan = n_elements(uscan)
  nwav = n_elements(uwav)
                                ;
                                ; Ordered arrays
                                ;
  idx = lonarr(nwav, nscan) -1L
  flag = bytarr(nwav, nscan)
  ofiles = strarr(nwav, nscan)
  owav = strarr(nwav, nscan)
  ostate = strarr(nwav, nscan)
                                ;
                                ; Assign elements
                                ;
  k = 0L
  for ss = 0L, nscan - 1 do begin
     for ww = 0L, nwav - 1 do begin
        pos = where((scan eq uscan[ss]) AND (wav eq uwav[ww]), count)
        if(count eq 1) then begin
           idx[ww,ss] = ord[pos]
           ofiles[ww,ss] = f[pos]
           owav[ww,ss] = wav[pos]
           ostate[ww,ss] = state[pos]
        endif else flag[ww,ss] = 1B
     endfor
  endfor
                                ;
  return, {files:f, ofiles:ofiles, flag:flag, ostate:ostate, state:state, idx:idx, uwav:uwav, $
           uscan:uscan, uiwav:uiwav, nscan:nscan, nwav:nwav}
end
