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
;    tstates : 
;   
;   
;   
;    rstates : 
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
;   2013-12-12 : PS. Move struct definition to separate file
;                (pol__define.pro)
;
;   2018-10-04 : MGL. Adapt to new code base. Use states rather than
;                file names as input parameters.
; 
;-
function red_getstates_polarim, tstates, rstates, fdir, camt = camt, camr = camr, camwb = camwb, newflats=newflats

  inam = 'red_getstates_polarim : '
  Nt = n_elements(tstates)
  Nr = n_elements(rstates)

  if Nt eq Nr then Nfiles = Nt else stop

  ;; Create arrays
;  statt = strarr(Nfiles)        ; e.g. 00044.6302.6302_-1120.lc1
;  statr = strarr(Nfiles)        ; e.g. 00044.6302.6302_-1120.lc1
;  stattnlc = strarr(Nfiles)     ; e.g. 00043.6302.6302_-1120

;  pref = strarr(Nfiles)
;  scant = strarr(Nfiles)
;  scanr = strarr(Nfiles)
;  wavt = strarr(Nfiles)
;  wavr = strarr(Nfiles)
;  lct = strarr(Nfiles)
;  lcr = strarr(Nfiles)
;  p = '.'
  
  pref = tstates.prefilter
  statnlc = string(tstates.scannumber,format='(i05)') $
            + '.' + tstates.prefilter $
            + '.' + tstates.tuning
  lct = 'lc' + string(tstates.lc,format='(i1)')
  lcr = 'lc' + string(rstates.lc,format='(i1)')
  statt = statnlc + '.' + lct
  statr = statnlc + '.' + lcr

  
;     ;; T-cam
;     for ifile = 0L, Nfiles -1 do begin
;       tmp = strsplit(file_basename(tstates[ifile].filename), '._', /extract)
;   ;    statt[ifile] = tmp[1] + p + tmp[2] + p +tmp[3] + p +tmp[4]
;   ;    scant[ifile] = tmp[1]
;   ;    wavt[ifile] = tmp[3]
;       lct[ifile] = tmp[4]
;       pref[ifile] = tmp[2]
;       statnlc[ifile] = tmp[1] + p + tmp[2] + p +tmp[3]
;     endfor                        ; ifile
;   
;     ;; R-cam
;     for ifile = 0L, Nfiles -1 do begin
;       tmp = strsplit(file_basename(rstates[ifile].filename), '.', /extract)
;   ;    statr[ifile] = tmp[1] + p + tmp[2] + p +tmp[3] + p +tmp[4]
;   ;    scanr[ifile] = tmp[1]
;   ;    wavr[ifile] = tmp[3]
;       lcr[ifile] = tmp[4]
;     endfor                        ; ifile

  ;; Get unique states and compare
  ustatt = statt[uniq(statt, sort(statt))]
  ustatr = statr[uniq(statr, sort(statr))]
  ulct = lct[uniq(lct, sort(lct))]
                                ;
  Nn = n_elements(ustatt)
  star = bytarr(nn)
  ord = replicate(-1L, Nn) 

  ;; Use T-cam as reference (frames have to exist on both)
  for ii = 0L, Nn - 1 do begin
    idx = where(ustatr eq ustatt[ii], count)
    if count ne 1 then star[ii] = 1B else ord[ii] = idx
  endfor
  posr = where(ord ne -1, Npr)
  post = where(star eq 0B, Npt)

  ;; keep only those states
  ustatt = ustatt[post]
  ustatr = ustatr[posr]

  ;; Create array of structures with all the states to demodulate
  Nlc = n_elements(ulct)
  Nstat = Npt / Nlc
  ustatnlc = strarr(n_elements(ustatt))
  for jj = 0L, n_elements(ustatt) - 1 do begin
    tmp = strsplit(ustatt[jj], '.',/extract)
    n = n_elements(tmp) - 2
    ustatnlc[jj] = strjoin(tmp[0:n], '.')
  endfor
  ustatnlc = ustatnlc(uniq(ustatnlc, sort(ustatnlc)))

  pol = objarr(Nstat)
  for istat = 0L, Nstat - 1 do pol[istat] = obj_new('pol')
  ;;pol = replicate(obj_new('pol'), nstat)

  for istat = 0L, Nstat - 1 do begin
    idx = where(statnlc eq ustatnlc[istat], count)
    pol[istat] -> assign_states $
       , ustatnlc[istat] $
       , tstates[idx].filename $
       , rstates[idx].filename $
       , (pref[idx])[0] $
       , fdir $
       , camt = camt $
       , camr = camr $
       , camwb = camwb $
       , newflats = newflats
  endfor

  return, pol

end
