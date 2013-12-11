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
;    tfiles : 
;   
;   
;   
;    rfiles : 
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
; 
;-
function red_getstates_polarim, tfiles, rfiles, fdir,camt = camt, camr = camr, camwb = camwb, newflats=newflats
                                ;
  inam = 'red_getstates_polarim : '
  nt = n_elements(tfiles)
  nr = n_elements(rfiles)
                                ;
                                ; Create arrays
                                ;
  statt = strarr(nt)
  statr = strarr(nr)
  stattnlc = strarr(nt)
                                ;
  pref = strarr(nt)
  scant = strarr(nt)
  scanr = strarr(nr)
  wavt = strarr(nt)
  wavr = strarr(nr)
  lct = strarr(nt)
  lcr = strarr(nr)
  p = '.'
  

                                ;
                                ; T-cam
  for ii = 0L, nt -1 do begin
     tmp = strsplit(file_basename(tfiles[ii]), '.', /extract)
     statt[ii] = tmp[1] + p + tmp[2] + p +tmp[3] + p +tmp[4]
     scant[ii] = tmp[1]
     wavt[ii] = tmp[3]
     lct[ii] = tmp[4]
     pref[ii] = tmp[2]
     stattnlc[ii] = tmp[1] + p + tmp[2] + p +tmp[3]
  endfor 
                                ;
                                ; R-cam
  for ii = 0L, nr -1 do begin
     tmp = strsplit(file_basename(rfiles[ii]), '.', /extract)
     statr[ii] = tmp[1] + p + tmp[2] + p +tmp[3] + p +tmp[4]
     scanr[ii] = tmp[1]
     wavr[ii] = tmp[3]
     lcr[ii] = tmp[4]
  endfor
                                ;
                                ; Get unique states and compare
                                ;
  ustatt = statt[uniq(statt, sort(statt))]
  ustatr = statr[uniq(statr, sort(statr))]
  ulct = lct[uniq(lct, sort(lct))]
                                ;
  nn = n_elements(ustatt)
  star = bytarr(nn)
  ord = lonarr(nn) - 1L
                                ;
                                ; Use T-cam as reference (frames have to exist on both)
                                ;
  for ii = 0L, nn - 1 do begin
     idx = where(ustatr eq ustatt[ii], count)
     if count ne 1 then star[ii] = 1B else ord[ii] = idx
  endfor
  posr = where(ord ne -1, npr)
  post = where(star eq 0B, npt)
                                ;
                                ; keep only those states
                                ;
  ustatt = ustatt[post]
  ustatr = ustatr[posr]
                                ;
                                ; Create array of structures with all the states to demodulate
                                ; 
  nlc = n_elements(ulct)
  nstat = npt / nlc
  ustattnlc = strarr(n_elements(ustatt))
  for jj = 0L, n_elements(ustatt) - 1 do begin
     tmp = strsplit(ustatt[jj], '.',/extract)
     n = n_elements(tmp) - 2
     ustattnlc[jj] = strjoin(tmp[0:n], '.')
  endfor
  ustattnlc = ustattnlc(uniq(ustattnlc, sort(ustattnlc)))
                                ;
  states = {pol, tfiles:strarr(nlc), rfiles:strarr(nlc),state:' ', $
            timg:ptrarr(nlc,/allocate_heap), rimg:ptrarr(nlc,/allocate_heap), $
            immt:ptr_new(), immr:ptr_new(), pref:' ', destretch:0B, wb:' ',$
            wbfiles:strarr(nlc), camt:' ', camr:' ', camwb:' ',$
            x0:0L, x1:0L, y0:0L, y1:0L, telog:' ', ftfiles:strarr(nlc), $
            frfiles:strarr(nlc), utflat:' ', urflat:' ', scan:' '}
                                ;
  pol = objarr(nstat)
  for ii = 0L, nstat - 1 do pol[ii] = obj_new('pol')
                                ;pol = replicate(obj_new('pol'), nstat)
                                ;
  for ii = 0L, nstat - 1 do begin
     

     idx = where(stattnlc eq ustattnlc[ii], count)
     pol[ii] -> assign_states, ustattnlc[ii], tfiles[idx], rfiles[idx], $
        (pref[idx])[0], fdir,camt = camt, camr = camr, camwb = camwb, newflats=newflats
  endfor
                                ;
  return, pol
end
