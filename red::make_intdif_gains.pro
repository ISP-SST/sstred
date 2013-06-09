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
; 
; :Keywords:
; 
;    pref  : 
;   
;   
;   
;    scan  : 
;   
;   
;   
;    cam  : 
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
pro red::make_intdif_gains, pref = pref, scan = scan, cam = cam
  inam = 'red::make_intdif_gains : '
  outdir = self.out_dir + '/gaintables/'
  file_mkdir, outdir

                                ;
                                ; Get int differences
                                ;
  root = self.out_dir + '/cmap_intdif/cam*intdiff'
  f = file_search(root, count = ct)
  if(ct eq 0) then begin
     print, inam +'ERROR, no files found in ' + self.out_dir + '/cmap_intdif/'
  endif else print, inam + red_stri(ct)+' files found'

                                ;
                                ; get states
                                ;
  mpref = strarr(ct)
  mcam = strarr(ct)
  mscan = strarr(ct)
  mwav = strarr(ct)
  mlc = strarr(ct)
  mstate = strarr(ct)
                                ;
  for ii = 0L, ct  -1 do begin
     tmp = strsplit(file_basename(f[ii]),'.',/extract)
     mcam[ii] = tmp[0]
     mscan[ii] = tmp[1]
     mpref[ii] = tmp[2]
     mwav[ii] = tmp[3]
     mlc[ii] = tmp[4]
     mstate[ii] = strjoin([tmp[0],tmp[1],tmp[2],tmp[3],tmp[4]], '.')
  endfor
  
                                ;
                                ; Unique states
                                ;
  ucam = mcam[uniq(mcam, sort(mcam))]
  upref = mpref[uniq(mpref, sort(mpref))]
  uwav = mwav[uniq(mwav, sort(mwav))]
  uscan = mscan[uniq(mscan, sort(mscan))]
  ulc = mlc[uniq(mlc, sort(mlc))]
                                ;
  ncam = n_elements(ucam)
  nscan = n_elements(uscan)
  nwav = n_elements(uwav)
  nlc = n_elements(ulc)
  npref = n_elements(upref)
  
                                ;
                                ; Loop
                                ;
  for cc = 0L, ncam -1 do begin
                                ;
     if(keyword_set(cam)) then begin
        if(ucam[cc] ne cam) then begin
           print, inam + 'skipping ' + ucam[cc]+' != '+cam
           continue
        endif
     endif
                                ;
     for pp = 0L, npref - 1 do begin
        if(keyword_set(pref)) then begin
           if(upref[pp] ne pref) then begin
              print, inam + 'skipping ' + pref[pp] + ' != '+pref
              continue
           endif
        endif
        
        for ww = 0L, nwav - 1 do begin
           fstate = strjoin([ucam[cc], upref[pp], uwav[ww]],'.')
           
           ff = self.out_dir + 'flats/'+fstate+'.unpol.flat'
           if(~file_test(ff)) then begin
              print, inam + 'file not found -> '+ff
              print, inam + 'skipping!'
              continue
           endif else begin
              print, inam + 'loading flat -> '+ff
              flat = f0(ff)
           endelse

           for ss = 0L, nscan-1 do begin
              if(keyword_set(scan)) then begin
                 if(uscan[ss] ne scan) then begin
                    print, inam + 'skipping ' + scan[ss] + ' != '+scan
                    continue
                 endif
              endif
              
              for ll = 0L, nlc - 1 do begin
                 istate = strjoin([ucam[cc],uscan[ss], upref[pp],uwav[ww],ulc[ll]],'.')
                 idx = where(mstate eq istate, nfiles)
                 if(nfiles eq 0) then continue
                 
                                ;
                                ; Construct gains
                                ;
                 print, inam + 'loading -> '+f[idx]
                 rat = f0(f[idx]) * 1.e-4
                 gain = red_flat2gain(flat*rat) 
                 
                                ;
                                ; Save result
                                ;
                 outname = outdir+'/'+istate+'.gain'
                 print, inam + 'saving -> ' + outname
                 fzwrite, gain, outname, ' '

              endfor
           endfor
        endfor                  ; ww
     endfor                     ; pp

  endfor                        ; cc
  
end
