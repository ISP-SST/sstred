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
;    cam  : 
;   
;   
;   
;    pref  : 
;   
;   
;   
;    descatter  : 
;   
;   
;   
;    nthreads  : 
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
pro red::polcalcube, cam = cam, pref = pref, descatter = descatter, nthreads = nthreads
  inam = 'red::polcalcube : '
  cams = [self.camt, self.camr]

                                ;
                                ; Check polcal_sums
                                ;
  if(~file_test(self.out_dir + '/polcal_sums', /directory)) then begin
     print, inam + 'ERROR, folder not found : '+ self.out_dir + '/polcal_sums'
     stop
  endif 
  
                                ;
                                ; Loop camera
                                ;
  first = 1B
  for cc = 0, 1 do begin
     if(keyword_set(cam)) then begin
        if(cams[cc] ne cam) then begin
           print, inam + 'skipping cam -> '+cams[cc]+' != '+cam
           continue
        endif
     endif
     print, inam + 'processing '+cams[cc]


                                ;
                                ; Files and states
                                ;
     f = file_search(self.out_dir + 'polcal_sums/'+cams[cc]+'/camX*', count = count)

     icam = (strsplit(file_basename(f[0]),'.',/extract))[0]
     stat = red_getstates_polcal_out(f)
     upref = stat.pref[uniq(stat.pref, sort(stat.pref))]
     uqw = stat.qws[uniq(stat.qw, sort(stat.qw))]
     ulp = stat.lps[uniq(stat.lp, sort(stat.lp))]
     ulc = stat.lcs[uniq(stat.lcs, sort(stat.lcs))]
     npref = n_elements(upref)
     nqw = n_elements(uqw)
     nlp = n_elements(ulp)
     nlc = n_elements(ulc)
     


                                ;
                                ; Load dark
                                ;
     df = self.out_dir + '/darks/'+icam+'.dark'
     if(~file_test(df)) then begin
        print, inam + 'ERROR, dark file not found -> '+df
        return
     endif
     dd = f0(df)

                                ;
                                ; Loop prefilters
                                ;
     if(first) then begin
        first = 0
        if(keyword_set(pref)) then begin
           idx= where(upref eq pref, count)
           if(count eq 0) then begin
              print, inam + 'ERROR, user provided prefilter is not on the list -> '+pref
              print, inam + 'Available prefilters are:'
              for ii = 0, npref-1 do print, ii, +' -> '+upref[ii], FORMAT='(I3,A)'
              read, ii, prompt = 'Select prefilter number: '
              ipref = upref[ii]
              
           endif else ipref = pref
           
        endif else begin
           if(npref gt 1) then begin
              print, inam + 'Found prefilters:'
              for ii = 0, npref-1 do print, ii, +' -> '+upref[ii], FORMAT='(I3,A)'
              read, ii, prompt = 'Select prefilter number: '
              ipref = upref[ii]
           endif else begin
              ipref = upref[0]
           endelse 
        endelse
     endif

     print, inam + 'Selected prefilter -> '+ipref
     if(keyword_set(descatter) AND (ipref eq '8542' OR ipref eq '7772')) then begin
        print, inam + 'loading descatter data for '+icam
        bg =  f0(self.descatter_dir + '/' + icam + '.backgain.f0')
        psf = f0(self.descatter_dir + '/' + icam + '.psf.f0')
     endif
                                ;
                                ; read data
                                ;
     dim = size(f0(f[0]), /dimension)
     nx = dim[0]
     ny = dim[1]
     d = fltarr(nlc, nqw, nlp, nx, ny)
     d1d = fltarr(nlc, nqw, nlp)
                                ;
     for pp = 0, nlp - 1 do for qq = 0, nqw - 1 do for ll = 0, nlc-1 do begin
        istate = ulp[pp]+'.'+uqw[qq]+'.'+ipref+'.'+ulc[ll]
        idx = where(stat.state eq istate, count)
        if count ne 1 then begin
           print, inam + 'ERROR, irregular state -> '+ istate
           stop
        endif else print, inam + 'loading -> '+cams[cc]+'.'+istate
        d[ll,qq,pp,*,*] = f0(f[idx]) - dd
        if(keyword_set(descatter) AND (ipref eq '8542' or ipref eq '7772')) then d[ll,qq,pp,*,*] = red_cdescatter(reform(d[ll,qq,pp,*,*]), bg, psf, /verbose, nthreads = nthreads)

        d1d[ll,qq,pp] = mean(red_fillnan(d[ll,qq,pp,100:nx-101,100:ny-101]))
     endfor

                                ;
                                ; Save data
                                ;
     outdir = self.out_dir + '/polcal_cubes/'
     file_mkdir, outdir
     print, inam + 'saving '+outdir+icam+'.'+ipref+'.3d.f0'
     fzwrite, temporary(d), outdir+icam+'.'+ipref+'.3d.f0',' '


                                ;
                                ; 1D data and states
                                ;
     print, inam + 'saving '+outdir+icam+'.'+ipref+'.1d.f0'
     fzwrite, d1d, outdir+icam+'.'+ipref+'.1d.f0',' '
     qw = float(strmid(uqw,2))
     lp = float(strmid(ulp,2))
     print, inam + 'saving '+outdir+icam+'.'+ipref+'.qw.f0'
     fzwrite, qw, outdir+icam+'.'+ipref+'.qw.f0', ' '
     print, inam + 'saving '+outdir+icam+'.'+ipref+'.lp.f0'
     fzwrite, lp, outdir+icam+'.'+ipref+'.lp.f0', ' '
     

  endfor
  
end
