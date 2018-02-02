; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
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
;    no_descatter : in, optional, type=boolean 
;   
;      Don't do back-scatter compensation.
;   
;   
;   
;    nthreads  : 
;   
;   
;   
; 
; 
; :History:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2013-08-27 : MGL. Added support for logging. Let the subprogram
;                find out its own name.
; 
;   2016-02-15 : MGL. Use loadbackscatter. Remove keyword descatter,
;                new keyword no_descatter.
; 
; 
;-
pro red::polcalcube, cam = cam, pref = pref, no_descatter = no_descatter, nthreads = nthreads

  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  cams = [self.camt, self.camr]

  ;; Check polcal_sums
  if(~file_test(self.out_dir + '/polcal_sums', /directory)) then begin
    print, inam + ' : ERROR, folder not found : '+ self.out_dir + '/polcal_sums'
    stop
  endif 
  
  ;; Loop camera
  first = 1B
  for cc = 0, 1 do begin
    if(keyword_set(cam)) then begin
      if(cams[cc] ne cam) then begin
        print, inam + ' : skipping cam -> '+cams[cc]+' != '+cam
        continue
      endif
    endif
    print, inam + ' : processing '+cams[cc]

    ;; Files and states
    f = file_search(self.out_dir + 'polcal_sums/'+cams[cc]+'/camX*', count = count)

    icam = (strsplit(file_basename(f[0]),'.',/extract))[0]
    stat = red_getstates_polcal_out(f)
    upref = stat.pref[uniq(stat.pref, sort(stat.pref))]
    uqw = stat.qws[uniq(stat.qw, sort(stat.qw))]
    ulp = stat.lps[uniq(stat.lp, sort(stat.lp))]
    ulc = stat.lcs[uniq(stat.lcs, sort(stat.lcs))]
    Npref = n_elements(upref)
    nqw = n_elements(uqw)
    nlp = n_elements(ulp)
    nlc = n_elements(ulc)
    
    ;; Load dark
    df = self.out_dir + '/darks/'+icam+'.dark'
    if(~file_test(df)) then begin
      print, inam + ' : ERROR, dark file not found -> '+df
      return
    endif
    dd = f0(df)

    ;; Loop prefilters
    if(first) then begin
      first = 0
      if(keyword_set(pref)) then begin
        idx= where(upref eq pref, count)
        if(count eq 0) then begin
          print, inam + ' : ERROR, user provided prefilter is not on the list -> '+pref
          print, inam + ' : Available prefilters are:'
          for ipref = 0, Npref-1 do print, ipref, +' -> '+upref[ipref], FORMAT='(I3,A)'
          read, ipref, prompt = 'Select prefilter number: '
          upref = upref[ipref]
        endif else upref = upref[idx]
      endif 
    endif

    Npref = n_elements(upref)
    for ipref = 0, Npref-1 do begin
      print, inam + ' : Processing prefilter -> '+upref[ipref]
      if (~keyword_set(no_descatter) AND (upref[ipref] eq '8542' OR upref[ipref] eq '7772')) then begin
        self -> loadbackscatter, icam, upref[ipref], bg, psf
      endif

      ;; Read data
      dim = size(f0(f[0]), /dimension)
      nx = dim[0]
      ny = dim[1]
      d = fltarr(nlc, nqw, nlp, nx, ny)
      d1d = fltarr(nlc, nqw, nlp)

      for pp = 0, nlp - 1 do for qq = 0, nqw - 1 do for ll = 0, nlc-1 do begin
        istate = ulp[pp]+'.'+uqw[qq]+'.'+upref[ipref]+'.'+ulc[ll]
        idx = where(stat.state eq istate, count)
        if count ne 1 then begin
          print, inam + ' : ERROR, irregular state -> '+ istate
          stop
        endif else print, inam + ' : loading -> '+cams[cc]+'.'+istate
        d[ll,qq,pp,*,*] = f0(f[idx]) - dd
        if (~keyword_set(no_descatter) and (upref[ipref] eq '8542' or upref[ipref] eq '7772')) then $
           d[ll,qq,pp,*,*] = red_cdescatter(reform(d[ll,qq,pp,*,*]), bg, psf, /verbose, nthreads = nthreads)
        d1d[ll,qq,pp] = mean(red_fillnan(d[ll,qq,pp,100:nx-101,100:ny-101]))
      endfor                    ; pp

      ;; Save data
      outdir = self.out_dir + '/polcal_cubes/'
      file_mkdir, outdir
      print, inam + ' : saving '+outdir+icam+'.'+upref[ipref]+'.3d.fits'
      ;; fzwrite, temporary(d), outdir+icam+'.'+upref[ipref]+'.3d.fits',' '
      writefits,  outdir+icam+'.'+upref[ipref]+'.3d.fits', temporary(d)

      ;; 1D data and states
      print, inam + ' : saving '+outdir+icam+'.'+upref[ipref]+'.1d.fits'
                                ;fzwrite, d1d, outdir+icam+'.'+upref[ipref]+'.1d.fits',' '
      writefits,outdir+icam+'.'+upref[ipref]+'.1d.fits', d1d
      qw = float(strmid(uqw,2))
      lp = float(strmid(ulp,2))
      print, inam + ' : saving '+outdir+icam+'.'+upref[ipref]+'.qw.fits'
                                ;fzwrite, qw, outdir+icam+'.'+upref[ipref]+'.qw.fits', ' '
      writefits, outdir+icam+'.'+upref[ipref]+'.qw.fits', qw
      print, inam + ' : saving '+outdir+icam+'.'+upref[ipref]+'.lp.fits'
                                ;     fzwrite, lp, outdir+icam+'.'+upref[ipref]+'.lp.fits', ' '
      writefits,  outdir+icam+'.'+upref[ipref]+'.lp.fits', lp

    endfor                      ; ipref
  endfor                        ; cc
  
end
