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

  ;; Check polcal_sums
  if(~file_test(self.out_dir + '/polcal_sums', /directory)) then begin
    print, inam + ' : ERROR, folder not found : '+ self.out_dir + '/polcal_sums'
    stop
  endif

  cams = [self.camt, self.camr]
  ;; Loop camera
  first = 1B
  for icam = 0, 1 do begin
    if(keyword_set(cam)) then begin
      if(cams[icam] ne cam) then begin
        print, inam + ' : skipping cam -> '+cams[icam]+' != '+cam
        continue
      endif
    endif
    print, inam + ' : processing '+cams[icam]

    ;; Files and states
    files = file_search(self.out_dir + 'polcal_sums/'+cams[icam]+'/camX*', count = count)

    detector = (strsplit(file_basename(files[0]),'.',/extract))[0]
    stat = red_getstates_polcal_out(files)
    upref = stat.pref[uniq(stat.pref, sort(stat.pref))]
    uqw = stat.qws[uniq(stat.qw, sort(stat.qw))]
    ulp = stat.lps[uniq(stat.lp, sort(stat.lp))]
    ulc = stat.lcs[uniq(stat.lcs, sort(stat.lcs))]
    Npref = n_elements(upref)
    Nqw = n_elements(uqw)
    Nlp = n_elements(ulp)
    Nlc = n_elements(ulc)
    
    ;; Load dark
    df = self.out_dir + '/darks/'+detector+'.dark'
    if(~file_test(df)) then begin
      print, inam + ' : ERROR, dark file not found -> '+df
      return
    endif
    dd = f0(df)

    ;; Loop prefilters
    if(first) then begin
      first = 0
      if keyword_set(pref) then begin
        indx = where(upref eq pref, count)
        if(count eq 0) then begin
          print, inam + ' : ERROR, user provided prefilter is not on the list -> '+pref
          print, inam + ' : Available prefilters are:'
          for ipref = 0, Npref-1 do print, ipref, +' -> '+upref[ipref], FORMAT='(I3,A)'
          read, ipref, prompt = 'Select prefilter number: '
          upref = upref[ipref]
        endif else upref = upref[indx]
      endif 
    endif

    Npref = n_elements(upref)
    for ipref = 0, Npref-1 do begin
      
      print, inam + ' : Processing prefilter -> '+upref[ipref]

      dodescatter = ~keyword_set(no_descatter) $
                    and (upref[ipref] eq '8542' OR upref[ipref] eq '7772')
      
      if dodescatter then begin
        self -> loadbackscatter, detector, upref[ipref], bg, psf
      endif

      ;; Read data
      dim = size(f0(f[0]), /dimension)
      Nx = dim[0]
      Ny = dim[1]
      d = fltarr(Nlc, Nqw, Nlp, Nx, Ny)
      d1d = fltarr(Nlc, Nqw, Nlp)

      for ilp = 0, Nlp - 1 do begin
        for iqw = 0, Nqw - 1 do begin
          for ilc = 0, Nlc-1 do begin
            istate = ulp[ilp]+'.'+uqw[iqw]+'.'+upref[ipref]+'.'+ulc[ilc]
            indx = where(stat.state eq istate, count)
            if count ne 1 then begin
              print, inam + ' : ERROR, irregular state -> '+ istate
              stop
            endif else print, inam + ' : loading -> '+cams[icam]+'.'+istate
            d[ilc,iqw,ilp,*,*] = f0(files[indx]) - dd
            if dodescatter then $
               d[ilc,iqw,ilp,*,*] = red_cdescatter(reform(d[ilc,iqw,ilp,*,*]) $
                                                   , bg, psf, /verbose, nthreads = nthreads)
            d1d[ilc,iqw,ilp] = mean(d[ilc,iqw,ilp,100:Nx-101,100:Ny-101], /nan)
          endfor                ; ilc
        endfor                  ; iqw
      endfor                    ; ilp

      ;; Save data
      outdir = self.out_dir + '/polcal_cubes/'
      file_mkdir, outdir
      print, inam + ' : saving '+outdir+detector+'.'+upref[ipref]+'.3d.fits'
      writefits, outdir+detector+'.'+upref[ipref]+'.3d.fits', temporary(d)

      ;; 1D data and states
      print, inam + ' : saving '+outdir+detector+'.'+upref[ipref]+'.1d.fits'
      writefits,outdir+detector+'.'+upref[ipref]+'.1d.fits', d1d
      qw = float(strmid(uqw,2))
      lp = float(strmid(ulp,2))
      print, inam + ' : saving '+outdir+detector+'.'+upref[ipref]+'.qw.fits'
      writefits, outdir+detector+'.'+upref[ipref]+'.qw.fits', qw
      print, inam + ' : saving '+outdir+detector+'.'+upref[ipref]+'.lp.fits'
      writefits, outdir+detector+'.'+upref[ipref]+'.lp.fits', lp

    endfor                      ; ipref
  endfor                        ; icam
  
end
