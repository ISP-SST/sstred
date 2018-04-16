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
;   2018-02-02 : MGL. Adapt to new codebase.
; 
;   2018-04-16 : MGL. Write single FITS files with extensions.
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

  ;; Files and states
  files = file_search(self.out_dir + 'polcal_sums/*/*.fits', count = count)
  self -> extractstates, files, states, /polcal

  ;; Get cameras
  cams = (states[uniq(states.camera, sort(states.camera))]).camera

  ;; Loop camera
  for icam = 0, n_elements(cams)-1 do begin
    
    if(keyword_set(cam)) then begin
      if(cams[icam] ne cam) then begin
        print, inam + ' : skipping cam -> '+cams[icam]+' != '+cam
        continue
      endif
    endif
    print, inam + ' : processing '+cams[icam]
    
    self -> selectfiles, files = files, states = states $
                         , cam = cams[icam], sel = sel
    selstates = states[sel]
    selfiles = files[sel]

    detector = selstates[0].detector

    upref = (states[uniq(selstates.prefilter, sort(states.prefilter))]).prefilter
    uqw = (selstates[uniq(selstates.qw, sort(selstates.qw))]).qw
    ulp = (selstates[uniq(selstates.lp, sort(selstates.lp))]).lp
    ulc = (selstates[uniq(selstates.lc, sort(selstates.lc))]).lc

    Npref = n_elements(upref)
    Nqw = n_elements(uqw)
    Nlp = n_elements(ulp)
    Nlc = n_elements(ulc)
    
    ;; Load dark
    self -> get_calib, selstates[0], darkdata = dd, status = status

    ;; Take pref keyword into account
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

    Npref = n_elements(upref)
    ;; Loop prefilters
    for ipref = 0, Npref-1 do begin
      
      print, inam + ' : Processing prefilter -> '+upref[ipref]

      dodescatter = ~keyword_set(no_descatter) $
                    and (upref[ipref] eq '8542' OR upref[ipref] eq '7772')
      
      if dodescatter then begin
        self -> loadbackscatter, detector, upref[ipref], bg, psf
      endif

      ;; Read data
      dim = fxpar(headfits(files[0]), 'NAXIS*')
      Nx = dim[0]
      Ny = dim[1]
      d = fltarr(Nlc, Nqw, Nlp, Nx, Ny)
      d1d = fltarr(Nlc, Nqw, Nlp)

      iloop = 0
      Nloop = Nlp*Nqw*Nlc
      for ilp = 0, Nlp - 1 do begin
        for iqw = 0, Nqw - 1 do begin
          for ilc = 0, Nlc-1 do begin

            undefine, fullstate_list
            red_append, fullstate_list, 'lp'+string(round(ulp[ilp]), format = '(i03)')
            red_append, fullstate_list, 'qw'+string(round(uqw[iqw]), format = '(i03)')
            red_append, fullstate_list, upref[ipref]
            red_append, fullstate_list, '*_*'
            red_append, fullstate_list, 'lc'+strtrim(long(ulc[ilc]), 2)
            statestring = strjoin(fullstate_list, '_')
            
            indx = where(strmatch(selstates.fullstate,statestring ),count)
            if count ne 1 then begin
              print, inam + ' : ERROR, irregular state -> '+ statestring
              stop
            endif 
            red_progressbar, iloop, Nloop, /predict, cams[icam]+' : '+file_basename(selfiles[indx])
            d[ilc,iqw,ilp,*,*] = red_readdata(selfiles[indx], /silent) - dd
            if dodescatter then $
               d[ilc,iqw,ilp,*,*] = red_cdescatter(reform(d[ilc,iqw,ilp,*,*]) $
                                                   , bg, psf, nthreads = nthreads)
            d1d[ilc,iqw,ilp] = mean(d[ilc,iqw,ilp,100:Nx-101,100:Ny-101], /nan)

            iloop++
            
          endfor                ; ilc
        endfor                  ; iqw
      endfor                    ; ilp

      ;; Save data
      outdir = self.out_dir + '/polcal_cubes/'
      file_mkdir, outdir
      
      pname = outdir+detector+'_'+upref[ipref]+'_polcalcube.fits'
      print, inam + ' : saving '+pname
      writefits, pname, d

      mkhdr, ehdr, d1d, /image
      red_fitsaddkeyword, anchor = anchor, ehdr, 'EXTNAME', 'D1D', '1D polcal data'
      writefits, pname, d1d, ehdr, /append

      mkhdr, ehdr, uqw, /image
      red_fitsaddkeyword, anchor = anchor, ehdr, 'EXTNAME', 'QW', 'Quarter wave plate angles'
      writefits, pname, uqw, ehdr, /append

      mkhdr, ehdr, ulp, /image
      red_fitsaddkeyword, anchor = anchor, ehdr, 'EXTNAME', 'LP', 'Linear polarrizer angles'
      writefits, pname, ulp, ehdr, /append

    endfor                      ; ipref
  endfor                        ; icam
  
end
