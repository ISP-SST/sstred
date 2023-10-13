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
;   2023-10-13 : MGL. Do gain correction. Optionally with nearest
;                tuning if polcal tuning gains not available. 
; 
;-
pro red::polcalcube, cam = cam, pref = pref, no_descatter = no_descatter, nthreads = nthreads

  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)  
  
  ;; Check polcal_sums
  if(~file_test(self.out_dir + '/polcal_sums', /directory)) then begin
    print, inam + ' : ERROR, folder not found : '+ self.out_dir + '/polcal_sums'
    stop
  endif

  files = file_search(self.out_dir + 'polcal_sums/*/*.fits', count = count)
  self -> extractstates, files, states, /polcal

  upref = (states[uniq(states.prefilter, sort(states.prefilter))]).prefilter
  cams = (states[uniq(states.camera, sort(states.camera))]).camera

  ;; Take pref keyword into account
  if keyword_set(pref) then begin
    indx = where(upref eq pref, count)
    if count eq 0 then begin
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

    undefine, use_gain
    
    dodescatter = ~keyword_set(no_descatter)  AND self.dodescatter $
                  and (upref[ipref] eq '8542' OR upref[ipref] eq '7772')
    
;    ;; Files and states
;    files = file_search(self.out_dir + 'polcal_sums/*/*_'+upref[ipref]+'_*.fits', count = count)
;    self -> extractstates, files, states, /polcal
 

    ;; Loop cameras
    for icam = 0, n_elements(cams)-1 do begin
      
      if(keyword_set(cam)) then begin
        if(cams[icam] ne cam) then begin
          print, inam + ' : skipping cam -> '+cams[icam]+' != '+cam
          continue
        endif
      endif
      print, inam + ' : processing '+cams[icam]
      
      self -> selectfiles, files = files, states = states $
                           , cam = cams[icam], pref = upref[ipref], sel = sel
      selstates = states[sel]
      selfiles = files[sel]

      detector = selstates[0].detector

      uqw = (selstates[uniq(selstates.qw, sort(selstates.qw))]).qw
      ulp = (selstates[uniq(selstates.lp, sort(selstates.lp))]).lp
      ulc = (selstates[uniq(selstates.lc, sort(selstates.lc))]).lc

      Npref = n_elements(upref)
      Nqw = n_elements(uqw)
      Nlp = n_elements(ulp)
      Nlc = n_elements(ulc)

      ;; Load dark
      self -> get_calib, selstates[0] $
                         , darkdata = dd, darkstatus  = darkstatus $
                         , gainname = gn 
      if darkstatus ne 0 then stop

      if file_test(gn) then begin
        gg = red_readdata(gn)
      endif else begin
        ;; Do we want to use gains from a nearby wavelength?
        gfiles = file_search('gaintables/' $
                             + selstates[0].detector + '_' $
                             + selstates[0].cam_settings + '_' $
                             + selstates[0].prefilter + '_' $
                             + '*' $
                             + 'lc'+strtrim(long(selstates[0].lc), 2) $
                             + '.gain.fits' $
                             , count = Ng )
                                ; get tuning from each of gfiles
        red_extractstates, gfiles, dwav = gtun_wavelength, wav = gtuning
        mn = min(abs(gtun_wavelength - selstates[0].tun_wavelength * 1e10),minloc)
        if n_elements(use_gain) eq 0 then begin
          print
          print, inam + ' : There are no flat field data for the polcal wavelength.'
          print, '   Tuning for the polcal data: '+selstates[0].tuning
          print, '   Nearest flats tuning: '+gtuning[minloc]
          s = ''
          read, '    Do you want to use it to flat-field the polcal data [Y/n]? ', s
          if strlen(s) eq 0 || ~(strlowcase(strmid(s, 0, 1)) eq 'n') then use_gain = 1 else use_gain = 0
        endif
      endelse 

      if dodescatter then begin
        self -> loadbackscatter, detector, upref[ipref], bg, psf
      endif

      ;; Read data
      dim = fxpar(headfits(files[0]), 'NAXIS*')
      Nx = dim[0]
      Ny = dim[1]
      d = fltarr(Nlc, Nqw, Nlp, Nx, Ny)
      d1d = fltarr(Nlc, Nqw, Nlp)

      gains = fltarr(Nx, Ny, Nlc) + 1.0 ; Default unity
      if use_gain then begin
        gstate = selstates[0]
        gstate.tun_wavelength = gtun_wavelength[minloc] * 1e-10
        gstate.tuning = gtuning[minloc]
        gstate.fpi_state = gtuning[minloc]
        for ilc = 0, Nlc-1 do begin
          gstate.lc = ulc[ilc]
          self -> get_calib, gstate $
                             , gaindata = gg, gainstatus  = gainstatus
          if gainstatus ne 0 then stop
          gains[*, *, ilc] = gg
        endfor
      endif 

      
      iloop = 0
      Nloop = Nlp*Nqw*Nlc
      for ilp = 0, Nlp - 1 do begin
        for iqw = 0, Nqw - 1 do begin
          for ilc = 0, Nlc-1 do begin

            undefine, fullstate_list
            red_append, fullstate_list, 'lp'+string(round(ulp[ilp]), format = '(i03)')
            red_append, fullstate_list, 'qw'+string(round(uqw[iqw]), format = '(i03)')
            red_append, fullstate_list, upref[ipref]
            red_append, fullstate_list, '*lc'+strtrim(long(ulc[ilc]), 2)
            statestring = '*'+strjoin(fullstate_list, '_')
            
            indx = where(strmatch(selstates.fullstate,statestring ),count)
            if count ne 1 then begin
              print, inam + ' : ERROR, irregular state -> '+ statestring
              stop
            endif 
            red_progressbar, iloop, Nloop, /predict, cams[icam]+' : '+file_basename(selfiles[indx])
            d[ilc,iqw,ilp,*,*] = red_readdata(selfiles[indx], /silent) - dd
            if dodescatter then $
               d[ilc,iqw,ilp,*,*] = rdx_descatter(reform(d[ilc,iqw,ilp,*,*]) $
                                                  , bg, psf, nthreads = nthreads) * gains[*, *, ilc]
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

    endfor                      ; icam
  endfor                        ; ipref
  
end
