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
;   2024-11-21 : JdlCR. Implemented flatfielding of the modulation
;                matrix for alternative demodulation done through. To
;                revert to old behaviour use /no_polcal_flatfielding.
;
;-
pro red::polcalcube, cam = cam, pref = pref, no_descatter = no_descatter, nthreads = nthreads $
                     , no_polcal_flatfielding = no_polcal_flatfielding $
                     , force_polcal_flatfielding = force_polcal_flatfielding
  
  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)  

  ;; We started taking proper flats at the polcal wavelength in late
  ;; 2023. Flat-field the polcal data by default if the date is later
  ;; than that.
  if (self.isodate LT red_dates(tag = 'polcal flats')) $
     && ~keyword_set(force_polcal_flatfielding) then no_polcal_flatfielding = 1B

  ;; Alternatively, check wavelength of polcal data and see if we have
  ;; flats at that wavelength. Then base the default processing on
  ;; that! (This check is done later so we might be able to move that
  ;; up to here. Or to the beginning of the prefilter loop.)
  
  ;; Check polcal_sums
  if(~file_test(self.out_dir + '/polcal_sums', /directory)) then begin
    print, inam + ' : ERROR, folder not found : '+ self.out_dir + '/polcal_sums'
    stop
  endif

  search_str = 'polcal_sums/'
  if keyword_set(cam) then search_str += cam + '/' else search_str += '*/'
  if ~keyword_set(nthreads) then nthreads = 4L
  if keyword_set(pref) then $
     search_str += '*' + pref + '*.fits' $
  else $
      search_str += '*.fits'
  files = file_search(self.out_dir + search_str, count = count)
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
    
    dodescatter = ~keyword_set(no_descatter) AND self.dodescatter $
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
      if dodescatter then begin
        self -> loadbackscatter, detector, upref[ipref], bg, psf
      endif

      uqw = (selstates[uniq(selstates.qw, sort(selstates.qw))]).qw
      ulp = (selstates[uniq(selstates.lp, sort(selstates.lp))]).lp
      ulc = (selstates[uniq(selstates.lc, sort(selstates.lc))]).lc

      Nqw = n_elements(uqw)
      Nlp = n_elements(ulp)
      Nlc = n_elements(ulc)

      dim = fxpar(headfits(files[0]), 'NAXIS*')
      Nx = dim[0]
      Ny = dim[1]

      ;; Load dark
      gstate = selstates[0]
      self -> get_calib, gstate $
                         , darkdata = dd, darkstatus  = darkstatus $
                         , gainname = gn, flatname = fn
      if darkstatus ne 0 then stop

      gains = replicate(1., Nx, Ny, Nlc) ; Just unity until gain correction is sufficiently tested
      mask = fltarr(Nx,Ny,Nlc) + 1.0 ;;red_taper([Nx,round(Nx*(1-1/sqrt(2))/2.),0])
      
      flat_failed = 1
      
      if(~keyword_set(no_polcal_flatfielding) && (Nlc gt 1)) then begin
        flat_failed = 0
        if(file_test(fn)) then begin
          for ilc=0,Nlc-1 do begin
            self -> get_calib, selstates[ilc] $
                               , darkdata = dd, darkstatus  = darkstatus $
                               , gainname = gn, flatname = fn
            if(file_test(fn)) then begin
              print, inam + " : reading polcal flat-field file and making gain -> "+fn
              ff = red_readdata(fn)
              gains[*,*,ilc] = red_flat2gain(ff, bad=1.0, smooth=3.0, min=0.085, max=5.0)
              if((gains[0,0,ilc] eq 0) && $
                 (gains[0,-1,ilc] eq 0) && $
                 (gains[-1,0,ilc] eq 0) && $
                 (gains[-1,-1,ilc] eq 0)) then begin
                mask[*,*,ilc] = red_cleanmask(gains[*,*,ilc] eq 0, /circ)
              endif else mask[*,*,ilc] = 1.0
            endif else begin
              print, inam + "something is wrong, did we take flats in all LC states at polcal wavelength?"
              flat_failed = 1
              break 
            endelse
          endfor
        endif else flat_failed = 1
        
        if(flat_failed) then begin
        
          ;; No flat at polcal wavelength, check if there is another
;          ffiles = file_search('flats/' $
;                               + selstates[0].detector + '_' $
;                               + selstates[0].cam_settings + '_' $
;                               + selstates[0].prefilter + '_' $
;                               + '*' $
;                               + 'lc'+strtrim(long(selstates[0].lc), 2) $
;                                        + '.flat.fits' $
;                                        , count = Ng )
          search_str = self -> filenames('flat', selstates[0], /wild_tuning)
          ffiles = file_search(search_str, count = Ng)

          self -> extractstates, ffiles, fstates
          mn = min(abs(fstates.tun_wavelength - selstates[0].tun_wavelength),minloc)
;          red_extractstates, ffiles, dwav = gtun_wavelength, wav = gtuning
;          mn = min(abs(gtun_wavelength - selstates[0].tun_wavelength * 1e10),minloc)
          print
          print, inam + ' : There are no flat field data for the polcal wavelength.'
          print, '   Tuning for the polcal data: '+selstates[0].fpi_state
          print, '   Nearest flats tuning: '+fstates[minloc].fpi_state
          s = ''
          read, '    Do you want to use it to flat-field the polcal data [Y/n]? ', s
          
          if strlen(s) eq 0 || ~(strlowcase(strmid(s, 0, 1)) eq 'n') then begin
            gstate.tun_wavelength =fstates[minloc].tun_wavelength
            gstate.tuning = fstates[minloc].tuning
            gstate.fpi_state = fstates[minloc].fpi_state
;            gstate.tun_wavelength = gtun_wavelength[minloc] * 1e-10
;            gstate.tuning = gtuning[minloc]
;            gstate.fpi_state = gtuning[minloc]
            ;;gains = fltarr(Nx, Ny, Nlc)
            for ilc = 0, Nlc-1 do begin
              gstate.lc = ulc[ilc]
              self -> get_calib, gstate, gainname = gainname $
                                 , flatstatus  = flatstatus $
                                 , flatname = fn, flatdata = ff
              
              
              
              if flatstatus ne 0 then stop
              print, inam+" : loading flat data and making gain -> "+fn 
              gains[*, *, ilc] = red_flat2gain(ff, bad=1.5, smooth=3.0, min=0.085, max=5.0)
              if((gains[0,0,ilc] eq 0) && $
                 (gains[0,-1,ilc] eq 0) && $
                 (gains[-1,0,ilc] eq 0) && $
                 (gains[-1,-1,ilc] eq 0)) then begin
                mask[*,*,ilc] = red_cleanmask(gains[*,*,ilc] eq 0, /circ)
              endif else mask[*,*,ilc] = 1.0

            endfor              ; ilc              
          endif
        endif
      endif
      
      mask = total(mask,3) ne 0
      
;      if file_test(gn) then begin
;        print, inam + ' : Gains exist for the polcal wavelength, using them.'
;        gains = fltarr(Nx, Ny, Nlc)
;        for ilc = 0, Nlc-1 do begin
;          gstate.lc = ulc[ilc]
;          self -> get_calib, gstate, gainname = gainname $
;                             , gaindata = gg, gainstatus  = gainstatus
;          print, gainname
;          if gainstatus ne 0 then stop
;          gains[*, *, ilc] = gg
;        endfor
;        mask = total(gains,3) ne 0 
;      endif else begin
;        ;; Do we want to use gains from a nearby tuning? Check what
;        ;; tunings are available.
;        gfiles = file_search('gaintables/' $
;                             + selstates[0].detector + '_' $
;                             + selstates[0].cam_settings + '_' $
;                             + selstates[0].prefilter + '_' $
;                             + '*' $
;                             + 'lc'+strtrim(long(selstates[0].lc), 2) $
;                             + '.gain.fits' $
;                             , count = Ng )
;        red_extractstates, gfiles, dwav = gtun_wavelength, wav = gtuning
;        mn = min(abs(gtun_wavelength - selstates[0].tun_wavelength * 1e10),minloc)
;        if n_elements(use_gain) eq 0 then begin
;          print
;          print, inam + ' : There are no flat field data for the polcal wavelength.'
;          print, '   Tuning for the polcal data: '+selstates[0].tuning
;          print, '   Nearest flats tuning: '+gtuning[minloc]
;          s = ''
;          read, '    Do you want to use it to flat-field the polcal data [Y/n]? ', s
;          if strlen(s) eq 0 || ~(strlowcase(strmid(s, 0, 1)) eq 'n') then use_gain = 1 else use_gain = 0
;        endif
;        
;        if use_gain then begin
;          gstate.tun_wavelength = gtun_wavelength[minloc] * 1e-10
;          gstate.tuning = gtuning[minloc]
;          gstate.fpi_state = gtuning[minloc]
;          gains = fltarr(Nx, Ny, Nlc)
;          for ilc = 0, Nlc-1 do begin
;            gstate.lc = ulc[ilc]
;            self -> get_calib, gstate, gainname = gainname $
;                               , gaindata = gg, gainstatus  = gainstatus
;            if gainstatus ne 0 then stop
;            print, gainname
;            gains[*, *, ilc] = gg
;          endfor                ; ilc
;          mask = total(gains,3) ne 0
;        endif else begin
;          gains = replicate(1., Nx, Ny, Nlc) ; Just unity then
;          if Nx ne Ny then stop              ; Not implemented yet
;          mask = red_taper([Nx,round(Nx*(1-1/sqrt(2))/2.),0])
;        endelse
;
;      endelse

      ;; Read data
      d = fltarr(Nlc, Nqw, Nlp, Nx, Ny)
      d1d = fltarr(Nlc, Nqw, Nlp)
      
      total_mask = total(mask)
      
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
            tmp = red_readdata(selfiles[indx], /silent) - dd
            if dodescatter then $
               tmp = rdx_descatter(reform(tmp) , bg, psf, nthreads = nthreads) 
            
            tmp = red_fillpix(tmp * gains[*, *, ilc], nthreads=max([nthreads,4]), mask=mask)

            idx=where(tmp lt 0.00001, count, complement=idx1)
            if(count gt 0) then tmp[idx] = 0.0

            
            
            d[ilc,iqw,ilp,*,*] = tmp
            ;;d1d[ilc,iqw,ilp] = mean(d[ilc,iqw,ilp,100:Nx-101,100:Ny-101], /nan)
            d1d[ilc,iqw,ilp] = total(mask*reform(d[ilc,iqw,ilp,*, *]), /nan) / total_mask
            
            ;; To be completely correct, total_mask would need to take
            ;; the number of NaNs into account. Otherwise we would
            ;; either have to make new mask each time and this is deep
            ;; inside multiple loops. However, this is only an
            ;; initialization for optimization done by polcal so we
            ;; bet on that number to be small within the mask. We
            ;; could also use both the mask and median(), which treats
            ;; NaNs as missing data.
            
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

      mkhdr, ehdr, [flat_failed], /image
      red_fitsaddkeyword, anchor = anchor, ehdr, 'EXTNAME', 'PF', 'Polcal flatfielding'
      writefits, pname, [~flat_failed], ehdr, /append
      
      
    endfor                      ; icam
  endfor                        ; ipref

end
