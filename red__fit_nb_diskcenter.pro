; docformat = 'rst'

;+
; Measure median disk center narrowband intensities and fit to function
; of time.
;
; Based on red::fit_wb_diskcenter.pro
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, ISP
; 
; 
; :Keywords:
;
;    dirs : in, optional, type=strarr
;
;        Set this to the time-stamp directories to use to limit the
;        automatic selection. Useful for excluding, e.g., data with a
;        large sunspot in the field of view.
;
;    exclude_dirs : in, optional, type=strarr
;
;        Set this to the time-stamp directories to exclude from
;        the automatically selected directories, or from the
;        directories specified with the dirs keyword.
;
;    fitexpr : in, optional, type="string or integer", default="depends on number of data points"
;
;        The expression that is fit to the intensities by use of
;        mpfitexpr. An integer is interpreted as shorthand for a
;        polynomial of that degree.
;
;    limb_darkening : in, optional, type=boolean
;
;        Correct intensities for distance from disk center.
;
;    mu_limit : in, optional, type=float, default="0.97 or 0.50"
;
;        WB data collected farther from disk center than this are
;        excluded automatically. The default depends on whether
;        correction for limb_darkening is in effect. 
;
;    pref : in, optional, type=string
;
;        Process data only for this prefilter.
;
;    tmax : in, optional, type=string
;
;        Do not use datasets started after this time, specified using
;        the format "HH[[:MM]:SS]". Useful if, e.g., there are both AM
;        and PM data, and the fit is not good if both are included.
; 
;    tmin : in, optional, type=string
;
;        Do not use datasets started before this time, specified using
;        the format "HH[[:MM]:SS]".
; 
; :History:
; 
;    2021-09-09 : MGL. First version.
; 
;    2021-10-19 : MGL. New keyword /demodulate.
; 
;-
pro red::fit_nb_diskcenter, demodulate = demodulate $
                            , dirs = dirs $
                            , exclude_dirs = exclude_dirs $
                            , fitexpr = fitexpr_in $
                            , limb_darkening = limb_darkening $
                            , mu_limit = mu_limit $
                            , pref = pref $
                            , tmin = tmin $
                            , tmax = tmax
  
  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)
  
  if keyword_set(limb_darkening) then begin 
    if n_elements(mu_limit) eq 0 then mu_limit = 0.50d
    if keyword_set(pref) then begin
      if pref eq '8542' then  begin        
        if n_elements(mu_limit) ne 0 and mu_limit lt 0.97d then begin
          print, inam, ' : For 8542 prefilter mu_limit should be >= 0.97'
          return
        endif
        mu_limit = 0.97d
      endif
    endif
  endif else begin
    if n_elements(mu_limit) eq 0 then begin
      mu_limit = 0.97d
    endif else begin
      if mu_limit lt 0.97d then begin
        print, inam, ' : You should use /limb_darkening with mu_limit < 0.97'
        return
      endif
    endelse
  endelse
    
  if n_elements(tmin) eq 0 then tmin = '00:00:00'
  if n_elements(tmax) eq 0 then tmax = '24:00:00'

  red_make_prpara, prpara, dirs
  red_make_prpara, prpara, exclude_dirs
  red_make_prpara, prpara, fitexpr_in
  red_make_prpara, prpara, limb_darkening
  red_make_prpara, prpara, mu_limit
  red_make_prpara, prpara, pref
  red_make_prpara, prpara, tmin  
  red_make_prpara, prpara, tmax  

  cams = *self.cameras
  camWB = (cams[where(strmatch(cams,'*-W'))])[0]
  camsNB = cams[where(strmatch(cams,'*-[NTR]'))]
  camNB = camsNB[0]

  instrument = strlowcase((strsplit(cams[0],'-',/extract))[0])
  
  if keyword_set(demodulate) then begin
    camT = 'Crisp-T'    
    camR = 'Crisp-R'    
  endif

  if n_elements(dirs) eq 0 then begin

    ;; Directories not provided, consider all.

    if ptr_valid(self.flat_dir)  then red_append, dirs, *self.flat_dir
    if ptr_valid(self.data_dirs) then red_append, dirs, *self.data_dirs

  endif else begin
    match2, dirs, file_basename(*self.data_dirs), suba, subb 
    mindx = where(subb ne -1, Nwhere)
    if Nwhere eq 0 then begin
       print, inam + ' : No data directories match the dirs keyword.'
       return
    endif
    dirs = (*self.data_dirs)[mindx]
  endelse

  if n_elements(exclude_dirs) gt 0 then begin

    match2, file_basename(dirs), file_basename(exclude_dirs), suba, subb
    dirs = dirs[where(suba eq -1)] 

  endif

  indx = where(strmatch(dirs,'*??:??:??'), Nwhere)
  if Nwhere eq 0 then stop
  dirs = dirs[indx]
  
  Ndirs = n_elements(dirs)
  if Ndirs eq 0 then stop

  
  ;; Figure out what prefilters are available and what tunings to use
  ;; by looking at fitprefilter output.

  if n_elements(pref) eq 0 then searchpref = '????' else searchpref = pref
  
  case instrument of
    'crisp'   : pfiles = file_search('prefilter_fits/'+camnb+'_'+searchpref+'_prefilter.idlsave', count = Nprefs)
    'chromis' : pfiles = file_search('prefilter_fits/chromis_'+searchpref+'_prefilter.idlsave', count = Nprefs)
    else : stop
  endcase

  pprefs = strarr(Nprefs)
  for ip = 0, Nprefs-1 do pprefs[ip] = (strsplit(file_basename(pfiles[ip]),'_',/extract))[1]
  if where(strmatch(pprefs, '8542')) ne -1 then begin
    if mu_limit lt 0.97d then begin
      print, inam, ' : You are using mu_limit < 0.97 with 8542 prefilter.'
      print,'Please change settings and rerun the program.'
      return
    endif
  endif
  ptunings = strarr(Nprefs)

  for ip = 0, Nprefs-1 do begin
    
    restore, pfiles[ip]  
    
    tmp = min(abs(red_time2double(file_basename(dirs)) - prf.time_avg),minloc)
    cdir = dirs[minloc]

    cfiles = self -> raw_search(cdir+'/'+camnb, pref = pprefs[ip], scannos = 0)
    self -> extractstates, cfiles, cstates
    tunindx = uniq(cstates.tun_wavelength, sort(cstates.tun_wavelength))
    utunwvl = cstates[tunindx].tun_wavelength
    ufpi = cstates[tunindx].fpi_state

    tmp = max(prf.spec/prf.pref, mloc)
    tmp = min(abs(utunwvl - prf.wav[mloc]*1e-10), minloc)

    ptunings[ip] = ufpi[minloc]

    cgwindow
    cgplot, /add, prf.wav/10., prf.spec/prf.pref, psym = 16, color = 'blue' $
            , title = self.isodate + ' ' + ptunings[ip] + ' ' + camnb $
            , xtitle = '$\lambda$ / 1 nm' $
            , ytitle = 'Intensity'
    
    cgplot, /add, /over, prf.wav[mloc]/10., prf.spec[mloc]/prf.pref[mloc], psym = 9, symsize = 2, color = 'red'

  endfor                        ; ip

  ;; CRISP tunings needs zero padding in raw dir!
  if instrument eq 'crisp' then begin
    for ipref = 0, Nprefs-1 do begin
      tunsplt = strsplit(ptunings[ipref], '_', /extract)
      tunlength = strlen(tunsplt[1]) ; length of the [+-]123 part
      if tunlength lt 5 then begin
        ;; Needs zero padding!
        tunsgn = strmid(tunsplt[1], 0, 1)
        ptunings[ipref] = tunsplt[0] + '_' + tunsgn + '*' + string(abs(tunsplt[1]), format = '(i03)')
      endif
    endfor                      ; ipref
  endif
  
  ;; Prefilter/tuning pairings are now in pprefs and ptunings.

  
  ;; Get mu and zenith angle
  timeregex = '[0-2][0-9]:[0-5][0-9]:[0-6][0-9]'
  times = dblarr(Ndirs)
  for idir = 0, Ndirs-1 do begin
    times[idir] = red_time2double(stregex(dirs[idir], timeregex, /extract))
  endfor                        ; idir
  red_logdata, self.isodate, times, mu = mu, zenithangle = za

;  ;; If mu isn't available in logs, assume flats are at large enough mu
;  indx = where(~finite(mu)  and strmatch(dirs,'*lats*'), Ndirs)  
;  if Ndirs gt 0 then mu[indx] = mu_limit+1d-3

  
  ;; Is mu large enough? 
  indx = where(mu gt mu_limit, Ndirs) ; Excludes mu=NaN
  if Ndirs eq 0 then begin
    print, inam + ' : No WB data with mu > ' + strtrim(mu_limit, 2)
    return
  endif else begin
    dirs = dirs[indx]
    times = times[indx]
    mu = mu[indx]
    za = za[indx]
  endelse
  
  ;; In specified time range?
  indx = where(times gt red_time2double(tmin) and times lt red_time2double(tmax), Ntime)
  if Ntime eq 0 then begin
    print, inam + ' : No WB data in time range ['+tmin+','+tmax+']'
    return
  endif else begin
    dirs = dirs[indx]
    times = times[indx]
    mu = mu[indx]
    za = za[indx]
  endelse 

  
;  wbindx = where(file_test(dirs+'/'+camwb), Nwb)
  nbindx = where(file_test(dirs+'/'+camnb), Nnb)
  
  if Nnb gt 0 then begin

    ;; Measure DC NB intensities

    nbdir = self.out_dir+'/nb_intensities/'
    if keyword_set(tmax) then $
      if tmax le '13:00:00' then $
        nbdir += 'morning/'
    if keyword_set(tmin) then $
      if tmin ge '13:00:00' then $
        nbdir += 'afternoon/'
    file_mkdir, nbdir
    
    nbdirs = dirs[nbindx]

    for ipref = 0, Nprefs-1 do begin

      print, pprefs[ipref]

      if keyword_set(demodulate) then begin
        file_mkdir, self.out_dir + '/nb_intensities/'+pprefs[ipref]+'/'
        self -> inverse_modmatrices, pprefs[ipref] $
                                     , self.out_dir + '/nb_intensities/'+pprefs[ipref]+'/' $
                                     , camr = 'Crisp-R', immr = immr $
                                     , camt = 'Crisp-T', immt = immt
      endif
      
      for inb = 0, Nnb-1 do begin
        print, nbdirs[inb]
        
        ;; Read every 20th scan
        for iscan = 0L, 1000, 20 do begin

          if keyword_set(demodulate) then begin

            fnamesT = self -> raw_search(nbdirs[inb]+'/'+camt, scannos = iscan $
                                         , count = NfilesT $
                                         , pref = pprefs[ipref], tunings = ptunings[ipref])
            
            fnamesR = self -> raw_search(nbdirs[inb]+'/'+camr, scannos = iscan $
                                         , count = NfilesR $
                                         , pref = pprefs[ipref], tunings = ptunings[ipref])

            print, nbdirs[inb]+'/'+camt            
            print, nbdirs[inb]+'/'+camr            
            print, NfilesT, NfilesR
            if NfilesR ne NfilesT then stop

            if NfilesT eq 0 then break
            
            self -> extractstates, fnamesT, statesT, /nondb             
            self -> extractstates, fnamesR, statesR, /nondb             

            ulc = statesT[uniq(statesT.lc, sort(statesT.lc))]
            Nlc = n_elements(ulc)

            if Nlc ne 4 then begin
              print, inam + ' : Non-polarimetric data, please call without /demodulate'
              print, pprefs[ipref]
              retall
            endif

            self -> get_calib, statesT[0], darkdata = darkT            
            self -> get_calib, statesR[0], darkdata = darkR          

            if self.dodescatter and (statesT[0].prefilter eq '8542' $
                                     or statesT[0].prefilter eq '7772') then begin
              self -> loadbackscatter, statesT[0].detector $              
                                       , statesT[0].prefilter, bgainT, bpsfT              
              self -> loadbackscatter, statesR[0].detector $              
                                       , statesR[0].prefilter, bgainR, bpsfR              
            endif
            
            ;; Demodulate

            ims_lcT = fltarr(1024, 1024, Nlc) ;; CRISP Sarnoffs!                        
            ims_lcR = fltarr(1024, 1024, Nlc) ;; CRISP Sarnoffs!                        
            
            indx0 = where(statesT.lc eq 0)         
            indx1 = where(statesT.lc eq 1)          
            indx2 = where(statesT.lc eq 2)          
            indx3 = where(statesT.lc eq 3)          

            ims_lcT[*,*,0] = (red_readdata(statesT[indx0[0]].filename, header = h0T) - darkT) ;* gain0          
            ims_lcT[*,*,1] = (red_readdata(statesT[indx1[0]].filename, header = h1T) - darkT) ;* gain1         
            ims_lcT[*,*,2] = (red_readdata(statesT[indx2[0]].filename, header = h2T) - darkT) ;* gain2          
            ims_lcT[*,*,3] = (red_readdata(statesT[indx3[0]].filename, header = h3T) - darkT) ;* gain3          

            ims_lcR[*,*,0] = (red_readdata(statesR[indx0[0]].filename, header = h0R) - darkR) ;* gain0          
            ims_lcR[*,*,1] = (red_readdata(statesR[indx1[0]].filename, header = h1R) - darkR) ;* gain1         
            ims_lcR[*,*,2] = (red_readdata(statesR[indx2[0]].filename, header = h2R) - darkR) ;* gain2          
            ims_lcR[*,*,3] = (red_readdata(statesR[indx3[0]].filename, header = h3R) - darkR) ;* gain3          
            
            if self.dodescatter and (pprefs[ipref] eq '8542' $
                                     or pprefs[ipref] eq '7772') then begin
              ims_lcT[*,*,0] = rdx_descatter(ims_lcT[*,*,0], bgainT, bpsfT, nthreads = nthread)            
              ims_lcT[*,*,1] = rdx_descatter(ims_lcT[*,*,1], bgainT, bpsfT, nthreads = nthread)            
              ims_lcT[*,*,2] = rdx_descatter(ims_lcT[*,*,2], bgainT, bpsfT, nthreads = nthread)            
              ims_lcT[*,*,3] = rdx_descatter(ims_lcT[*,*,3], bgainT, bpsfT, nthreads = nthread)            

              ims_lcR[*,*,0] = rdx_descatter(ims_lcR[*,*,0], bgainR, bpsfR, nthreads = nthread)            
              ims_lcR[*,*,1] = rdx_descatter(ims_lcR[*,*,1], bgainR, bpsfR, nthreads = nthread)            
              ims_lcR[*,*,2] = rdx_descatter(ims_lcR[*,*,2], bgainR, bpsfR, nthreads = nthread)            
              ims_lcR[*,*,3] = rdx_descatter(ims_lcR[*,*,3], bgainR, bpsfR, nthreads = nthread)            

            endif

            red_fitspar_getdates, h0T, date_avg = date_avg0          
            red_fitspar_getdates, h1T, date_avg = date_avg1          
            red_fitspar_getdates, h2T, date_avg = date_avg2          
            red_fitspar_getdates, h3T, date_avg = date_avg3          

            time_avg0 = red_time2double((strsplit(date_avg0, 'T', /extract))[1])          
            time_avg1 = red_time2double((strsplit(date_avg1, 'T', /extract))[1])          
            time_avg2 = red_time2double((strsplit(date_avg2, 'T', /extract))[1])          
            time_avg3 = red_time2double((strsplit(date_avg3, 'T', /extract))[1])

            time_avg = (time_avg0+time_avg1+time_avg2+time_avg3)/4.

            
            
            ims_stokesT = red_demodulate_images(ims_lcT, reform(immt, [4, 4, 1024,1024]), self.isodate, pprefs[ipref], time_avg)            
            ims_stokesR = red_demodulate_images(ims_lcR, reform(immr, [4, 4, 1024,1024]), self.isodate, pprefs[ipref], time_avg)            

            ;; Should really use the projective transform here! But
            ;; misalignment errors should average out in the median
            ;; intensity. 
            
            imsT = ims_stokesT[*, *, 0] ; Stokes I           
            imsR = ims_stokesR[*, *, 0]             

            mediant = median(imsT)
            medianr = median(imsR)
            aver = (mediant + medianr) / 2.
            sct = aver / mediant
            scr = aver / medianr
            
            ims = (sct * (imsT) + scr * (imsR)) / 2.

            date_beg = self.isodate+'T'+red_timestring(time_avg)
            hdr = h0T

          endif else begin
            fnamesN = self -> raw_search(nbdirs[inb]+'/'+camnb, scannos = iscan $
                                         , count = NfilesN $
                                         , pref = pprefs[ipref], tunings = ptunings[ipref])
            print, nbdirs[inb]+'/'+camnb
            print, NfilesN
            if NfilesN eq 0 then break

            self -> extractstates, fnamesN, states, /nondb 

;            ulc = states[uniq(states.lc, sort(states.lc))]
;            Nlc = n_elements(ulc)
            
;        pindx = uniq(states.prefilter, sort(states.prefilter))
;        upref = states[pindx].prefilter
;        Npref = n_elements(upref)

            
            self -> get_calib, states, darkdata = dark, gaindata = gain
;          if keyword_set(sum_all_frames) then begin
;
;            ;; Not for CHROMIS data!
;
;            aindx = where(states.prefilter eq upref[ipref])
;            Nframes = n_elements(aindx)
;            imeans = fltarr(Nframes)
;            tmeans = dblarr(Nframes)
;            maskindx = where(gain ne 0)
;            for iframe = 0, Nframes-1 do begin
;              ims = red_readdata(fnamesW[aindx[iframe]], head = hdr, /silent)
;              ;;  help, ims
;              ims = ims[*, *, 0] - dark ; Select first frame in multi-frame file.
;              ;; Descatter...
;              imeans[iframe] = mean(ims[maskindx])
;              red_fitspar_getdates, hdr, date_beg = date_beg
;              tmeans[iframe] = red_time2double((strsplit(date_beg, 'T', /extract))[1]) 
;            endfor              ; iframe
;
;            red_append, wbintensity, max(imeans) ; This is how momfbd does it
;            red_append, wbtimes, mean(tmeans)
;            
;            
;          endif else begin
            
            ims = red_readdata(fnamesN[0], head = hdr, /silent) ; Select first file available
            
            help, ims
            ims = ims[*, *, 0] - dark ; Select first frame in multi-frame file.
            
            print, fnamesN[0]

            if self.dodescatter and (states[0].prefilter eq '8542' $
                                     or states[0].prefilter eq '7772') then begin
              self -> loadbackscatter, states[0].detector $
                                       , states[0].prefilter, bgain, bpsf
              ims = rdx_descatter(temporary(ims), bgain, bpsf, nthreads = nthread)
            endif
            
            red_fitspar_getdates, hdr, date_beg = date_beg

          endelse

          if keyword_set(limb_darkening) then begin
            red_append, nbintensity, median(ims) / red_limb_darkening(float(pprefs[ipref])*1e-10, mu[nbindx[inb]])
            red_append, nbintensity_orig, median(ims)
          endif else begin
            red_append, nbintensity, median(ims)
          endelse

          red_append, nbtimes, red_time2double((strsplit(date_beg, 'T', /extract))[1])
;          endelse
          
          red_append, nbprefs, pprefs[ipref]
          red_append, nbexpt, fxpar(hdr, 'XPOSURE')
          red_append, nbmu, mu[nbindx[inb]]
          red_append, nbza, za[nbindx[inb]]
        endfor                  ; iscan
      endfor                    ; inb
    endfor                      ; ipref
    
    ;; Sort
    indx = sort(nbtimes)
    nbintensity = nbintensity[indx]
    nbprefs = nbprefs[indx]
    nbtimes = nbtimes[indx]
    nbexpt = nbexpt[indx]
    nbmu = nbmu[indx]
    nbza = nbza[indx]
    if keyword_set(limb_darkening) then nbintensity_orig = nbintensity_orig[indx]
    
    upref = nbprefs[uniq(nbprefs, sort(nbprefs))]
    if n_elements(pref) ne 0 then begin
      ;; Limit to specified prefilters 
      match2, upref, pref, suba
      mindx = where(suba ne -1, Nmatch)
      if Nmatch ne 0 then upref = upref[mindx]
    endif
    
    Nprefs = n_elements(upref)
    
    fitexpr_used = strarr(Nprefs)
    coeffs_str = strarr(Nprefs)
    
    ;; Prepare for plotting the results
    ;;colors = ['blue', 'red', 'green', 'plum', 'cyan', 'darkkhaki']
    ;; Colors from cgPickColorName()
    colors = ['RED', 'BLU', 'GRN', 'ORG', 'PUR', 'YGB', 'PBG', 'BLK'] + '5'
    if Nprefs le n_elements(colors) then begin
      colors = colors[0:Nprefs-1]
    endif else begin
      ;; Not likely to be needed:
      colors = red_wavelengthtorgb(float(upref)/10., /num)
      ;;colors = distinct_colors(n_colors = Nprefs, /num)
    endelse
    cgwindow
    
    ;; Loop over prefilters
    for ipref = 0, Nprefs-1 do begin
      indx = where(nbprefs eq pprefs[ipref])
       
      if ipref eq 0 then begin
        ;; Set up the plot
        red_timeplot, /add $
                      , psym = 16, color = colors[ipref] $
                      , nbtimes[indx] $
                      , nbintensity[indx] / nbexpt[indx] / 1000. $
                      , xtitle = 'time [UT]' $
                      , ytitle = 'NB median intensity/exp time [counts/ms]' $
                      , xrange = [min(nbtimes), max(nbtimes)] + [-1., 1.]*60*20 $
                      , yrange = [0, max(nbintensity/nbexpt)*1.05]/1000.
      endif else begin
        cgplot, /add, /over, psym = 16, color = colors[ipref] $
                , nbtimes[indx] $
                , nbintensity[indx] / nbexpt[indx] / 1000.
      endelse

      if keyword_set(limb_darkening) then $
         cgplot, /add, /over, psym = 9, color = colors[ipref] $
                 , nbtimes[indx] $
                 , nbintensity_orig[indx] / nbexpt[indx] / 1000.
      
      ;; Set up for fitting
      if n_elements(fitexpr_in) eq 0 then begin
        ;; Defaults based on range of data
        case n_elements(indx) of
          1    : fitexpr_used[ipref] = ''
          2    : fitexpr_used[ipref] = 'P[0] + X*P[1]'
          else : fitexpr_used[ipref] = 'P[0] + X*P[1] + X*X*P[2]'
        endcase
      endif else if size(fitexpr_in, /tname) eq 'STRING' then begin
        ;; Use an actual string
        fitexpr_used[ipref] = fitexpr_in
      endif else begin
        ;; Assume integer if not string
        case fitexpr_in of
          1 : fitexpr_used[ipref] = 'P[0] + X*P[1]'
          2 : fitexpr_used[ipref] = 'P[0] + X*P[1] + X*X*P[2]'          
          3 : fitexpr_used[ipref] = 'P[0] + X*P[1] + X*X*P[2] + X*X*X*P[3]'          
          4 : fitexpr_used[ipref] = 'P[0] + X*P[1] + X*X*P[2] + X*X*X*P[3] + X*X*X*X*P[4]'                    
          5 : fitexpr_used[ipref] = 'P[0] + X*P[1] + X*X*P[2] + X*X*X*P[3] + X*X*X*X*P[4] + X*X*X*X*X*P[5]'
          else : stop
        endcase
      endelse
      

;      fitexpr_used[ipref] = fitexpr
      
      if fitexpr_used[ipref] ne '' then begin

        ;; Do the fit
        err = 1. / nbmu[indx]   ; Large mu is better!
        pp = mpfitexpr(fitexpr_used[ipref] $
                       , nbtimes[indx]/3600. $
                       , nbintensity[indx] / nbexpt[indx] $
                       , err $
                      )

        ;; Plot it
        tt = (findgen(1000)/1000 * (max(nbtimes[indx])-min(nbtimes[indx])) + min(nbtimes[indx])) 
        nbint = red_evalexpr(fitexpr_used[ipref], tt/3600, pp)
        cgplot, /add, /over, color = colors[ipref], tt, nbint / 1000.

        

        if 0 then begin         ; 6302

          ;; Compare with stokes data
          stokfiles = file_search('momfbd_nopd/??:??:??/6302/cfg/results/stokes_sbs0_sbk5/stokesIQUV_00000_6302_6302_-290.fits', count = Nstokes)

          stokint = fltarr(Nstokes)
          stoktime = dblarr(Nstokes)
          for istok = 0, Nstokes-1 do begin
            red_fitscube_getframe, stokfiles[istok], stokim, istokes = 0
            stokint[istok] = median(stokim)
            stokhdr = headfits(stokfiles[istok])
            red_fitspar_getdates, stokhdr, date_avg = date_avg
            stoktime[istok] = red_time2double((strsplit(date_avg,'T',/extract))[1])
          endfor

          cubfiles = file_search('cubes_scan_none/nb_6302_2016-09-19T??:??:??_scan=0_stokes_corrected.fits', count = Ncub)
          cubint = fltarr(Ncub)
          cubtime = dblarr(Ncub)
          for icub = 0, Ncub-1 do begin
            red_fitscube_getframe, cubfiles[icub], cubim, istokes = 0, itun = 9
            cubint[icub] = median(cubim)
            cubhdr = headfits(cubfiles[icub])
            red_fitspar_getdates, cubhdr, date_avg = date_avg
            cubtime[icub] = red_time2double((strsplit(date_avg,'T',/extract))[1])
          endfor

          ;; WB camera is camXX so not included here:
          momfiles = file_search('momfbd_nopd/??:??:??/6302/cfg/results/cam???_*_00000_*_-290_*.momfbd', count = Nmom)
          momint = fltarr(Nmom)
          momtime = dblarr(Nmom)
          for imom = 0, Nmom-1 do begin
            momim = red_readdata(momfiles[imom], h = momhdr)
            momint[imom] = median(momim)
            red_fitspar_getdates, momhdr, date_avg = date_avg
            momtime[imom] = red_time2double((strsplit(date_avg,'T',/extract))[1])
          endfor



          
          red_timeplot,nbtimes,nbintensity/max(nbintensity),psym=16, yrange =[0, 1.1], trange = [8, 11.5]*3600
          cgplot, /over, momtime, momint/max(momint), psym = 2, color = 'green'
          cgplot, /over, stoktime, stokint/max(stokint), psym = 16, color = 'red'          
          cgplot, /over, cubtime, cubint/max(cubint), psym = 16, color = 'blue'
          cgplot, /over, tt, nbint / max(nbint)

          stop
          
        endif
        
        
        
        ;; Save it
        red_mkhdr, phdr, pp
        anchor = 'DATE'
        red_fitsaddkeyword, anchor = anchor, phdr $
                            , 'FITEXPR', fitexpr_used[ipref], 'mpfitexpr fitting function'
        red_fitsaddkeyword, anchor = anchor, phdr $
                            , 'TIME-BEG', red_timestring(min(nbtimes[indx])) $
                            , 'Begin fit data time interval'
        red_fitsaddkeyword, anchor = anchor, phdr $
                            , 'TIME-END', red_timestring(max(nbtimes[indx])) $
                            , 'End fit data time interval'
        
        ;; Add info about this step
        self -> headerinfo_addstep, phdr $
                                    , prstep = 'CALIBRATION' $
                                    , prpara = prpara $
                                    , prproc = inam

        writefits, nbdir+'nb_fit_'+upref[ipref]+'.fits', pp, phdr

        coeffs_str[ipref] = strjoin(strtrim(pp, 2), ',')
        
      endif

      ;; Store the data
      openw, lun, /get_lun, nbdir+'nb_calibration_'+upref[ipref]+'.txt'
      printf, lun, '# time, intensity, exp time, mu, zenith angle'
      for i = 0, n_elements(indx)-1 do $
         printf, lun, nbtimes[indx[i]], nbintensity[indx[i]], nbexpt[indx[i]], nbmu[indx[i]], nbza[indx[i]]
      if keyword_set(limb_darkening) then begin
        printf, lun, ''        
        printf, lun, 'Without limb darkening correction:'        
        for i = 0, n_elements(indx)-1 do $
           printf, lun, nbtimes[indx[i]], nbintensity_orig[indx[i]], nbexpt[indx[i]], nbmu[indx[i]], nbza[indx[i]]
      endif
      free_lun, lun

    endfor                      ; ipref


    ;; Transform the fitted expression for printing
    for ipref = 0, Nprefs-1 do begin
      fitexpr_used[ipref] = strlowcase(fitexpr_used[ipref])
      fitexpr_used[ipref] = red_strreplace(fitexpr_used[ipref],'x*x*x*x*','x$\exp4$')
      fitexpr_used[ipref] = red_strreplace(fitexpr_used[ipref],'x*x*x*','x$\exp3$')
      fitexpr_used[ipref] = red_strreplace(fitexpr_used[ipref],'x*x*','x$\exp2$')
      fitexpr_used[ipref] = red_strreplace(fitexpr_used[ipref],'x*','x')            
      fitexpr_used[ipref] = red_strreplace(fitexpr_used[ipref],'p[','p$\sub',n=10)  
      fitexpr_used[ipref] = red_strreplace(fitexpr_used[ipref],']','$',n=10)        
    endfor                      ; ipref
    
    ;; Finish the plot
    cglegend, /add $
              , titles = upref + ' : ' + fitexpr_used $
              ;;, titles = upref + ' : ' + coeffs_str $
              , location = [0.9, 0.12], align = 2 $
              , colors = colors, psym = 16, length = 0, vspace = 2
    outfile = 'nb_intensities'
    for ipref=0,Nprefs-1 do outfile += '_' + upref[ipref]
    outfile += '.pdf'
    if file_test(nbdir+outfile) then begin
      print,'File ', nbdir+outfile, ' exists.'
      print
      ans=''
      read,'Would you like to overwrite it? (Y/n): ', ans
      if strupcase(ans) eq 'N' then return
    endif
    cgcontrol, output = nbdir+outfile

  endif

  
end


;a = chromisred(/dev)
a = crispred(/dev)

a -> fit_nb_diskcenter
a -> fit_nb_diskcenter, pref = '4846'


END
