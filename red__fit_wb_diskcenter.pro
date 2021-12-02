; docformat = 'rst'

;+
; Measure median disk center wideband intensities and fit to function
; of time.
;
; The use of the median means it is OK to include data with pores or
; even a small spot in the FOV.
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
;    2019-07-15 : MGL. First version.
; 
;    2020-01-15 : MGL. New keywords tmin and tmax. If mu unavailable
;                 in logs, assume flats are at large enough mu.
; 
;    2020-01-20 : MGL. Use data from every 20th scan. 
; 
;    2021-08-27 : MGL. New keyword exclude_dirs. 
; 
;    2021-09-06 : MGL. New keyword limb_darkening. 
; 
;-

pro red::fit_wb_diskcenter, dirs = dirs $
                            , exclude_dirs = exclude_dirs $
                            , fitexpr = fitexpr_in $
                            , limb_darkening = limb_darkening $
                            , mu_limit = mu_limit $
                            , pref = pref $
                            , sum_all_frames = sum_all_frames $
                            , tmin = tmin $
                            , tmax = tmax
  
  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)

  if keyword_set(limb_darkening) and pref ne '8542' then begin
    if n_elements(mu_limit) eq 0 then mu_limit = 0.50d    
  endif else begin
    if n_elements(mu_limit) eq 0 then mu_limit = 0.97d
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
  camNB = (cams[where(strmatch(cams,'*-[NT]'))])[0]

  if n_elements(dirs) eq 0 then begin

    ;; Directories not provided, consider all.

    if ptr_valid(self.flat_dir)  then red_append, dirs, *self.flat_dir
    if ptr_valid(self.data_dirs) then red_append, dirs, *self.data_dirs

  endif

  if n_elements(exclude_dirs) gt 0 then begin

    match2, file_basename(dirs), file_basename(exclude_dirs), suba, subb
    dirs = dirs[where(suba eq -1)] 

  endif
  
  indx = where(strmatch(dirs,'*??:??:??'), Nwhere)
  if Nwhere eq 0 then stop
  dirs = dirs[indx]
  
  Ndirs = n_elements(dirs)
  if Ndirs eq 0 then stop

  ;; First get mu and zenith angle
  timeregex = '[0-2][0-9]:[0-5][0-9]:[0-6][0-9]'
  times = dblarr(Ndirs)
  for idir = 0, Ndirs-1 do begin
    times[idir] = red_time2double(stregex(dirs[idir], timeregex, /extract))
  endfor                        ; idir
  red_logdata, self.isodate, times, mu = mu, zenithangle = za
  ;; If mu isn't available in logs, assume flats are at large enough mu
;  indx = where(~finite(mu)  and strmatch(dirs,'*lats*'), Ndirs)
;  if Ndirs gt 0 then mu[indx] = mu_limit+1d-3

  
  ;; Idea: for flats directories, mu might vary during the data
  ;; collection. So find the data where mu peaks!
  
  ;; Is mu large enough? 
  indx = where(mu gt mu_limit, Ndirs)
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

  
  wbindx = where(file_test(dirs+'/'+camwb), Nwb)
  nbindx = where(file_test(dirs+'/'+camnb), Nnb)
  if Nwb eq 0 then noabsunits = 1
  
  if Nwb gt 0 then begin

    ;; Measure DC WB intensities

    wbdir = self.out_dir+'/wb_intensities/'
    file_mkdir, wbdir
    wbdirs = dirs[wbindx]
    for iwb = 0, Nwb-1 do begin
      print, wbdirs[iwb]

      ;; Read every 20th scan
      for iscan = 0L, 1000, 20 do begin

        ;;fnamesW = file_search(wbdirs[iwb]+'/'+camwb+'/*[._][0-9][0-9][0-9][02468]0[._]*', count = NfilesW)      
        ;;fnamesW = file_search(wbdirs[iwb]+'/'+camwb+'/*[._]'+string(iscan, format = '(i05)')+'[._]*' $
        ;;                    , count = NfilesW)      
        fnamesW = red_raw_search(wbdirs[iwb]+'/'+camwb, scannos = iscan, count = NfilesW)
        print, wbdirs[iwb]+'/'+camwb
        print, NfilesW
        if NfilesW eq 0 then break

        self -> extractstates, fnamesW, states, /nondb 

        pindx = uniq(states.prefilter, sort(states.prefilter))
        upref = states[pindx].prefilter
        Npref = n_elements(upref)

        for ipref = 0, Npref-1 do begin
          
          self -> get_calib, states[pindx[ipref]], darkdata = dark, gaindata = gain

          
          if keyword_set(sum_all_frames) then begin

            ;; Not for CHROMIS data!

            aindx = where(states.prefilter eq upref[ipref])
            Nframes = n_elements(aindx)
            imeans = fltarr(Nframes)
            tmeans = dblarr(Nframes)
            maskindx = where(gain ne 0)
            for iframe = 0, Nframes-1 do begin
              ims = red_readdata(fnamesW[aindx[iframe]], head = hdr, /silent)
              ;;  help, ims
              ims = ims[*, *, 0] - dark ; Select first frame in multi-frame file.
              ;; Descatter...
              imeans[iframe] = mean(ims[maskindx])
              red_fitspar_getdates, hdr, date_beg = date_beg
              tmeans[iframe] = red_time2double((strsplit(date_beg, 'T', /extract))[1]) 
            endfor              ; iframe

            red_append, wbintensity, max(imeans) ; This is how momfbd does it
            red_append, wbtimes, mean(tmeans)
            
            
          endif else begin
            
            ims = red_readdata(fnamesW[pindx[ipref]], head = hdr, /silent)
            
            help, ims
            ims = ims[*, *, 0] - dark ; Select first frame in multi-frame file.
            
            print, fnamesW[pindx[ipref]]          

            if self.dodescatter and (states[pindx[ipref]].prefilter eq '8542' $
                                     or states[pindx[ipref]].prefilter eq '7772') then begin
              self -> loadbackscatter, states[pindx[ipref]].detector $
                                       , states[pindx[ipref]].prefilter, bgain, bpsf
              ims = rdx_descatter(temporary(ims), bgain, bpsf, nthreads = nthread)
            endif
            
            red_fitspar_getdates, hdr, date_beg = date_beg
            

            if keyword_set(limb_darkening) then begin
              red_append, wbintensity, median(ims) / red_limb_darkening(float(upref[ipref])*1e-10, mu[wbindx[iwb]])
              red_append, wbintensity_orig, median(ims)
            endif else begin
              red_append, wbintensity, median(ims)
            endelse

            red_append, wbtimes, red_time2double((strsplit(date_beg, 'T', /extract))[1])
          endelse
          
          red_append, wbprefs, upref[ipref]
          red_append, wbexpt, fxpar(hdr, 'XPOSURE')
          red_append, wbmu, mu[wbindx[iwb]]
          red_append, wbza, za[wbindx[iwb]]
        endfor                  ; ipref
      endfor                    ; iscan
    endfor                      ; iwb
    
    ;; Sort
    indx = sort(wbtimes)
    wbintensity = wbintensity[indx]
    wbprefs = wbprefs[indx]
    wbtimes = wbtimes[indx]
    wbexpt = wbexpt[indx]
    wbmu = wbmu[indx]
    wbza = wbza[indx]
    if keyword_set(limb_darkening) then wbintensity_orig = wbintensity_orig[indx]
    
    upref = wbprefs[uniq(wbprefs, sort(wbprefs))]
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
      indx = where(wbprefs eq upref[ipref])
       
      if ipref eq 0 then begin
        ;; Set up the plot
        red_timeplot, /add $
                      , psym = 16, color = colors[ipref] $
                      , wbtimes[indx] $
                      , wbintensity[indx] / wbexpt[indx] / 1000. $
                      , xtitle = 'time [UT]' $
                      , ytitle = 'WB median intensity/exp time [counts/ms]' $
                      , xrange = [min(wbtimes), max(wbtimes)] + [-1., 1.]*60*20 $
                      , yrange = [0, max(wbintensity/wbexpt)*1.05]/1000.
      endif else begin
        cgplot, /add, /over, psym = 16, color = colors[ipref] $
                , wbtimes[indx] $
                , wbintensity[indx] / wbexpt[indx] / 1000.
      endelse

      if keyword_set(limb_darkening) then $
         cgplot, /add, /over, psym = 9, color = colors[ipref] $
                 , wbtimes[indx] $
                 , wbintensity_orig[indx] / wbexpt[indx] / 1000.
      
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
      
      if fitexpr_used[ipref] ne '' then begin

        ;; Do the fit
        err = 1. / wbmu[indx]   ; Large mu is better!
        pp = mpfitexpr(fitexpr_used[ipref] $
                       , wbtimes[indx]/3600. $
                       , wbintensity[indx] / wbexpt[indx] $
                       , err $
                      )

        ;; Plot it
        tt = (findgen(1000)/1000 * (max(wbtimes[indx])-min(wbtimes[indx])) + min(wbtimes[indx])) 
        wbint = red_evalexpr(fitexpr_used[ipref], tt/3600, pp)
        cgplot, /add, /over, color = colors[ipref], tt, wbint / 1000.

        ;; Save it
        red_mkhdr, phdr, pp
        anchor = 'DATE'
        red_fitsaddkeyword, anchor = anchor, phdr $
                            , 'FITEXPR', fitexpr_used[ipref], 'mpfitexpr fitting function'
        red_fitsaddkeyword, anchor = anchor, phdr $
                            , 'TIME-BEG', red_timestring(min(wbtimes[indx])) $
                            , 'Begin fit data time interval'
        red_fitsaddkeyword, anchor = anchor, phdr $
                            , 'TIME-END', red_timestring(max(wbtimes[indx])) $
                            , 'End fit data time interval'
        writefits, wbdir+'wb_fit_'+upref[ipref]+'.fits', pp, phdr

        coeffs_str[ipref] = strjoin(strtrim(pp, 2), ',')
        
      endif

      ;; Store the data
      openw, lun, /get_lun, wbdir+'wb_calibration_'+upref[ipref]+'.txt'
      printf, lun, '# time, intensity, exp time, mu, zenith angle'
      for i = 0, n_elements(indx)-1 do $
         printf, lun, wbtimes[indx[i]], wbintensity[indx[i]], wbexpt[indx[i]], wbmu[indx[i]], wbza[indx[i]]
      if keyword_set(limb_darkening) then begin
        printf, lun, ''        
        printf, lun, 'Without limb darkening correction:'        
        for i = 0, n_elements(indx)-1 do $
           printf, lun, wbtimes[indx[i]], wbintensity_orig[indx[i]], wbexpt[indx[i]], wbmu[indx[i]], wbza[indx[i]]
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
    if Nprefs eq 1 then $
       graph_name = 'wb_intensities_'+upref+'.pdf' $
    else $
       graph_name = 'wb_intensities.pdf'
    cgcontrol, output = wbdir+graph_name

  endif

  
end


;a = chromisred(/dev)
a = crispred(/dev)

a -> fit_wb_diskcenter
a -> fit_wb_diskcenter, pref = '4846'

pp=readfits('prefilter_fits/wb/wb_fit_6563.fits',h)

print, h, format = '(a0)'

END
