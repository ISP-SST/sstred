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
; :Returns:
; 
; 
; :Params:
; 
; 
; :Keywords:
;
;    dirs : in, optional, type=string
;
;        Set this to the time-stamp directories to use to limit the
;        automatic selection.
;
;    fitexpr : in, optional, type=string, default="depends on number of data points"
;
;        The function fit to the intensities by use of mpfitexpr.
;
;    pref : in, optional, type=string
;
;        Process data only for this prefilter.
; 
; :History:
; 
;    2019-07-15 : MGL. First version.
; 
;-
pro red::fit_wb_diskcenter, dirs = dirs $
                            , fitexpr = fitexpr_in $
                            , mu_limit = mu_limit $
                            , pref = pref

  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)

  if n_elements(mu_limit) eq 0 then mu_limit = 0.97
  
  cams = *self.cameras
  camWB = (cams[where(strmatch(cams,'*-W'))])[0]
  camNB = (cams[where(strmatch(cams,'*-[NT]'))])[0]

  if n_elements(dirs) eq 0 then begin

    ;; Directories not provided, consider all.

    if ptr_valid(self.flat_dir)  then red_append, dirs, *self.flat_dir
    if ptr_valid(self.data_dirs) then red_append, dirs, *self.data_dirs

  endif
  
  Ndirs = n_elements(dirs)
  if Ndirs eq 0 then stop

  ;; First get mu and zenith angle
  timeregex = '[0-2][0-9]:[0-5][0-9]:[0-6][0-9]'
  times = dblarr(Ndirs)
  for idir = 0, Ndirs-1 do begin
    times[idir] = red_time2double(stregex(dirs[idir], timeregex, /extract))
  endfor                        ; idir
  red_logdata, self.isodate, times, mu = mu, zenithangle = za

  ;; Idea: for flats directories, mu might vary during the data
  ;; collection. So find the data where mu peaks!
  
  
  indx = where(mu gt mu_limit, Ndirs)
  if Ndirs eq 0 then stop
  dirs = dirs[indx]
  times = times[indx]
  mu = mu[indx]
  za = za[indx]

  wbindx = where(file_test(dirs+'/'+camwb), Nwb)
  nbindx = where(file_test(dirs+'/'+camnb), Nnb)
  if Nwb eq 0 then noabsunits = 1

  if Nwb gt 0 then begin

    ;; Measure DC WB intensities

    ;;wbdir = self.out_dir+'/prefilter_fits/wb/'
    wbdir = self.out_dir+'/wb_intensities/'
    file_mkdir, wbdir
    wbdirs = dirs[wbindx]
    for iwb = 0, Nwb-1 do begin
      print, wbdirs[iwb]
      ;;fnamesW = file_search(wbdirs[iwb]+'/Chromis-W/*_00000_*fits', count = NfilesW)      
      fnamesW = file_search(wbdirs[iwb]+'/'+camwb+'/*[._]00000[._]*', count = NfilesW)      
      self -> extractstates, fnamesw, states

      pindx = uniq(states.prefilter, sort(states.prefilter))
      upref = states[pindx].prefilter
      Npref = n_elements(upref)

      for ipref = 0, Npref-1 do begin

        self -> get_calib, states[pindx[ipref]], darkdata = dark
        ims = red_readdata(fnamesW[pindx[ipref]], head = hdr, /silent) - dark

        if self.dodescatter and (states[pindx[ipref]].prefilter eq '8542' $
                                 or states[pindx[ipref]].prefilter eq '7772') then begin
          self -> loadbackscatter, states[pindx[ipref]].detector $
                                   , states[pindx[ipref]].prefilter, bgain, bpsf
          ims = rdx_descatter(temporary(ims), bgain, bpsf, nthreads = nthread)
        endif
        
        
        

        red_fitspar_getdates, hdr, date_beg = date_beg
        
        red_append, wbintensity, median(ims)
        red_append, wbprefs, upref[ipref]
        red_append, wbtimes, red_time2double((strsplit(date_beg, 'T', /extract))[1])
        red_append, wbexpt, fxpar(hdr, 'XPOSURE')
        red_append, wbmu, mu[wbindx[iwb]]
        red_append, wbza, za[wbindx[iwb]]
      endfor                    ; ipref
    endfor                      ; iwb

    ;; Sort
    indx = sort(wbtimes)
    wbintensity = wbintensity[indx]
    wbprefs = wbprefs[indx]
    wbtimes = wbtimes[indx]
    wbexpt = wbexpt[indx]
    wbmu = wbmu[indx]
    wbza = wbza[indx]
    
    upref = wbprefs[uniq(wbprefs, sort(wbprefs))]
    Nprefs = n_elements(upref)
    
    fitexpr_used = strarr(Nprefs)
    coeffs_str = strarr(Nprefs)
    
    ;; Prepare for plotting the results
    if Nprefs le n_elements(colors) then begin
      colors = ['blue', 'red', 'green', 'plum', 'cyan', 'darkkhaki']
      ;; Colors from cgPickColorName()
      colors = ['RED', 'BLU', 'GRN', 'ORG', 'PUR', 'YGB', 'PBG', 'BLK'] + '5'
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

      if n_elements(pref) ne 0 then begin
        ;; Limit to specified prefilters 
        match2, upref[ipref], pref, suba
        mindx = where(suba ne -1, Nmatch)
        if Nmatch eq 0 then continue
      endif 

      if ipref eq 0 then begin
        ;; Set up the plot
        red_timeplot, /add $
                      , psym = 16, color = colors[ipref] $
                      , wbtimes[indx] $
                      , wbintensity[indx] / wbexpt[indx] / 1000. $
                      , xtitle = 'time [UT]', ytitle = 'WB median intensity/exp time [counts/ms]' $
                      , xrange = [min(wbtimes), max(wbtimes)] + [-0.1, 0.1] $
                      , yrange = [0, max(wbintensity/wbexpt)*1.05]/1000.
      endif else begin
        cgplot, /add, /over, psym = 16, color = colors[ipref] $
                , wbtimes[indx] $
                , wbintensity[indx] / wbexpt[indx] / 1000.
      endelse
      
      ;; Set up for fitting
      if n_elements(fitexpr_in) eq 0 then begin
        case n_elements(indx) of
          1    : fitexpr = ''
          2    : fitexpr = 'P[0] + X*P[1]'
          else : fitexpr = 'P[0] + X*P[1] + X*X*P[2]'
        endcase
      endif else fitexpr = fitexpr_in

      fitexpr_used[ipref] = fitexpr
      
      if fitexpr ne '' then begin

        ;; Do the fit
        err = 1. / wbmu[indx]   ; Large mu is better!
        pp = mpfitexpr(fitexpr $
                       , wbtimes[indx]/3600. $
                       , wbintensity[indx] / wbexpt[indx] $
                       , err $
                      )

        ;; Plot it
        tt = (findgen(1000)/1000 * (max(wbtimes[indx])-min(wbtimes[indx])) + min(wbtimes[indx])) 
        wbint = red_evalexpr(fitexpr, tt/3600, pp)
        cgplot, /add, /over, color = colors[ipref], tt, wbint / 1000.

        ;; Save it
        red_mkhdr, phdr, pp
        anchor = 'DATE'
        red_fitsaddkeyword, anchor = anchor, phdr $
                            , 'FITEXPR', fitexpr, 'mpfitexpr fitting function'
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
      free_lun, lun

      
    endfor                      ; ipref

    ;; Finish the plot
    cglegend, /add $
              , titles = upref + ' : ' + fitexpr_used $
              ;;, titles = upref + ' : ' + coeffs_str $
              , location = [0.9, 0.12], align = 2 $
              , colors = colors, psym = 16, length = 0, vspace = 2
    cgcontrol, output = wbdir+'wb_intensities.pdf'

  endif

  
end


;a = chromisred(/dev)
a = crispred(/dev)

a -> fit_wb_diskcenter
a -> fit_wb_diskcenter, pref = '4846'

pp=readfits('prefilter_fits/wb/wb_fit_6563.fits',h)

print, h, format = '(a0)'

END
