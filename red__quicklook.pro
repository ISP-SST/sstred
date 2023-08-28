; docformat = 'rst'

;+
; Make quick-look movies and images from raw data.
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;     Mats Löfdahl, ISP
;     (Using code from an earlier version by Jaime.)
; 
; :Keywords:
;
;    align : in, optional, type=boolean
;
;      Set this to align the cube. Implies derotate.
;
;    bit_rate : in, optional, type=integer, default=40000 
;
;      Target bit_rate used when calling write_video.
;
;    cam : in, optional, type=string, default="A narrowband camera"
;   
;      Make quicklook for this camera only. 
;   
;    clip : in, optional, type="integer or intarr(4)", default="[50,50,10,10]"
;   
;      A margin (in pixels) to apply to the FOV edges . The order is
;      [left, right, top, bottom]. If a scalar, clip this many pixels
;      from all edges.
;
;    core_and_wings : in, optional, type=boolean
;
;      Automatically select states for which to make quicklook movies.
;      The selection includes the core and extreme wing tunings of
;      each prefilter. If this keyword is set, use_states is ignored.
;      Cannot be combined with /filter_change.
;
;    cube : out, optional, type=fltarr
;
;      If present, return without writing movies, after constructing
;      the first cube. (So use this together with other keywords to
;      specify a single movie.)
;
;    datasets : in, optional, type=strarr
;
;      Timestamp strings that identify datasets to process. Selection
;      menu for data sets bypassed if given.
;   
;    derotate : in, optional, type=boolean
;
;      Set this to derotate the cube to compensate for the field
;      rotation of the telescope.
;
;    destretch : in, optional, type=boolean
;
;      Set this to destretch the cube to compensate for geometrical
;      effects from anisoplanatism. Implies derotate and align.
;
;    filter_change : in, optional, type=boolean
;
;      Automatically select states that are directly after a change of
;      prefilter. Then output only the highest contrast frames but no
;      videos. If this keyword is set, use_states is ignored. Cannot
;      be combined with /core_and_wings.
;
;    format : in, optional, string, default='mp4'
;   
;      The container format of the movies. If 'mp4' or 'avi',
;      write_video will make such a file. If 'mov' (Mac-friendly), a
;      spawned ffmpeg command will convert write_video's mp4 output to
;      the desired format.
;   
;    maxshift : in, optional, type=integer, default=6
;
;      When aligning, filter out isolated shifts larger than this.
;
;    mtf_deconvolve  : in, optional, type=boolean
;
;      Deconvolve data with respect to the diffraction limited MTF of
;      the telescope.
;
;    neuralnet : in, optional, type=boolean
;
;      Do neural net MFBD deconvolution. (So far only in Stockholm.)
;
;    nexp : in, optional, type=integer, default=1
;
;      Number of exposures to base neural net deconvolutions on. 
;
;    no_calib : in, optional, type=boolean
;
;      Do not use dark and flat calibration.
;
;    no_descatter : in, optional, type=boolean
;
;      Do not perform descattering for 7777 and 8542 data.
;
;    no_normalize : in, optional, type=boolean
;
;      Do not normalize intensty levels to median.
;
;    no_plot_r0 : in, optional, type=boolean
;
;      Do not plot per-scan r0 statistics for the datasets.   
;
;    nthreads : in, optional, type=integer, default=6
;
;      Number of threads used in a couple of steps.
;
;    only_plot_r0 : in, optional, type=boolean
;
;      Only plot per-scan r0 statistics for the datasets, no video.
;
;    overwrite :  in, optional, type=boolean
;
;      Overwrite existing movies.
;   
;    scannos : in, optional, type=string
;
;      Choose scan numbers to include in the movies with an array of
;      scan numbers or a comma-and-dash delimited string, like
;      '2-5,7-20,22-30'.
;
;    textcolor : in, optional, type=string, default='yellow'
;   
;      The color to use for text annotations in the movie.
;   
;    use_states : in, optional, type=strarr
;
;      Skip state selection menu, instead make quicklook for states
;      that match any state in this list.    
;   
;    verbose : in, optional, type=boolean
;   
;      Print more screen output.
;   
;    video_codec : in, optional, type=string
;   
;      The codec to use when making the video. Used when calling
;      write_video.
;   
;    video_fps : in, optional, type=integer, default=8
;   
;      Frames per second in the movie. Used when calling write_video.
;   
;    x_flip  : in, optional, type=boolean 
;   
;      Flip the images in the X direction.
;   
;    y_flip  : in, optional, type=boolean 
;   
;      Flip the images in the Y direction.
;   
; 
; :History:
; 
;   2018-07-10 : MGL. First version.
; 
;   2018-07-26 : MGL. Various improvements.
; 
;   2018-08-22 : MGL. Redefine the clip keyword.
; 
;   2018-09-03 : MGL. New keywords no_calib. 
; 
;   2018-09-11 : MGL. New keywords no_descatter, video_codec, format.
;                Change keyword fps to video_fps, bitrate to bit_rate.
; 
;   2019-04-09 : MGL. New keyword core_and_wings.
; 
;   2019-05-22 : MGL. Reorganized the code to reduce reading of file
;                headers. Implemented neuralnet and telescope-MTF
;                deconvolutions. Implemented derotation and
;                destretching. Bugfixes in core_and_wings state
;                selection.
;
;   2022-03-22 : OA. Reorganized the code for selecting spectral
;                points. Changed core_and_wings part to select points
;                in wings instead of first and last ones.
;
;   2022-03-25 : OA. Changed derotation code to avoid FOV clipping 
;                (imported from other part of the pipeline).
;
;   2023-05-30 : MGL. New keyword filter_change.
;
;   2023-08-28 : MGL. Adaptations for mosaic data.
; 
;-
pro red::quicklook, align = align $
                    , bit_rate = bit_rate $
                    , cam = cam $
                    , clip = clip $
                    , core_and_wings = core_and_wings $
                    , cube = cube $
                    , dark = dark $
                    , datasets = datasets $
                    , derotate = derotate $
                    , destretch = destretch $
                    , filter_change = filter_change $
                    , format = format $
                    , gain =  gain $
                    , maxshift = maxshift $
                    , mtf_deconvolve = mtf_deconvolve $
                    , neuralnet = neuralnet $
                    , nexp = nexp $
                    , no_calib = no_calib $
                    , no_descatter = no_descatter $
                    , no_histo_opt = no_histo_opt $
                    , no_plot_r0 = no_plot_r0 $
                    , nthreads = nthreads $
                    , only_plot_r0 = only_plot_r0 $
                    , overwrite = overwrite $
                    , remote_dir = remote_dir $
                    , remote_login = remote_login $
                    , scannos = scannos $
                    , ssh_find = ssh_find $
                    , textcolor = textcolor $
                    , use_states = use_states $
                    , verbose = verbose $
                    , video_codec = video_codec $
                    , video_fps = video_fps $
                    , x_flip = x_flip $
                    , y_flip = y_flip $
                    , min_nscan = min_nscan $
                    , cube_save = cube_save 
  
  inam = red_subprogram(/low, calling = inam1)

  if n_elements(min_nscan) eq 0 then min_nscan=5
  if keyword_set(no_plot_r0) then no_display = 1B else no_display = 0B

  ;; We don't want to use our database for quicklooks on La Palma.
  self.getproperty, 'db_present', db_present
  if db_present then begin
    red_currentsite, site=site
    if site eq 'La Palma' then self.setproperty,'db_present', 0B
  endif

  if n_elements(format) eq 0 then format = 'mp4'
  case format of
    'avi' : extension = format
    'mp4' : extension = format
    else  : extension = 'mp4'
  end

  if keyword_set(mtf_deconvolve) and keyword_set(neuralnet) then begin
    print, inam + ' : Please use only one of the /neuralnet and /mtf_deconvolve keywords.'
    return
  endif
  
  if n_elements(Nexp) eq 0 or ~keyword_set(neuralnet) then Nexp = 1
  if n_elements(scannos) eq 1 && size(scannos, /tname) eq 'STRING' then scannos = rdx_str2ints(scannos)
  if keyword_set(core_and_wings) or keyword_set(filter_change) then undefine, use_states

  ;; The r0 log file is not available until the day after today 
;  if self.isodate eq (strsplit(red_timestamp(/utc,/iso),'T',/extract))[0] then no_plot_r0 = 1
  
  if ~keyword_set(cam) then begin
    cindx = where(~strmatch(*self.cameras,'*-[WD]'), Ncam)
    if Ncam eq 0 then begin
      print, inam + ' : No NB camera in self.cameras?'
      print, *self.cameras
      retall
    endif
    ;; If more than one, pick the first
    cam = (*self.cameras)[cindx[0]]
  endif
  detector = self -> getdetector(cam)
  
  if ~ptr_valid(self.data_dirs) then begin
    print, inam+' : ERROR : undefined data_dir'
    return
  endif

  if n_elements(textcolor) eq 0 then textcolor = 'yellow'
  if n_elements(maxshift) eq 0 then maxshift = 6
  if n_elements(nthreads) eq 0 then nthreads = 6
  if n_elements(bit_rate) eq 0 then bit_rate = 40000
  if n_elements(video_fps) eq 0 then video_fps = 8

  if n_elements(datasets) eq 0 then begin
    ;; Select data sets in a menu.
    tmp = red_select_subset(file_basename(*self.data_dirs) $
                            , qstring = inam + ' : Select data sets' $
                            , count = Nsets, indx = ichoice)
    dirs = (*self.data_dirs)[ichoice] 
  endif else begin
    if max(datasets eq '*') eq 1 then begin
      ;; If called with datasets='*' (or, if an array, any element is
      ;; '*') then process all datasets.
      dirs = *self.data_dirs
    endif else begin
      ;; Select datasets that match the ones listed in the datasets
      ;; keyword.
      match2, datasets, file_basename(*self.data_dirs), suba, subb 
      mindx = where(subb ne -1, Nwhere)
      if Nwhere eq 0 then begin
        print, inam + ' : No data directories match the datasets keyword.'
        return
      endif
      dirs = (*self.data_dirs)[mindx]
    endelse
    Nsets = n_elements(dirs)
  endelse
  
  ;; Now loop over datasets (timestamps)
  for iset = 0, Nsets-1 do begin

    timestamp = file_basename(dirs[iset])

    print, inam + ' : Working on '+timestamp    
    
    ;; Try to limit the number of files we need to extract states for.

    if strmid(cam, 0, 5) eq 'Crisp' then begin
      instrument = 'CRISP'
    endif else begin
      instrument = 'CHROMIS'
    endelse

    ;; Search file names for scan 0, use them to find out what states
    ;; are available.
    ;;files0 = red_file_search('*[_.]00000[_.]*', dirs[iset] + '/' + cam + '/', count = Nfiles)
    files0 = self -> raw_search(dirs[iset] + '/' + cam, scanno = 0, count = Nfiles)

    self -> extractstates, files0, states0
        
    upref = states0(uniq(states0.prefilter, sort(states0.prefilter))).prefilter
    Npref = n_elements(upref)
    indx = uniq(states0.tun_wavelength, sort(states0.tun_wavelength))
    ustat = states0[indx].fullstate        

    states_count = 0
    undefine, pat, ustat2, sel_in
    case 1 of

      keyword_set(filter_change) : begin

        tunings = self -> tunings_after_filterchange(dirs[iset] + '/' + cam, count = Ntunings)
        
        for ituning = 0, Ntunings-1 do begin
          ;; The _ added to ustat needed to match chromis data. Didn't
          ;; work to add the end-of-line character ($) to the regex.
          sindx = where(strmatch(ustat+'_', '*_'+tunings[ituning]+'_*'), Nmatch)
          red_append, ustat2, ustat[sindx[0]]
          red_append, sel_in, sindx[0]
        endfor                  ; ituning
        
        states_count = n_elements(ustat2)
        for istate = 0, states_count-1 do begin
          fn = states0[indx[sel_in[istate]]].filename
          prts = strsplit(fn,'[_.]',/extract)
          if strmatch(cam,'*W*') then begin
            if instrument eq 'CHROMIS' then $
              red_append, pat, '*' + prts[-3] + '*' $
            else $
              red_append, pat, '*' +prts[-5] + '*'
          endif else begin
            if instrument eq 'CHROMIS' then $
              red_append, pat, '*' + prts[-3] + '*' + prts[-2] + '*' $
            else $
              red_append, pat, '*' +prts[-5] + '*' + prts[-4] + '*' + prts[-3] + '*'
          endelse
        endfor                 ; istate
        use_states = ustat2
        
      end                       ; filter_change
      
      keyword_set(core_and_wings) : begin
        
        for ipref = 0, Npref-1 do begin
          sindx = where(strmatch(ustat, '*_'+upref[ipref]+'_*'), Nmatch)
          if Nmatch eq 1 then begin
            ;; If just one state for this prefilter, then use it!
            red_append, ustat2, ustat[sindx[0]]
            red_append, sel_in, 0
          endif else begin              
            imatch = where(strmatch(ustat[sindx], '*+0*'), Nmatch)
            if Nmatch gt 1 then begin
              print, 'There is more then one spectral line in a scan with ', upref[ipref], ' prefilter.'
              print, 'You have to choose spectral points manually.'
              undefine, ustat2
              break
            endif
            if Nmatch gt 0 then begin
              red_append, ustat2, ustat[sindx[imatch]]
              red_append, sel_in, sindx[imatch]
              wv = states0[indx[sindx]].tun_wavelength
              if n_elements(wv[0:imatch]) mod 2 ne 1 then ex=1 else ex=0
              wng = median(wv[0:imatch-ex])
              in = where(wv eq wng)
              red_append, sel_in, sindx[in]
              red_append, ustat2, ustat[sindx[in]]
              if n_elements(wv[imatch:*]) mod 2 ne 1 then ex=1 else ex=0
              wng = median(wv[imatch+ex:*])
              in = where(wv eq wng)
              red_append, sel_in, sindx[in]
              red_append, ustat2, ustat[sindx[in]]
           endif else begin
              np = n_elements(sindx)
              red_append, ustat2, ustat[sindx[np/4]]
              red_append, ustat2, ustat[sindx[np/2]]
              red_append, ustat2, ustat[sindx[np/2+np/4]]
              red_append, sel_in, [sindx[np/4],sindx[np/2],sindx[np/2+np/4]]
            endelse
          endelse
        endfor                 ; ipref

        states_count = n_elements(ustat2)
        for istate = 0, states_count-1 do begin
          fn = states0[indx[sel_in[istate]]].filename
          prts = strsplit(fn,'[_.]',/extract)
          if strmatch(cam,'*W*') then begin
            if instrument eq 'CHROMIS' then $
              red_append, pat, '*' + prts[-3] + '*' $
            else $
              red_append, pat, '*' +prts[-5] + '*'
          endif else begin
            if instrument eq 'CHROMIS' then $
              red_append, pat, '*' + prts[-3] + '*' + prts[-2] + '*' $
            else $
              red_append, pat, '*' +prts[-5] + '*' + prts[-4] + '*' + prts[-3] + '*'
          endelse
        endfor                 ; istate
        use_states = ustat2

      end                       ; core_and_wings

      n_elements(use_states) gt 0 : begin
        
        for istate = 0,n_elements(use_states)-1 do begin
          imatch = where(strmatch(ustat, '*'+use_states[istate]+'*'), Nmatch)
          if Nmatch ge 1 then begin
            states_count++
            red_append, ustat2, ustat[imatch] 
            fn = states0[indx[imatch]].filename
            prts = strsplit(fn,'[_.]',/extract)
            if strmatch(cam,'*W*') then begin
              if instrument eq 'CHROMIS' then $
                 red_append, pat, '*' + prts[-3] + '*' $
              else $
                 red_append, pat, '*' +prts[-5] + '*'
            endif else begin
              if instrument eq 'CHROMIS' then $
                 red_append, pat, '*' + prts[-3] + '*' + prts[-2] + '*' $
              else $
                 red_append, pat, '*' +prts[-5] + '*' + prts[-4] + '*' + prts[-3] + '*'
            endelse
          endif else print, 'There is no match for ', use_states[istate], ' state.'
        endfor                  ; istate
        if states_count eq 0 then print,'There are no matches for provided states. You have to choose states manually.'
        
      end                       ; use_states
      
      else : 
      
    endcase

    
    if states_count eq 0 then begin 
      tmp = red_select_subset(states0[indx].prefilter+'_'+states0[indx].fullstate $
                              , qstring = inam + ' : Select states' $
                              , count = Nstates, indx = ichoice)
      ustat2 = ustat[ichoice]
      for istate = 0,n_elements(ustat2)-1 do begin
        imatch = where(ustat2[istate] eq ustat)
        fn = states0[indx[imatch]].filename
        prts = strsplit(fn,'[_.]',/extract)
        if strmatch(cam,'*W*') then begin
          if instrument eq 'CHROMIS' then $
            red_append, pat, '*' + prts[-3] + '*' $
          else $
            red_append, pat, '*' +prts[-5] + '*'
        endif else begin
          if instrument eq 'CHROMIS' then $
            red_append, pat, '*' + prts[-3] + '*' + prts[-2] + '*' $
          else $
            red_append, pat, '*' +prts[-5] + '*' + prts[-4] + '*' + prts[-3] + '*'
        endelse
      endfor                 ; istate
      use_states = ustat2      
    endif

    if n_elements(scannos) gt 0 then begin
      sc = '*' + string(scannos, format='(I05)')
      for ii=0, Nstates-1 do $
        red_append, patt, sc + pat[ii]
      pat = patt
    endif
    files = red_file_search(pat, dirs[iset] + '/' + cam, count = Nfiles)

    if instrument eq 'CRISP' then begin
      ;; Check for lcd files!
      windx = where(~strmatch(files, '*.lcd.*'), Nwhere)
      if Nwhere eq 0 then files = '' else begin
        files = files[windx]
      endelse
      Nfiles = Nwhere      
    endif
    
    if files[0] eq '' then begin
      print, inam + ' : ERROR -> no frames found in '+dirs[iset]
      continue                  ; Goto next dataset
    endif
    
    self -> extractstates, files, states

    ;; If this is a mosaic directory we should 1) not select best
    ;; scan, because there is only one, 2) make images of all tiles,
    ;; 3) not make video, 4) make sure we get a sensibe r0 plot.
    ismos = strmatch(files[0], '*mos00*')
       
    nsc = max(states.scannumber)
    if ~keyword_set(filter_change) && nsc lt min_nscan && ~ismos then begin
      ;; Ignore short datasets (likely disk center intensity
      ;; calibrations) unless the filter_change keyword is set.
      ;; But allow short data sets in the "filter_change" mode and for
      ;; mosaic data.
      fn = states[0].filename
      date = stregex(fn, '20[0-2][0-9][.-][01][0-9][.-][0-3][0-9]', /extract)
      timestamp = stregex(fn, '[0-2][0-9]:[0-5][0-9]:[0-6][0-9]', /extract)
      print, 'Dataset ', date + ' ' + timestamp + ' is too short.' 
      print,"We do not want to bother with short datasets. Skipping it."
      continue
    endif

    outdir = self.out_dir +'/quicklook/'+timestamp+'/'
    file_mkdir, outdir

    if ~keyword_set(no_plot_r0) then begin

      ;; Plot r0

      cgcontrol, /delete, /all

      pname = states[0].camera
      pname += '_' + 'r0'
      pname += '_' + self.isodate
      pname += '_' + timestamp
      pname += '.pdf'

      if keyword_set(overwrite) || ~file_test(outdir+pname) then begin
        red_plot_r0_stats, states, pname = outdir+pname, ismos = ismos
      endif
      
    endif

    
    if keyword_set(only_plot_r0) then continue     

    ustat = ustat2
    Nstates = n_elements(ustat)
        
    for istate = 0, Nstates-1 do begin
      
      self -> selectfiles, files = files, states = states $
                           , ustat = ustat[istate] $
                           , selected = sel, count = Nsel

      if Nsel eq 0 then continue
      
      uscan = states[sel[uniq(states[sel].scannumber,sort(states[sel].scannumber))]].scannumber

      Nscans = n_elements(uscan)
      Nexp_available = Nsel/Nscans
      
      sel = sel[sort(states[sel].scannumber)] ; Make sure scans are in order!

      if keyword_set(no_calib) then begin
        dd = 0.
        gg = 1.
      endif else begin
        ;; Load dark and flat
        self -> get_calib, states[sel[0]] $
                           , darkdata = dd, darkstatus = darkstatus $
                           , gaindata = gg, gainstatus = gainstatus
        if darkstatus ne 0 then dd = 0.
        if gainstatus ne 0 then gg = 1.
      endelse
      
      print, inam + ' : found scans -> '+red_stri(Nscans)

      pref = states[sel[0]].prefilter

      DoBackscatter = 0
      if ~keyword_set(no_descatter) $
         && self.dodescatter $
         && (pref eq '8542' || pref eq '7772') then begin
        self -> loadbackscatter, detector, pref, bgt, Psft
        DoBackscatter = 1
      endif

      namout = states[sel[0]].camera
      namout += '_' + 'quick'
      namout += '_' + self.isodate
      namout += '_' + timestamp
      namout += '_' + pref $
                + '_' + states[sel[0]].tuning

      if ismos then begin

        ;; Mosaic directory: no video, save images for all tiles.
        Nmos = Nsel
        umos = reform((stregex(files[sel], 'mos([0-9][0-9])', /extract, /subexpr))[1, *])
        
        for imos = 0, Nmos-1 do begin

          hdr = red_readhead(files[sel[imos]])
          red_fitspar_getdates, hdr, date_avg = date_avg
          date_avg = red_strreplace(date_avg, 'T', ' ')
          date_avg = (strsplit(date_avg,'.',/extr))[0] ; No decimal seconds
          
          ims = red_readdata(files[sel[imos]])
          dim = size(ims, /dim)
          if n_elements(dim) eq 2 then Nframes = 1 else Nframes = dim[2]

          if (n_elements(gg) ne 1) && (gainstatus eq 0) then $
             mask = red_cleanmask(gg ne 0)

          ;; Do dark and gain  correction.
          for iframe = 0, Nframes-1 do begin
            ims[*, *, iframe] -= dd
            ims[*, *, iframe] *= gg

            if  strmatch(cam,'*W*') then $
               im = red_fillpix(ims[*, *, iframe], nthreads=nthreads) $
            else $
               im = red_fillpix(ims[*, *, iframe], mask=mask, nthreads=nthreads)

            idx1 = where(im eq 0.0, Nwhere, complement = idx, Ncomplement = Nnowhere)
            if Nwhere gt 0 && Nnowhere gt 0 then im[idx1] = median(im[idx])

;            if(keyword_set(x_flip)) then im = reverse(temporary(im), 1)
;            if(keyword_set(y_flip)) then im = reverse(temporary(im), 2)
;
            ims[*, *, iframe] = im

          endfor                ; iframe

          x0 = 0
          x1 = dim[0] - 1
          y0 = 0
          y1 = dim[1] - 1
          case n_elements(clip) of
            0 : begin
              x0 += 50
              x1 -= 50
              y0 += 10
              y1 -= 10
            end
            1 : begin
              x0 += clip
              x1 -= clip
              y0 += clip
              y1 -= clip
            end
            4 : begin
              x0 += clip[0]
              x1 -= clip[1]
              y0 += clip[2]
              y1 -= clip[3]
            end
            else : stop
          endcase

          im = ims[*, *, 0]     ; Select best image?
          ;; if ~keyword_set(no_histo_opt) then im = red_histo_opt(im)
          im = bytscl(im)

          ;; Add annotations and make the rgb cube
          thisDevice = !D.Name
          set_plot, 'Z'
          device, Set_Resolution = [x1-x0+1, y1-y0+1] $
                  , Z_Buffer = 0 $
                  , Set_Pixel_Depth = 24 $
                  , Decomposed=0 $
                  , set_font="Helvetica" $
                  , /tt_font
          !P.font = 1


          erase
          tv, im[x0:x1, y0:y1]



          
          ;; Annotation
          annstring = date_avg $ ;self.isodate $ ;+ ' ' + time_avg[iscan] $
                      + '   ' + ustat[istate] $
                      + '   tile : ' + umos[imos] 
          cgtext, [0.01], [0.95], annstring $
                  , /normal, charsize=3., color=textcolor, font=1
          
          rgbcube = tvrd(/true)
          set_plot,'X'
          
          ;; Make jpeg images.
          jname = outdir+namout+'_mos'+umos[imos] + '.jpg'

          file_delete, jname, /allow
          write_jpeg, jname, rgbcube, q = 100, /true

          print, jname


        endfor                  ; imos
        

      endif else begin
        
        print, inam + ' : Cube '+namout
        
        if keyword_set(neuralnet) then namout += '_NN' 
        if keyword_set(mtf_deconvolve) then namout += '_MTF'
        if keyword_set(cube_save) then begin
          cubnam = namout +  '.fits'
          contrastnam = namout + '_contrast.fits'
        endif
        namout += '.' + extension

        if ~keyword_set(overwrite) && file_test(outdir+namout) then continue

        print, inam + ' : saving to folder -> '+outdir 
        
        dim = (fxpar(red_readhead(files[sel[0]]),'NAXIS*'))[0:1]
        x0 = 0
        x1 = dim[0] - 1
        y0 = 0
        y1 = dim[1] - 1
        case n_elements(clip) of
          0 : begin
            x0 += 50
            x1 -= 50
            y0 += 10
            y1 -= 10
          end
          1 : begin
            x0 += clip
            x1 -= clip
            y0 += clip
            y1 -= clip
          end
          4 : begin
            x0 += clip[0]
            x1 -= clip[1]
            y0 += clip[2]
            y1 -= clip[3]
          end
          else : stop
        endcase

        ;; RGB cube
        cube = fltarr(dim[0], dim[1], Nexp <Nexp_available, Nscans)
        best_contrasts = fltarr(Nscans)
        time_avg = strarr(Nscans)
        
        if (n_elements(gg) ne 1) && (gainstatus eq 0) then $
           mask = red_cleanmask(gg ne 0)

        for iscan = 0L, Nscans -1 do begin

          red_progressbar, iscan, Nscans, 'Read and select files', /predict

          self -> selectfiles, files = files, states = states $
                               , ustat = ustat[istate] $
                               , scan = uscan[iscan] $
                               , selected = sel2

          docontinue = 0
          case n_elements(sel2) of
            0 : begin
              if iscan ne Nscans-1 then stop
              ;; The last scan does not include this state. Make the
              ;; movie one frame shorter!
              Nscans--
              cube = cube[*, *, 0:Nscans-1]
              best_contrasts = best_contrasts[0:Nscans-1]
              time_avg = time_avg[0:Nscans-1]
              docontinue = 1
            end
            1 : begin
              ;; Just a single file
              ims = float(red_readdata(states[sel2[0]].filename))
              Nframes = (fxpar(red_readhead(files[sel2[0]]),'NAXIS3') >1)
              red_fitspar_getdates, red_readhead(states[sel2[0]].filename), date_avg = date_avg
            end
            else : begin
              ;; Multiple files, e.g. CHROMIS WB. Or any CRISP camera.
              ;; (Ideally, you should be able to call red_readdata with
              ;; multiple file names and if would do the right thing.
              ;; Needs to produce a combined header, too!)
              Nframes = 0L
              red_fitspar_getdates, red_readhead(states[sel2[0]].filename), date_avg = date_avg
              for ifile = 0, n_elements(sel2)-1 do begin
                Nframes += (fxpar(red_readhead(files[sel2[ifile]]),'NAXIS3') >1)
              endfor            ; ifile
              ims = fltarr(dim[0], dim[1], Nframes)
              iframe = 0
              for ifile = 0, n_elements(sel2)-1 do begin
;              print, 'Read '+file_basename(files[sel2[ifile]])
                ims[0, 0, iframe] = red_readdata(files[sel2[ifile]], head = head)
                iframe += (fxpar(head,'NAXIS3') >1)
              endfor            ; ifile
            end
          endcase
          if docontinue then continue ; Not allowed in the case statement

          time_avg[iscan] = (strsplit((strsplit(date_avg, 'T', /extract))[1], '.', /extract))[0]

          contrasts = fltarr(Nframes)

          
          ;; Do dark and gain (and possibly backscatter) correction.
          for iframe = 0, Nframes-1 do begin
            ims[*, *, iframe] -= dd
            if DoBackscatter then begin ; Backscatter if needed
              ims[*, *, iframe] = rdx_descatter(ims[*, *, iframe] $
                                                , bgt, Psft $
                                                , verbose = verbose $
                                                , nthreads = nthreads)
            endif
            ims[*, *, iframe] *= gg

            if  strmatch(cam,'*W*') then $
               im = red_fillpix(ims[*, *, iframe], nthreads=nthreads) $
            else $
               im = red_fillpix(ims[*, *, iframe], mask=mask, nthreads=nthreads)

            idx1 = where(im eq 0.0, Nwhere, complement = idx, Ncomplement = Nnowhere)
            if Nwhere gt 0 && Nnowhere gt 0 then im[idx1] = median(im[idx])

            if(keyword_set(x_flip)) then im = reverse(temporary(im), 1)
            if(keyword_set(y_flip)) then im = reverse(temporary(im), 2)

            ims[*, *, iframe] = im

          endfor                ; iframe
          
          for iframe = 0, Nframes-1 do contrasts[iframe] $
             = stddev(ims[x0:x1, y0:y1, iframe])/mean(ims[x0:x1, y0:y1, iframe])

          if keyword_set(neuralnet) then begin

            ;; Deconvolve the data with a neural net.

            ;; Select the highest contrast frames
            nn_indx = reverse(sort(contrasts))
            cube[0, 0, 0, iscan] = ims[*, *, nn_indx[0:(Nexp <Nexp_available)-1]]
            
          endif else begin

            if size(ims, /n_dim) eq 3 then begin
              ;; Select the best frame
              ;; Possible improvement: base this selection on wb contrast?
              mx = max(contrasts, ml)
              im = ims[*, *, ml] 
              best_contrasts[iscan] = contrasts[ml]
            endif else begin
              im = ims 
              best_contrasts[iscan] = contrasts[0]
            endelse
            
          endelse
          
          cube[*, *, 0, iscan] = im

        endfor                  ; iscan

        if keyword_set(neuralnet) then begin

          if ~file_test('~reduc/mfbdNN/encdec_sst.py') then begin
            print, inam + ' : You do not seem to have the neural net software installed.'
            return          
          endif
          
          ;; Run the neural net
          ts = cgTimeStamp(11, RANDOM_DIGITS=12, /UTC)
          tmpoutfile = '/tmp/quick_raw'+ts+'.fits'
          tmpinfile = '/tmp/quick_deconvolved'+ts+'.fits'

          ;; Make the command to process the data with the NN.
          nn_cmd = 'cd ~reduc/mfbdNN/ ; python encdec_sst.py -i '+tmpoutfile+' -o '+tmpinfile 
          
          ;; Write the data to disk
          red_mkhdr, hhh, cube
          anchor = 'DATE'
          red_fitsaddkeyword, anchor = anchor, hhh $
                              , 'WAVELNTH' $
                              , float((strsplit(ustat[istate],'_',/extract))[0]) / 10. $
                              , 'Wavelength based on prefilter'
          red_fitsaddkeyword, anchor = anchor, hhh, 'WAVEUNIT', -9 $
                              , 'Wavelength unit 10^WAVEUNIT m = nm'
          red_fitsaddkeyword, anchor = anchor, hhh, 'INSTRUME', instrument $
                              , 'CRISP or CHROMIS data?'
          writefits, tmpoutfile, cube, hhh


          
          ;; Spawn the NN command and wait for it to terminate
          print, 'Do neural net deconvolution of '+strtrim(Nscans, 2)+' scans:'
          tic
          spawn, nn_cmd
          toc
          
          ;; Read data back into the same cube variable, which will now
          ;; have one dimension less.
          cube = readfits(tmpinfile, /silent)

          file_delete, tmpinfile, tmpoutfile

          ;; Caclulate the contrasts of the deconvolved images
          for iscan = 0, Nscans-1 do begin
            best_contrasts[iscan] = stddev(cube[x0:x1, y0:y1, iscan])/mean(cube[x0:x1, y0:y1, iscan])
          endfor
          
        endif else begin

          ;; Remove the extra dimension
          cube = reform(cube, dim[0], dim[1], Nscans)
          
        endelse
        
        if keyword_set(mtf_deconvolve) then begin

          ;; This part needs to 1) not use non-pipeline subprograms and
          ;; 2) be modified to support aspect ratios != 1.
          
          ;; caminfo = red_camerainfo(detector)
          
          lambda = states[sel2[0]].tun_wavelength ; Wavelength [m]
          telescope_d = 0.97d
          arcsecperpix = self.image_scale
          ;; pixelsize = caminfo.pixelsize
          sz = max(dim)
          
          ;;   F_number = pixelsize/telescope_d/(arcsecperpix*2d*!dpi/(360.*3600.))
          ;;   Q_number = F_number * lambda/pixelsize

          Q_number = lambda/telescope_d/(arcsecperpix*2d*!dpi/(360.*3600.)) 
          
          LimFreq = sz / Q_number
          rc = LimFreq/2.d
          r = round(rc)+2
          r = sz/4

          x_coord = (findgen(sz, sz) mod sz)-sz/2
          y_coord = transpose(x_coord)
          r_coord = sqrt(x_coord*x_coord+y_coord*y_coord)
          pupil = r_coord lt limfreq/2
          ap = r_coord lt limfreq
          
;        pupil = double(aperture(2*r, rc))
;        ap = double(aperture(2*r, 2*rc))

          psf = abs(fft(pupil))^2 / total(pupil)
          mtf = double(fft(psf, /inv))
          mtf = mtf/mtf[0, 0]

          ;; Make a window in an array the size of a padded (in case of
          ;; aspect ratio != 1) image.
          w = fltarr(sz, sz)
          w[0:dim[0]-1, 0:dim[1]-1] = hanning(dim[0], dim[1])
          dimdiff = dim[0] - dim[1] 

          ;; Deconvolve and caclulate the contrasts of the deconvolved images
          for iscan = 0, Nscans-1 do begin
            
            red_progressbar, iscan, Nscans, 'MTF deconvolution.'
            im = double(cube[*, *, iscan])
            md = median(im)
            im = im - md

            if dimdiff ne 0 then begin
              im2 = fltarr(sz, sz)
              im2[0:dim[0]-1, 0:dim[1]-1] = im
              im = im2
            endif
            
            fim = fft(w*im)     ; FFT of windowed image
            
            ;; Make a low-pass noise filter
            nl = red_noiselevel(fim, limfreq = limfreq, /Fourier) / sz
            filt = smooth(abs(shift(fim, sz/2, sz/2))^2,15) / (smooth(abs(shift(fim, sz/2, sz/2))^2,15) + nl^2*4)
            nl_filt = median(filt(where(~ap))) ; Level outside diffraction limit
            mask = filt gt nl_filt*2
            ;; Clean the filter from high-frequency noise contributions
            x = INDGEN(16*16) MOD 16 + sz/2
            y = INDGEN(16*16) / 16 + sz/2
            roiPixels = x + y * sz
            mindx = region_grow(mask,roipixels,threshold=[0.9,1.1])
            mask[mindx] = 255
            mask = mask gt 200
            filt = shift(smooth(filt*mask, 15), sz/2, sz/2)
            
            im = float(fft(filt*(fim/(mtf >1e-3)), /inv)) ; Deconvolved image
            im /= (w >3e-4)                               ; Undo the windowing
            
            cube[0, 0, iscan] = im[0:dim[0]-1, 0:dim[1]-1] + md ; Add the median back
            best_contrasts[iscan] = stddev(cube[x0:x1, y0:y1, iscan])/mean(cube[x0:x1, y0:y1, iscan])
            
          endfor
          
        endif
        
        if Nscans gt 3 && ~keyword_set(no_normalize) then begin
          ;; Normalize intensities
;        med = median(median(cube, dim = 1), dim = 1)
          med = fltarr(Nscans)
          for iscan = 0, Nscans-1 do med[iscan] = median(cube[*, *, iscan])
          mm = mean(med)
          cc = poly_fit(findgen(Nscans), med/mm, 2, yfit=yfit)
          for iscan = 0, Nscans-1 do cube[*,*,iscan] /= yfit[iscan]*mm
        endif
        
        if keyword_set(derotate) || keyword_set(align) || keyword_set(destretch) then begin
          ;; Derotate images
          ang = red_lp_angles(time_avg, self.isodate)
          mang = median(ang)
          ang -= mang
          case cam of
            'Crisp-W' : ang = -ang
;          'Crisp-T' : ang = -ang
            else :
          endcase
          for iscan = 0L, Nscans -1 do begin
            red_progressbar, iscan, Nscans, 'Derotating images.'
            cube[*,*,iscan] = red_rotation(cube[*,*,iscan], ang[iscan], nthreads = nthreads)
          endfor                ; iscan
        endif

        if keyword_set(align) || keyword_set(destretch) then begin
          ;; Align images
          
          ;; Measure image shifts
          shifts = red_aligncube(cube, 5, /center $ ;, cubic = -0.5 $
                                 , xbd = round(dim[0]*.9) $
                                 , ybd = round(dim[1]*.9) $
                                 , no_display = no_display $
                                 , nthreads = nthreads)

          if Nscans gt 3 then begin
            ;; Outliers?
            indx_included = where((abs(shifts[0,*] - median(reform(shifts[0,*]),3)) le maxshift) $
                                  and (abs(shifts[1,*] - median(reform(shifts[1,*]),3)) le maxshift) $
                                  , complement = indx_excluded, Nincluded, Ncomplement = Nexcluded)
            if Nexcluded gt 0 then begin
              shifts[0, indx_excluded] = interpol(shifts[0, indx_included], indx_included, indx_excluded)
              shifts[1, indx_excluded] = interpol(shifts[1, indx_included], indx_included, indx_excluded)
            endif
          endif
          
          ;; Align the cube
          for iscan = 0, Nscans-1 do begin
            red_progressbar, iscan, Nscans, 'Applying the shifts.'
            cube[*, *, iscan] = red_shift_im(cube[*, *, iscan] $
                                             , shifts[0, iscan] $
                                             , shifts[1, iscan] $
                                             , cubic = -0.5 $
                                             , nthreads = nthreads)
          endfor                ; iscan
        endif

        if keyword_set(destretch) then begin

          ;; Destretch images
          
          if n_elements(clips) eq 0 then clips = [12,  6,  3]
          if n_elements(tiles) eq 0 then tiles = [10, 20, 30]

          dts = red_time2double(time_avg)
          if n_elements(tstep) eq 0 then begin
            tstep = fix(round(180. / median(abs(dts[0:Nscans-2] - dts[1:*]))))
          endif
          tstep = tstep < (Nscans-1)

          ;; Calculate stretch vectors
          grid = red_destretch_tseries(cube, 1.0/float(self.image_scale), tiles, clips, tstep, nthreads = nthreads)

          for iscan = 0L, Nscans - 1 do begin
            red_progressbar, iscan, Nscans, 'Applying the stretches.'
            cube[*,*,iscan] = red_stretch(cube[*,*,iscan], reform(grid[iscan,*,*,*]))
          endfor                ; iscan
        endif

        if arg_present(cube) then return
        
        if ~keyword_set(no_histo_opt) then cube = red_histo_opt(temporary(cube))

        cube = bytscl(cube[x0:x1, y0:y1, *])
        rgbcube = bytarr(3, x1-x0+1, y1-y0+1, Nscans)
        
        ;; Add annotations and make the rgb cube
        thisDevice = !D.Name
        set_plot, 'Z'
        device, Set_Resolution = [x1-x0+1, y1-y0+1] $
                , Z_Buffer = 0 $
                , Set_Pixel_Depth = 24 $
                , Decomposed=0 $
                , set_font="Helvetica" $
                , /tt_font
        !P.font = 1
        ;; k = 0L
        for iscan=0, Nscans-1 do begin

          erase
          tv, cube[*,*,iscan]

          ;; Annotation
          annstring = self.isodate + ' ' + time_avg[iscan] $
                      + '   ' + ustat[istate] $
                      + '   scan :' + string(uscan[iscan],format='(I5)') 
          cgtext, [0.01], [0.95], annstring $
                  , /normal, charsize=3., color=textcolor, font=1

;        print, date_avg
;        print, time_avg[iscan]
;        print, annstring
;        stop
          
          ;; tvrd(/true) reads an RGB image [3,Nx,Ny]
          snap2 = tvrd(/true)

          rgbcube[0, 0, 0, iscan] = snap2

        endfor                  ; iscan

        if keyword_set(cube_save) then begin
          writefits,outdir+cubnam,cube
          writefits, outdir+contrastnam, best_contrasts
        endif
        set_plot,'X'

        if ~keyword_set(filter_change) then begin ; Skip making video if only checking for filter change
          
          ;; We could write metadata to the video file with the
          ;; write_video METADATA keyword:
          ;;
          ;; METADATA
          ;; 
          ;; Set this keyword to a [2, n] element string array denoting n
          ;; [key, value] pairs of metadata to be written to the file. 
          ;; 
          ;; Note: Metadata must be written before any video or audio
          ;; data.
          ;;
          ;; The keys are limited to this set: [album, artist, comment,
          ;; copyright, genre, title]. We could perhaps use comment,
          ;; genre, and title.
          ;;
          ;; Something like:
          ;; genre = "SST quicklook movie"
          ;; title = instrument + date + timestamp + state
          ;; comment = Things like date when made, pipeline version?
          metadata = strarr(2, 3)
          metadata[*, 0] = ['genre', 'SST quicklook movie']
          metadata[*, 1] = ['title', strjoin([cam $
                                              , self.isodate $
                                              , timestamp $
                                              , pref + '_' + states[sel[0]].tuning $
                                             ], ' ')]

          metadata[*, 2] = ['comment', 'Movie made ' $
                            + red_timestamp(/utc, /iso) $
                            + ' UTC' $
                           ]


;      dims = size(rgbcube, /dim)
;      rgbcube = byte(255*randomu(seed, [3, 100, 100, 13]))

;      help, bit_rate, video_fps

          file_delete, outdir+namout, /allow
          write_video, outdir+namout, rgbcube $
                       , bit_rate = bit_rate $
                       , metadata = metadata $
                       , video_fps = video_fps $
                       , video_codec = video_codec 
          
          if format eq 'mov' then begin
            ;; Convert to Mac-friendly (and smaller) .mov file using recipe from Tiago
            mname = outdir + red_strreplace(namout, '.'+extension,'.'+format)
            file_delete, mname, /allow_nonexist
            spawn, 'ffmpeg -n -i "' + outdir + namout $
                   + '" -c:v libx264 -preset slow -crf 26 -tune grain "' $
                   + mname + '"'
            file_delete, outdir + namout
;        spawn, 'rm "' + outdir + namout + '"'
;        find . -name '*mp4' -exec sh -c 'ffmpeg -n -i "$1" -c:v libx264 -preset slow -crf 26 -vf scale=-1:800  -tune grain "${1%.mp4}.mov"' sh {} \ ;
          endif

        endif
        
;      stop
        
        ;; Make a jpeg image of the best frame.
        if Nscans eq 1 then ml = 0 else mx = max(best_contrasts, ml)
        jname = outdir+red_strreplace(namout, '.'+extension, '_scan='+strtrim(uscan[ml], 2)+'.jpg')

        file_delete, jname, /allow
        write_jpeg, jname, rgbcube[*, *, *, ml], q = 100, /true

        print, outdir+namout
        print, jname
        print, strjoin(strtrim(size(cube, /dim), 2), ' x ')

      endelse
      
    endfor                      ; istate
  endfor                        ; iset
  
end

;; Todo: Change again the structure of the program. First find all
;; states from scan 0. Then use that to determine what states to
;; produce movies for, whether decided with use_states or
;; core_and_wings or the selection dialogue. Then loop over those
;; states, and in that loop search for files that match that state and
;; the list of scan numbers (if given). This should make it much
;; quicker to make movies with selected scans, and also to return the
;; cube for such movies as needed by make_raw_cube.
