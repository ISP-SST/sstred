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
;      Set this to align the cube.
;
;    bitrate : inm optional, type=integer, default=40000 
;
;      Target bitrate used when making movies. 
;
;    cam : in, optional, type=string, default="A narrowband camera"
;   
;      Make quicklook for this camera only. 
;   
;    clip : in, optional, type="integer or intarr(4)", default="[50,50,10,10]"
;   
;      A margin (in pixels) to apply to the FOV edges . The order is
;      [left, right, top, bottom]. If not array, clip this many pixels
;      from all edges.
;   
;    datasets : in, optional, type=strarr
;
;      Timestamp strings that identify datasets to process. Selection
;      menu for data sets buypassed if given.
;   
;    fps : in, optional, type=integer, default=8
;   
;      Frames per second in the movie. 
;   
;    maxshift : in, optional, type=integer, default=6
;
;      When aligning, filter out isolated shifts larger than this.
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
;    overwrite :  in, optional, type=boolean
;
;      Overwrite existing movies.
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
; 
;   2018-08-22 : MGL. Redefine the clip keyword.
; 
;-
pro red::quicklook, align = align $
                    , bitrate = bitrate $
                    , cam = cam $
                    , clip = clip $
                    , dark = dark $
                    , datasets = datasets $
                    , destretch = destretch $
                    , fps = fps $
                    , gain =  gain $
                    , maxshift = maxshift $
                    , no_histo_opt = no_histo_opt $
                    , nthreads = nthreads $
                    , overwrite = overwrite $
;                    , pattern = pattern $
                    , no_plot_r0 = no_plot_r0 $
                    , remote_dir = remote_dir $
                    , remote_login = remote_login $
                    , ssh_find = ssh_find $
                    , textcolor = textcolor $
                    , use_states = use_states $
                    , verbose = verbose $
                    , x_flip = x_flip $
                    , y_flip = y_flip 
  
  inam = red_subprogram(/low, calling = inam1)

  if ~keyword_set(cam) then begin
    cindx = where(~strmatch(*self.cameras,'*-[WD]'), Ncam)
    if Ncam eq 0 then begin
      print, inam + ' : No NB camera in self.cameras?'
      print, *self.cameras
      retall
    endif
    ;; If more than one, pick the first
    cam = (*self.cameras)[cindx[0]]
    detector = self -> getdetector(cam)
  endif
  
  if ~ptr_valid(self.data_dirs) then begin
    print, inam+' : ERROR : undefined data_dir'
    return
  endif

  if n_elements(textcolor) eq 0 then textcolor = 'yellow'
  if n_elements(maxshift) eq 0 then maxshift = 6
  if n_elements(nthreads) eq 0 then nthreads = 6
  if n_elements(bitrate) eq 0 then bitrate = 40000
  if n_elements(fps) eq 0 then fps = 8

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
  
  ;; If search pattern is given, use that, otherwise just use *
;  if keyword_set(pattern) then pat = "*" + pattern + "*" else
  
  if strmid(cam, 0, 5) eq 'Crisp' $
     && n_elements(use_states) gt 0 then begin
    ;; If use_states is provided, we base the file search pattern on
    ;; that. This works for CRISP, where the state is part of the file
    ;; name.
    
    ;; We need to do a bit of massaging here, because the CRISP states
    ;; in the file names have the tuning (after the sign) zero-padded
    ;; to 4 digits, while the states returned by extractstates (used
    ;; below) are not zero padded. We want this to work whether the
    ;; use_states are given zero padded or not. The pattern used for
    ;; file seaching needs the padding but any padding has to be gone
    ;; when we get to the state comparison later. To make it even more
    ;; complicated, we don't know if this tuning is even part
    ;; of the use_states! 
    pat = strarr(n_elements(use_states))
    for istate = 0, n_elements(use_states)-1 do begin
      st = stregex(use_states[istate], '(_|\.|^)([+-][0-9]*)(_|\.|$)' $
                   , /extract,/sub)
      if st[2] ne '' then begin
        ;; We had a match for the tuning part of the state
        tun = st[2]
        tun_padded    = strmid(tun, 0, 1) + string(long(strmid(tun, 1)), format='(i04)')
        tun_nonpadded = strmid(tun, 0, 1) + string(long(strmid(tun, 1)), format='(i0)')
        pat[istate] = '*' + red_strreplace(use_states[istate], tun, tun_padded) + '*'
        use_states[istate] = red_strreplace(use_states[istate], tun, tun_nonpadded)
      endif
    endfor                      ; istate
  endif else begin
    ;; CHROMIS file names do not include the state in a useful form.
    ;; We could activate a similar method to the CRISP one above, once
    ;; the acquisition program does the wheel/hrz translation and puts
    ;; the true tuing info in the file name. (Or could we do some sort
    ;; of reverse translation from state to wheel/hrz?)
    pat = "*"
  endelse
  
  for iset = 0, Nsets-1 do begin

    timestamp = file_basename(dirs[iset])

    print, inam + ' : Working on '+timestamp

    outdir = self.out_dir +'/quicklook/'+timestamp+'/'
    file_mkdir, outdir
    
;    searchstring = dirs[iset] + '/' + cam + '/' + pat
    files = red_file_search(pat, dirs[iset] + '/' + cam + '/', count = Nfiles)

    if strmid(cam, 0, 5) eq 'Crisp' then begin
      ;; Check for lcd files!
      windx = where(~strmatch(files, '*.lcd.*'), Nwhere)
      if Nwhere eq 0 then files = '' else begin
        files = files[windx]
      endelse
      Nfiles = Nwhere
    endif
    
    if files[0] eq '' then begin
      print, inam + ' : ERROR -> no frames found in '+dirs[iset]
      continue
    endif
    
    self -> extractstates, files, states
    indx = uniq(states.tun_wavelength, sort(states.tun_wavelength))
    ustat = states[uniq(states.tun_wavelength, sort(states.tun_wavelength))].fullstate
    upref = states(uniq(states.prefilter, sort(states.prefilter))).prefilter
    Npref = n_elements(upref)

    if ~keyword_set(no_plot_r0) then begin

      ;; Plot r0

      cgcontrol, /delete, /all

      pname = states[0].camera
      pname += '_' + 'r0'
      pname += '_' + self.isodate
      pname += '_' + timestamp
      pname += '.pdf'

      if keyword_set(overwrite) || ~file_test(outdir+pname) then begin
        red_plot_r0_stats, states, pname = outdir+pname
      endif
      
    endif

    if n_elements(use_states) gt 0 then begin
      for istate = 0, n_elements(use_states)-1 do begin
        ;; Might want to change this into something involving strmatch
        ;; or stregex!
        print, use_states[istate]
        imatch = where(strmatch(ustat, '*'+use_states[istate]+'*'), Nmatch)
        print, Nmatch
        if Nmatch eq 0 then continue ; This ustat didn't match
        red_append, ustat2, ustat[imatch]
      endfor                    ; istate
      Nstates = n_elements(ustat2)
      if Nstates eq 0 then continue ; Next dataset
      ustat = ustat2                ; Go with the matching states 
    endif else begin    
      if n_elements(ustat) gt 1 then begin
        ;; Select states.
;        tmp = red_select_subset(ustat $
        tmp = red_select_subset(states[indx].prefilter+'_'+states[indx].tuning $
                                , qstring = inam + ' : Select states' $
                                , count = Nstates, indx = ichoice)
        ustat = ustat[ichoice]
      endif else begin
        Nstates = 1
      endelse
    endelse

    for istate = 0, Nstates-1 do begin
      
      self -> selectfiles, files = files, states = states $
                           , ustat = ustat[istate] $
                           , selected = sel

      uscan = states[uniq(states.scannumber,sort(states.scannumber))].scannumber
      Nscans = n_elements(uscan)

      sel = sel[sort(states[sel].scannumber)] ; Make sure scans are in order!

      ;; Load dark and flat
      self -> get_calib, states[sel[0]] $
                         , darkdata = dd, darkstatus = darkstatus $
                         , gaindata = gg, gainstatus = gainstatus
      if darkstatus ne 0 then dd = 0.
      if gainstatus ne 0 then gg = 1.
    
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
      namout += '.mp4'

      if ~keyword_set(overwrite) && file_test(outdir+namout) then continue

      print, inam + ' : saving to folder -> '+outdir 

      dim = size(dd, /dim)
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
      cube = fltarr(dim[0], dim[1], Nscans)
      med = fltarr(Nscans)
      best_contrasts = fltarr(Nscans)
      if gainstatus eq 0 then mask = red_cleanmask(gg ne 0)
      
      for iscan = 0L, Nscans -1 do begin

        red_progressbar, iscan, Nscans, 'Cube '+namout, /predict

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
            med = med[0:Nscans-1]
            best_contrasts = best_contrasts[0:Nscans-1]
            docontinue = 1
          end
          1 : begin
            ;; Just a single file
            ims = float(red_readdata(states[sel[iscan]].filename))
            Nframes = 1
          end
          else : begin
            ;; Multiple files, e.g. CHROMIS WB. Or any CRISP camera.
            ;; (Ideally, you should be able to call red_readdata with
            ;; multiple file names and if would do the right thing.
            ;; Needs to produce a combined header, too!)
            Nframes = 0L
            for ifile = 0, n_elements(sel2)-1 do begin
              Nframes += (fxpar(red_readhead(files[sel2[ifile]]),'NAXIS3') >1)
            endfor              ; ifile
            ims = fltarr(dim[0], dim[1], Nframes)
            iframe = 0
            for ifile = 0, n_elements(sel2)-1 do begin
              ims[0, 0, iframe] = red_readdata(files[sel2[ifile]], head = head)
              iframe += (fxpar(head,'NAXIS3') >1)
            endfor              ; ifile
          end
        endcase
        if docontinue then continue ; Not allowed in the case statement
        
        ;; Do dark and gain (and possibly backscatter) correction.
        for iframe = 0, Nframes-1 do begin
          ims[*, *, iframe] -= dd
          if DoBackscatter then begin ; Backscatter if needed
            ims[*, *, iframe] = red_cdescatter(ims[*, *, iframe] $
                                               , bgt, Psft $
                                               , verbose = verbose $
                                               , nthreads = nthreads)
          endif
          ims[*, *, iframe] *= gg
        endfor                  ; iframe

        if size(ims, /n_dim) eq 3 then begin
          ;; Possible improvement: base this selection on wb contrast?
          Nframes = (size(ims, /dim))[2]
          contrasts = fltarr(Nframes)
          for iframe = 0, Nframes-1 do contrasts[iframe] $
             = stddev(ims[x0:x1, y0:y1, iframe])/mean(ims[x0:x1, y0:y1, iframe])
          mx = max(contrasts, ml)
          im = ims[*, *, ml] 
          ;;im = (ims[*, *, ml] - dd) * gg
          best_contrasts[iscan] = contrasts[ml]
        endif else begin
          im = ims 
          ;;im = (ims - dd) * gg
          best_contrasts[iscan] = stddev(im[x0:x1, y0:y1])/mean(im[x0:x1, y0:y1])
        endelse

        ;; This should be done before calculating the contrast!
        if gainstatus eq 0 then im = red_fillpix(im, mask=mask, nthreads=nthreads)
        
        idx1 = where(im eq 0.0, Nwhere, complement = idx, Ncomplement = Nnowhere)
        if Nwhere gt 0 && Nnowhere gt 0 then im[idx1] = median(im[idx])
        
        if(keyword_set(x_flip)) then im = reverse(temporary(im), 1)
        if(keyword_set(y_flip)) then im = reverse(temporary(im), 2)
        
        cube[*, *, iscan] = im
        med[iscan] = median(cube[*, *, iscan])
        
      endfor                    ; iscan

      ;; Normalize intensities
      if Nscans gt 3 && ~keyword_set(no_normalize) then begin
        mm = mean(med)
        cc = poly_fit(findgen(Nscans), med/mm, 2, yfit=yfit)
        for iscan = 0, Nscans-1 do cube[*,*,iscan] /= yfit[iscan]*mm
      endif
      
      if keyword_set(align) || keyword_set(destretch) then begin
        
        ;; Measure image shifts
        shifts = red_aligncube(cube, 5, /center $ ;, cubic = -0.5 $
                               , xbd = round(dim[0]*.9) $
                               , ybd = round(dim[1]*.9) )

        ;; Outliers?
        indx_included = where((abs(shifts[0,*] - median(reform(shifts[0,*]),3)) le maxshift) $
                              and (abs(shifts[1,*] - median(reform(shifts[1,*]),3)) le maxshift) $
                              , complement = indx_excluded, Nincluded, Ncomplement = Nexcluded)
        if Nexcluded gt 0 then begin
          shifts[0, indx_excluded] = interpol(shifts[0, indx_included], indx_included, indx_excluded)
          shifts[1, indx_excluded] = interpol(shifts[1, indx_included], indx_included, indx_excluded)
        endif
        
        ;; Align the cube
        for iscan = 0, Nscans-1 do begin
          cube[*, *, iscan] = red_shift_im(cube[*, *, iscan] $
                                           , shifts[0, iscan] $
                                           , shifts[1, iscan] $
                                           , cubic = -0.5)
        endfor                  ; iscan
      endif

      if keyword_set(destretch) then begin

        ;; Not implemented yet!
        stop
        
        if n_elements(clip) eq 0 then clip = [12,  6,  3]
        if n_elements(tile) eq 0 then tile = [10, 20, 30]

        ;; Needs: array with time; rotate before align and stretch

        dts = red_time2double(time)
        if n_elements(tstep) eq 0 then begin
          tstep = fix(round(180. / median(abs(dts[0:Nscans-2] - dts[1:*]))))
        endif
        tstep = tstep < (Nscans-1)

        ;; Calculate stretch vectors
        grid = red_destretch_tseries(cube, 1.0/float(self.image_scale), tile, clip, tstep)

        for iscan = 0L, Nscans - 1 do begin
          red_progressbar, iscan, Nscans, inam+' : Applying the stretches.'
          cube[*,*,iscan] = red_stretch(cube[*,*,iscan], reform(grid[iscan,*,*,*]))
        endfor                  ; iscan
      endif
      
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
        red_fitspar_getdates, red_readhead(states[sel[iscan]].filename), date_avg = date_avg
        time_avg = (strsplit((strsplit(date_avg, 'T', /extract))[1], '.', /extract))[0]
        annstring = self.isodate + ' ' + time_avg $
                    + '   ' + ustat[istate] $
                    + '   scan :' + string(uscan[iscan],format='(I5)') 
        cgtext, [0.01], [0.95], annstring $
                , /normal, charsize=3., color=textcolor, font=1

        ;; tvrd(/true) reads an RGB image [3,Nx,Ny]
        snap2 = tvrd(/true)

        rgbcube[0, 0, 0, iscan] = snap2

      endfor                    ; iscan

      set_plot,'X'
      
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

      write_video, outdir+namout, rgbcube $
                   ;; , bit_rate=bitrate $
                   , metadata = metadata $
                   , video_fps = fps 

      ;; Make a jpeg image of the best frame. 
      mx = max(best_contrasts, ml)
      jname = outdir+red_strreplace(namout, '.mp4', '_scan='+strtrim(uscan[ml], 2)+'.jpg')

      write_jpeg, jname, rgbcube[*, *, *, ml], q = 100, /true

      print, outdir+namout
      print, jname
      print, strjoin(strtrim(size(cube, /dim), 2), ' x ')
      
    endfor                      ; istate
  endfor                        ; iset
  
end
