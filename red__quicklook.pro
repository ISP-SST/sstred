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
;
;    datasets : in, optional, type=strarr
;
;      Timestamp strings that identify datasets to process. Selection
;      menu for data sets buypassed if given.
;   
;    destretch : in, optional, type=boolean
;
;      Set this to destretch the cube to compensate for geometrical
;      effects from anisoplanatism. Implies derotate and align.
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
;      Do neural net MFBD deconvolution. (So far only in Stockholm and
;      for CRISP data.)
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
;    overwrite :  in, optional, type=boolean
;
;      Overwrite existing movies.
;   
;    derotate : in, optional, type=boolean
;
;      Set this to derotate the cube to compensate for the field
;      rotation of the telescope.
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
;-
pro red::quicklook, align = align $
                    , bit_rate = bit_rate $
                    , cam = cam $
                    , clip = clip $
                    , core_and_wings = core_and_wings $
                    , dark = dark $
                    , datasets = datasets $
                    , destretch = destretch $
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
                    , overwrite = overwrite $
                    , remote_dir = remote_dir $
                    , remote_login = remote_login $
                    , derotate = derotate $
                    , ssh_find = ssh_find $
                    , textcolor = textcolor $
                    , use_states = use_states $
                    , verbose = verbose $
                    , video_codec = video_codec $
                    , video_fps = video_fps $
                    , x_flip = x_flip $
                    , y_flip = y_flip 
  
  inam = red_subprogram(/low, calling = inam1)

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
  
  ;; The r0 log file is not available until the day after today 
  if self.isodate eq (strsplit(red_timestamp(/utc,/iso),'T',/extract))[0] then no_plot_r0 = 1
  
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


  if strmid(cam, 0, 5) eq 'Crisp' then begin

    ;; CRISP data: we can't do this for CHROMIS data because
    ;; we need actual files for the conversion between wheel/hrz to
    ;; states. 
    
    if n_elements(use_states) gt 0 then begin
      ;; If use_states is provided, we base the file search pattern on
      ;; that. This works for CRISP, where the state is part of the file
      ;; name.
      
      ;; We need to do a bit of massaging here, because the CRISP states
      ;; in the file names have the tuning (after the sign) zero-padded
      ;; to 4 digits, while the states returned by extractstates (used
      ;; below) are not zero padded. We want this to work whether the
      ;; use_states are given zero padded or not. The pattern used for
      ;; file searching needs the padding but any padding has to be gone
      ;; when we get to the state comparison later. To make it even more
      ;; complicated, we don't know if this tuning is even part
      ;; of the use_states! 
      use_pat = strarr(n_elements(use_states))
      for istate = 0, n_elements(use_states)-1 do begin
        st = stregex(use_states[istate], '(_|\.|^)([+-][0-9]*)(_|\.|$)' $
                     , /extract,/sub)
        if st[2] ne '' then begin
          ;; We had a match for the tuning part of the state
          tun = st[2]
          tun_padded    = strmid(tun, 0, 1) + string(long(strmid(tun, 1)), format='(i04)')
          tun_nonpadded = strmid(tun, 0, 1) + string(long(strmid(tun, 1)), format='(i0)')
          use_pat[istate] = '*' + red_strreplace(use_states[istate], tun, tun_padded) + '*'
          use_states[istate] = red_strreplace(use_states[istate], tun, tun_nonpadded)
        endif
      endfor                    ; istate
    endif
  endif   

  ;; Now loop over datasets (timestamps)
  for iset = 0, Nsets-1 do begin

    timestamp = file_basename(dirs[iset])

    print, inam + ' : Working on '+timestamp

    outdir = self.out_dir +'/quicklook/'+timestamp+'/'
    file_mkdir, outdir
    

    ;; Try to limit the number of files we need to extract states for.

    if strmid(cam, 0, 5) eq 'Crisp' then begin

      ;; CRISP data

      ;; Search file names for scan 0, use them to find out what states
      ;; are available.
      files0 = red_file_search('*[_.]00000[_.]*', dirs[iset] + '/' + cam + '/', count = Nfiles)
      self -> extractstates, files0, states0
      indx = uniq(states0.tun_wavelength, sort(states0.tun_wavelength))
      ustat = states0[indx].fullstate
      upref = states0(uniq(states0.prefilter, sort(states0.prefilter))).prefilter
      Npref = n_elements(upref)

      
      if n_elements(use_pat) gt 0 then begin

        ;; use_pat constructed from use_states above
        files = red_file_search(use_pat, dirs[iset] + '/' + cam + '/', count = Nfiles)
        
      endif else begin

        if keyword_set(core_and_wings) then begin

          ustat_pat = reform((stregex(file_basename(states0[indx].filename) $
                                      , '\.([0-9][0-9][0-9][0-9]\.[0-9][0-9][0-9][0-9]_[+-][0-9]*\.lc[0-4])\.' $
                                      , /subex, /extract))[0,*])

          undefine, pat
          for ipref = 0, Npref-1 do begin
            sindx = where(strmatch(ustat, '*_'+upref[ipref]+'_*'), Nmatch)
            if Nmatch eq 1 then begin
              ;; If just one state for this prefilter, then use it!
              red_append, pat, ustat_pat[sindx[0]]
            endif else begin
              ;; Select red and blue wing points. The states are sorted in
              ;; wavelength order so we just have to pick the first and
              ;; last states for each prefilter.
              red_append, pat, ustat_pat[sindx[0]]
              red_append, pat, ustat_pat[sindx[-1]]
              ;; Find and select the core.
              imatch = where(strmatch(ustat[sindx], '*+0*'), Nmatch)
              if Nmatch gt 0 then red_append, pat, ustat_pat[sindx[imatch]]
            endelse
          endfor                ; ipref

          pat = '*'+pat+'*'
          
          files = red_file_search(pat, dirs[iset] + '/' + cam + '/', count = Nfiles)

        endif else begin        

          files = red_file_search('*', dirs[iset] + '/' + cam + '/', count = Nfiles)

        end
      end

      ;; Check for lcd files!
      windx = where(~strmatch(files, '*.lcd.*'), Nwhere)
      if Nwhere eq 0 then files = '' else begin
        files = files[windx]
      endelse
      Nfiles = Nwhere
      
    endif else begin

      ;; CHROMIS data

      if keyword_set(core_and_wings) then begin

        ;; Search file names for scan 0, use them to find out what states
        ;; are available.
        files0 = red_file_search('*[_.]00000[_.]*', dirs[iset] + '/' + cam + '/', count = Nfiles)
        self -> extractstates, files0, states0
        indx = uniq(states0.tun_wavelength, sort(states0.tun_wavelength))
        ustat = states0[uniq(states0.tun_wavelength, sort(states0.tun_wavelength))].fullstate
        upref = states0(uniq(states0.prefilter, sort(states0.prefilter))).prefilter
        Npref = n_elements(upref)

        undefine, ustat2
        for ipref = 0, Npref-1 do begin
          sindx = where(strmatch(ustat, '*_'+upref[ipref]+'_*'), Nmatch)
          if Nmatch eq 1 then begin
            ;; If just one state for this prefilter, then use it!
            red_append, ustat2, ustat[sindx[0]]
          endif else begin
            ;; Select red and blue wing points. The states are sorted in
            ;; wavelength order so we just have to pick the first and
            ;; last states for each prefilter.
            red_append, ustat2, ustat[sindx[ 0]]
            red_append, ustat2, ustat[sindx[-1]]
            ;; Find and select the core.
            imatch = where(strmatch(ustat[sindx], '*+0*'), Nmatch)
            if Nmatch gt 0 then red_append, ustat2, ustat[sindx[imatch]]
          endelse
        endfor                  ; ipref

        ;; CHROMIS file names do not include the state in a
        ;; human-readable form. But we can see what states are available
        ;; for scan 0, match the wanted states, and figure out what the
        ;; file name states are.
        for istate = 0, n_elements(ustat2)-1 do begin
          imatch = where(ustat2[istate] eq states0.fullstate, Nmatch)
          if Nmatch ge 1 then ustat2[istate] = states0[imatch[0]].fpi_state
        endfor                  ; istate
        
        pat = '*_'+ustat2+'.fits'
        files = red_file_search(pat, dirs[iset] + '/' + cam + '/', count = Nfiles)

      endif else begin

        files = red_file_search('*', dirs[iset] + '/' + cam + '/', count = Nfiles)

      endelse
      
    endelse
    
    if files[0] eq '' then begin
      print, inam + ' : ERROR -> no frames found in '+dirs[iset]
      continue                  ; Goto next dataset
    endif
    
    self -> extractstates, files, states

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

    if keyword_set(core_and_wings) then begin
      ;; Make an automatic selection of states
      undefine, ustat2
      for ipref = 0, Npref-1 do begin
        sindx = where(strmatch(ustat, '*_'+upref[ipref]+'_*'), Nmatch)
        if Nmatch eq 1 then begin
          ;; E.g., Chromis Ca II core should be included
          red_append, ustat2, ustat[sindx[0]]
        endif else begin
          ;; Select red and blue wing points. The states are sorted in
          ;; wavelength order so we just have to pick the first and
          ;; last states for each prefilter.
          red_append, ustat2, ustat[sindx[ 0]]
          red_append, ustat2, ustat[sindx[-1]]
          ;; Find and select the core.
          imatch = where(strmatch(ustat[sindx], '*+0*'), Nmatch)
          if Nmatch gt 0 then red_append, ustat2, ustat[sindx[imatch]]
        endelse
      endfor                    ; ipref
      Nstates = n_elements(ustat2)
      if Nstates eq 0 then continue ; Next dataset
      ustat = ustat2 
    endif else if n_elements(use_states) gt 0 then begin
      undefine, ustat2
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
      if keyword_set(neuralnet) then namout += '_NN' 
      if keyword_set(mtf_deconvolve) then namout += '_MTF' 
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

      print, inam + ' : Cube '+namout
      
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
            endfor              ; ifile
            ims = fltarr(dim[0], dim[1], Nframes)
            iframe = 0
            for ifile = 0, n_elements(sel2)-1 do begin
;              print, 'Read '+file_basename(files[sel2[ifile]])
              ims[0, 0, iframe] = red_readdata(files[sel2[ifile]], head = head)
              iframe += (fxpar(head,'NAXIS3') >1)
            endfor              ; ifile
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

          im = red_fillpix(ims[*, *, iframe], mask=mask, nthreads=nthreads)

          idx1 = where(im eq 0.0, Nwhere, complement = idx, Ncomplement = Nnowhere)
          if Nwhere gt 0 && Nnowhere gt 0 then im[idx1] = median(im[idx])

          if(keyword_set(x_flip)) then im = reverse(temporary(im), 1)
          if(keyword_set(y_flip)) then im = reverse(temporary(im), 2)

          ims[*, *, iframe] = im

        endfor                  ; iframe
        
        for iframe = 0, Nframes-1 do contrasts[iframe] $
           = stddev(ims[x0:x1, y0:y1, iframe])/mean(ims[x0:x1, y0:y1, iframe])

        if keyword_set(neuralnet) then begin

          ;; Deconvolve the data with a neural net.

          ;; Select the best frames
          nn_indx = reverse(sort(contrasts))
          cube[0, 0, 0, iscan] =ims[*, *, nn_indx[0:(Nexp <Nexp_available)-1]]
          
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

      endfor                    ; iscan

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
        
        caminfo = red_camerainfo(detector)
        
        lambda = states[sel2[0]].tun_wavelength ; Wavelength [m]
        telescope_d = 0.97d
        arcsecperpix = self.image_scale
        pixelsize = caminfo.pixelsize
        sz = max(dim)
        
        F_number = pixelsize/telescope_d/(arcsecperpix*2d*!dpi/(360.*3600.))
        Q_number = F_number * lambda/pixelsize

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
              
          fim = fft(w*im)       ; FFT of windowed image
          
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
          cube[*,*,iscan] = red_rotation(cube[*,*,iscan], ang[iscan])
        endfor                  ; iscan
      endif

      if keyword_set(align) || keyword_set(destretch) then begin
        ;; Align images
        
        ;; Measure image shifts
        shifts = red_aligncube(cube, 5, /center $ ;, cubic = -0.5 $
                               , xbd = round(dim[0]*.9) $
                               , ybd = round(dim[1]*.9) )

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
                                           , cubic = -0.5)
        endfor                  ; iscan
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
        grid = red_destretch_tseries(cube, 1.0/float(self.image_scale), tiles, clips, tstep)

        for iscan = 0L, Nscans - 1 do begin
          red_progressbar, iscan, Nscans, 'Applying the stretches.'
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


;      dims = size(rgbcube, /dim)
;      rgbcube = byte(255*randomu(seed, [3, 100, 100, 13]))

;      help, bit_rate, video_fps

      file_delete, outdir+namout, /allow
      write_video, outdir+namout, rgbcube $
                   , bit_rate = bit_rate $
                   , metadata = metadata $
                   , video_fps = video_fps $
                   , video_codec = video_codec 
      

;      stop
      
      ;; Make a jpeg image of the best frame. 
      mx = max(best_contrasts, ml)
      jname = outdir+red_strreplace(namout, '.'+extension, '_scan='+strtrim(uscan[ml], 2)+'.jpg')

      file_delete, jname, /allow
      write_jpeg, jname, rgbcube[*, *, *, ml], q = 100, /true

      print, outdir+namout
      print, jname
      print, strjoin(strtrim(size(cube, /dim), 2), ' x ')

      if format eq 'mov' then begin
        ;; Convert to Mac-friendly (and smaller) .mov file using recipe from Tiago
        mname = outdir + red_strreplace(namout, '.'+extension,'.'+format)
        spawn, 'ffmpeg -n -i "' + outdir + namout $
               + '" -c:v libx264 -preset slow -crf 26 -vf scale=-1:800  -tune grain "' $
               + mname + '"'
        spawn, 'rm "' + outdir + namout + '"'
;        find . -name '*mp4' -exec sh -c 'ffmpeg -n -i "$1" -c:v libx264 -preset slow -crf 26 -vf scale=-1:800  -tune grain "${1%.mp4}.mov"' sh {} \ ;
      endif
      
    endfor                      ; istate
  endfor                        ; iset
  
end
