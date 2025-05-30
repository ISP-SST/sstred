; docformat = 'rst'

;+
;   Find the transformation matrices that maps the reference channel
;   onto the other ones.
;
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;     Tomas Hillberg, Institute for Solar Physics, 2015
;   
; 
; :Keywords:
;    
;    threshold : in, optional, type=float, default=0.25
;
;       Threshold for identifying a strong enough pinhole.
;
;    margin : in, optional, type=integer, default=30
;
;       A margin (in pixels) to remove from the detector FOV.
;
;    max_shift : in, optional, type=int, default=200
;
;       Only consider mappings with a linear shift < max_shift pixels.
;
;    nref : in, optional, type=integer, default=5
;
;      How many of the strongest pinholes to use for finding the
;      approximate transform. Afterwards a refinement is made using
;      >80% of the detected pinholes.
;
;    pref : in, optional, type=string
;
;      Indicate the prefilter you want to calculate the clips for,
;      Default is to do it for all prefilters there is data for.
;
;    dir : in, optional, type=strarr
;
;      Restrict to pinholes in the listed directories. Default
;      is to use all directories listed in *(self.pinh_dirs).
;
;    cams : in, optional, type=strarr
;
;      Restrict to specified cameras. Default is to use all
;      cameras listed in *(self.cameras).
;
;    no_polywarp ; in, optional, type=boolean
;
;      Do not use the polywarp model for camera alignments.   
;
;    refcam : in, optional, type=integer, default=0
;
;      Select reference-channel.
;    
;    verbose : in, optional, type=integer, default=0
;
;      Provide more output
;
;    verify : in, optional, type=boolean
;
;      Turn on tool for visually confirming the calibration results. 
; 
; :History:
;
;   2015-12-01 : New implementation that uses OpenCV functionality for
;                finding and aligning the pinholes. The offset files
;                are now directly computed instead of fitted.
;
;   2016-06-13 : MGL. Various cosmetic edits. Change bunch of if
;                statements to a case statement. Move make_corners to
;                the red_ namespace and its own file. Bugfix.
;
;   2016-08-23 : THI. Rename camtag to detector and channel to camera,
;                so the names match those of the corresponding
;                SolarNet keywords.
;
;   2016-09-21 : THI. Re-write to support matching of prefilters
;                between NB/WB channels. Some cleanup. Move generation
;                of clips/offset-files to a separate method
;                (red::getalignment)
;
;   2017-03-07 : MGL. Removed "_thi" suffix. This is now the default
;                pinhole calibration method.
;
;   2017-04-10 : MGL. Changed threshold default to 0 due to
;                reimplementation of rdx_img_align.
;
;   2017-07-19 : THI. Change pinholecalib parameter nref defaut value
;                to 10.
;
;   2017-10-04 : MGL. Add progressbar.
;
;   2017-12-06 : THI. Added a simple tool for verifying the calibrations.
;
;   2017-12-20 : THI. Add keyword smooth to pinholecalib, which is
;                passed on to rdx_img_align. Pass the pref keyword to
;                selectfiles.
;
;   2017-12-21 : THI. Print a warning, and some suggestions how to
;                proceed, if the pinhole calibration fails for some
;                state(s).
;
;   2018-10-05 : MGL. New keyword margin.
;
;   2019-06-18 : THI. Don't overwrite existing results when processing a
;                different prefilter. New keyword overwrite.
;
;   2020-10-29 : JdlCR propagate more keywords to phverify, allow the
;                user to provide an external initialization
;
;   2025-02-18 : MGL. Call pinholecalib2 for polywarp method for data
;                with new pinhole array. Also new keyword no_polywarp.
;
;-
pro red::pinholecalib, cams = cams $
                       , dir = dir $
                       , init = init $
                       , margin = margin $
                       , max_shift = max_shift $
                       , nref = nref $
                       , no_polywarp = no_polywarp $
                       , overwrite = overwrite $
                       , pref = pref $
                       , refcam = refcam $
                       , smooth = smooth $
                       , threshold = threshold $
                       , verbose = verbose $
                       , verify = verify 

  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)  

  if self.isodate gt red_dates(tag = 'pinhole array with L') && $
     ~keyword_set(no_polywarp) then begin
    ;; When the new pinhole array is installed, switch to using the
    ;; new polywarp model by default.
    self -> pinholecalib2 $
       , avg_fpistates = avg_fpistates $ ; Keyword not used
       , cams=cams $
       , dir=dir $
       , manual=manual $
       , margin=margin $
       , overwrite=overwrite $
       , pref=pref $
       , refcam=refcam $
       , threshold=threshold $
       , verbose=verbose $      ; Keyword not used
       , verify=verify $        ; Keyword not used
       , warp_dim=warp_dim
    return
  endif
  
  ;; Logging
  help, /obj, self, output = selfinfo 
  red_writelog, selfinfo = selfinfo, logfile = logfile

  if(n_elements(threshold) eq 0) then threshold = 0.0
  if(n_elements(nref) eq 0) then nref = 10
  if( n_elements(dir) gt 0 ) then dir = [dir] $
  else if ptr_valid(self.pinh_dirs) then dir = *self.pinh_dirs
  self -> getdetectors, dir = dir
  if( n_elements(cams) gt 0 ) then cams = [cams] $
  else if ptr_valid(self.cameras) then cams = [*self.cameras]

  if(n_elements(refcam) eq 0) THEN refcam = self.refcam
  if(n_elements(verbose) eq 0) THEN verbose = 0
  if(n_elements(max_shift) eq 0) then max_shift = 200
  if n_elements(margin) eq 0 then margin = 30
  
  Ncams = n_elements(cams)

  case n_elements(extraclip) of
    0 : extraclip = [0L, 0L, 0L, 0L]
    1 : extraclip = replicate(extraclip, 4)
    2 : extraclip = [ replicate(extraclip[0], 2), $
                      replicate(extraclip[1], 2) ]
    4 :                         ; Leave as it is.
    else : begin
      print, inam + "ERROR: Don't know how to use keyword extraclip with " $
             + strtrim(n_elements(extraclip), 2) + ' elements.'
      stop
    end
  endcase

  ;; Search summed pinh images and camtag
  ph_dir = self.out_dir+'/pinhs/'
  ph_dir = red_strreplace(ph_dir,'//','/')
  
  output_dir = self.out_dir+'/calib/'
  output_dir = red_strreplace(output_dir,'//','/')
  output_file = output_dir+'alignments.sav'
  file_mkdir, output_dir
  
  if Ncams lt 2 then begin
    print, inam, ' : Need at least 2 cameras.'
    print, inam, ' : dir: ',dir
    print, inam, ' : cams: ',cams
    print, inam, ' : self.cameras: ',*self.cameras
    red_writelog, /add, logfile = logfile, top_info_strings = ' : Need at least 2 cameras.'
    return
  endif

  if refcam ge Ncams then begin
    print, inam, ' : index of reference camera out of range: ', refcam, ' >= ', Ncams
    return
  endif
  
  if( file_test(output_file) ) then begin
    print, inam, ' : loading existing alignments: ', output_file
    restore, output_file
  endif


  all_files = file_search( ph_dir + '*.pinh.fits' )
  nf = n_elements( all_files )
  
  self->selectfiles, files=all_files, states=states, cam = cams[refcam] $
                     , /strip_settings, selected=selection, prefilter=pref
  ref_states = states[selection]
  ref_states_unique = ref_states.fpi_state
  ;; We discard multiple pinhole files for the same state
  ref_states_unique = ref_states[ uniq(ref_states_unique, sort(ref_states_unique)) ]
  for icam = 0, Ncams-1 do begin
    if icam EQ refcam then continue
    self->selectfiles, files=all_files, states=states, cam = cams[icam] $
                       , /strip_settings, selected=selection, prefilter=pref, count = Nselect
    if Nselect eq 0 then continue

    cam_states = states[selection]

    cam_states_unique = cam_states.fpi_state
    ;; We discard multiple pinhole files for the same state
    cam_states_unique = cam_states[ uniq(cam_states_unique, sort(cam_states_unique)) ] 
    
    last_prefilter = ''
    for iref=0, n_elements(ref_states_unique)-1 do begin

      red_progressbar, iref, n_elements(ref_states_unique), /predict, ''
      
      cam_idx = where(cam_states.fpi_state eq ref_states_unique[iref].fpi_state)
      if n_elements(cam_idx) ne 1 || max(cam_idx) lt 0 then begin
        print, inam, ' : No (unique) pair of pinhole images.'
        continue
      endif
      
      if self->match_prefilters( ref_states_unique[iref].prefilter $
                                 , cam_states[cam_idx].prefilter ) eq 0 then continue

      this_pair = ref_states_unique[iref].fullstate+cam_states[cam_idx].fullstate
      this_pair += ref_states_unique[iref].camera+cam_states[cam_idx].camera
      if n_elements( processed_pairs ) gt 0 then begin
        found = where( processed_pairs eq this_pair )
        if max(found) ge 0 then begin
          ;; pair is already processed during this run of pinholecalib
          continue
        endif
      endif
      
      align_idx = -1
      if n_elements( alignments ) gt 0 then begin
        align_idx = where( (alignments.state1.fullstate eq ref_states_unique[iref].fullstate) and $
                           (alignments.state2.fullstate eq cam_states[cam_idx].fullstate) and $
                           (alignments.state1.camera eq ref_states_unique[iref].camera) and $
                           (alignments.state2.camera eq cam_states[cam_idx].camera) )
      endif
      if max(align_idx) lt 0 then begin
        ;; Pair does not exist in alignments, append it.
        align_idx = n_elements( alignments )
        red_append, alignments, { state1:ref_states_unique[iref], state2:cam_states[cam_idx], map:fltarr(3,3) }
      endif else begin
        ;; pair already exists: Run again if requested, or if the existing map is invalid
        if( ~keyword_set(overwrite) and alignments[align_idx].map[2,2] eq 1 ) then continue
      endelse
      
      this_prefilter = cam_states[cam_idx].prefilter
      if this_prefilter ne last_prefilter then begin
        undefine, this_init
        if(n_elements(init) eq 9) then this_init = reform(init, [3,3])
      endif
      
      ref_fn = ref_states_unique[iref].filename
      red_progressbar, iref, n_elements(ref_states_unique), /predict $
                       , 'load ref im '+file_basename(ref_fn)
;      print, 'loading reference file: ', ref_fn
      ref_img = red_readdata( ref_fn, /silent )
      
      cam_fn = cam_states[cam_idx].filename
      red_progressbar, iref, n_elements(ref_states_unique), /predict $
                       , 'load cam im '+file_basename(cam_fn)
;      print, 'loading camera file: ', cam_fn
      cam_img = red_readdata( cam_fn, /silent )
      
      rdims = size( ref_img, /dim )
      cdims = size( cam_img, /dim )
      if( n_elements(rdims) ne 2 || n_elements(cdims) ne 2 ) then begin
        print, inam, ' : pinhole data is not 2-dimensional: ', ref_fn, cam_fn
        continue
      endif
      
      this_map = rdx_img_align( ref_img $
                                , cam_img $
                                , nref=nref $
                                , h_init=this_init $
                                , threshold=threshold $
                                , smooth=smooth $
                                , verbose=verbose $
                                , max_shift=max_shift $
                                , margin = margin )
      
      if (last_prefilter ne this_prefilter) or (this_map[2,2] ne 1) then begin
        if keyword_set(verify) then begin
          this_map = red_phverify( ref_img, cam_img, this_map, nref=nref, threshold=threshold $
                                   , max_shift=max_shift, margin=margin)
          this_init = this_map
        endif
      endif else begin
        ; TODO: sanity check for maps within the same prefilter.
        ; e.g. check how different this_map is from this_init ??
      endelse
      alignments[align_idx].map = this_map
      last_prefilter = this_prefilter
      red_append, processed_pairs, this_pair
    endfor                      ; iref
    
    okmaps     = where( (alignments.state2.camera eq cams[icam]) and (alignments.map[2,2] eq 1) )
    failedmaps = where( (alignments.state2.camera eq cams[icam]) and (alignments.map[2,2] ne 1) )
    if (max(failedmaps) ge 0) and (max(okmaps) ge 0) then begin
      avg_map = alignments[okmaps].map
      if n_elements(okmaps) gt 1 then avg_map = total(avg_map,3)/n_elements(okmaps)
      for ifailed=0, n_elements(failedmaps)-1 do begin
        this_init = avg_map
        ref_img = red_readdata( alignments[failedmaps[ifailed]].state1.filename, /silent )
        img2 = red_readdata( alignments[failedmaps[ifailed]].state2.filename, /silent )
        alignments[failedmaps[ifailed]].map = rdx_img_align( ref_img $
                                                             , img2 $
                                                             , nref=nref $
                                                             , h_init=this_init $
                                                             , threshold=threshold $
                                                             , smooth=smooth $
                                                             , verbose=verbose $
                                                             , max_shift=max_shift $
                                                             , margin = margin )
      endfor                    ; ifailed
    endif

    failedmaps = where( (alignments.state2.camera eq cams[icam]) and (alignments.map[2,2] ne 1) )
    if max(failedmaps) ge 0 then begin
      LF = string(10b)
      print, LF + 'Failed to match pinholes for the following states:'
      for ifailed=0, n_elements(failedmaps)-1 do begin
        msg = '    alignments[' + strtrim(failedmaps[ifailed],2) + ']: '
        msg += alignments[failedmaps[ifailed]].state1.detector + ':' + alignments[failedmaps[ifailed]].state1.tuning + ' <-> '
        msg += alignments[failedmaps[ifailed]].state2.detector + ':' + alignments[failedmaps[ifailed]].state2.tuning
        print, msg
      endfor                    ; ifailed
      print, LF + 'Have a look at the pinholes/mapping that they look sane. For example:' + LF
      print, 'restore,"calib/alignments.sav"'
      print, 'ind = ' + strtrim( max(failedmaps), 2 )
      print, 'ph1 = red_readdata( alignments[ind].state1.filename )'
      print, 'ph2 = red_readdata( alignments[ind].state2.filename )'
      print, 'print, alignments[ind].map'
      print, 'sz = size(ph2,/dim)'
      print, 'window, 0, xs=sz[0], ys=sz[1]'
      print, 'tvscl, ph2'
      print, 'window, 1, xs=sz[0], ys=sz[1]'
      print, 'tvscl, rdx_img_project( alignments[ind].map, ph1 )'
      print, 'blink, [0,1]'
      print, LF + 'If the images look fine, verify that the mapping-matrix looks approximately like:'
      print, ' ~\pm{1}     ~0     x-origin' 
      print, '   ~0     ~\pm{1}   y-origin' 
      print, '   ~0        ~0        1' + LF
      print, 'If not, re-run pinholecalib with the /verify keyword.'
      print, 'This will open a simple GUI that can be used to tweak the parameters until the calibration looks ok (i.e. pinholes overlap).'
    endif
  endfor                        ; icam

  save, file = output_file, alignments
  
end
