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
;    max_shift : in, optional, type=int, default=100
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
;    refcam : in, optional, type=integer, default=0
;
;      Select reference-channel.
;    
;    verbose : in, optional, type=integer, default=0
;
;      Provide more output
;
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
;   2016-06-13 : MGL. Removed "_thi" suffix. This is now the default
;                pinhole calibration method.
;
;-
pro red::pinholecalib, threshold = threshold $
                           , max_shift = max_shift $
                           , nref = nref $
                           , pref = pref $
                           , dir = dir $
                           , cams = cams $
                           , refcam = refcam $
                           , verbose = verbose

  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])
  
  ;; Logging
  help, /obj, self, output = selfinfo 
  red_writelog, selfinfo = selfinfo, logfile = logfile

  if(n_elements(threshold) eq 0) then threshold = 0.25
  if(n_elements(nref) eq 0) then nref = 5
  if( n_elements(dir) gt 0 ) then dir = [dir] $
  else if ptr_valid(self.pinh_dirs) then dir = *self.pinh_dirs
  self -> getdetectors, dir = dir
  if( n_elements(cams) gt 0 ) then cams = [cams] $
  else if ptr_valid(self.cameras) then cams = [*self.cameras]

  if(n_elements(refcam) eq 0) THEN refcam = self.refcam
  if(n_elements(verbose) eq 0) THEN verbose = 0
  if(n_elements(max_shift) eq 0) then max_shift = 200
  
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

  all_files = file_search( ph_dir + '*.pinh' )
  nf = n_elements( all_files )
  
  self->selectfiles, files=all_files, states=states, cam = cams[refcam] $
                     , /strip_settings, selected=selection
  ref_states = states[selection]
  ref_states_unique = ref_states.fpi_state
  ;; We discard multiple pinhole files for the same state
  ref_states_unique = ref_states[ uniq(ref_states_unique, sort(ref_states_unique)) ]
  for icam = 0, Ncams-1 do begin
    if icam EQ refcam then continue
    self->selectfiles, files=all_files, states=states, cam = cams[icam] $
                       , /strip_settings, selected=selection
    cam_states = states[selection]

    cam_states_unique = cam_states.fpi_state
    ;; We discard multiple pinhole files for the same state
    cam_states_unique = cam_states[ uniq(cam_states_unique, sort(cam_states_unique)) ] 
    
    last_prefilter = ''
    for iref=0, n_elements(ref_states_unique)-1 do begin
      cam_idx = where(cam_states.fpi_state eq ref_states_unique[iref].fpi_state)
      if n_elements(cam_idx) ne 1 || max(cam_idx) lt 0 then begin
        print, inam, ' : No (unique) pair of pinhole images.'
        continue
      endif
      if self->match_prefilters( ref_states_unique[iref].prefilter $
                                 , cam_states[cam_idx].prefilter ) eq 0 then continue

      this_prefilter = cam_states[cam_idx].prefilter
      if this_prefilter ne last_prefilter then undefine, this_init
      
      ref_fn = ref_states_unique[iref].filename
      print, 'loading reference file: ', ref_fn
      ref_img = red_readdata( ref_fn, /silent )
      
      cam_fn = cam_states[cam_idx].filename
      print, 'loading camera file: ', cam_fn
      cam_img = red_readdata( cam_fn, /silent )
      
      rdims = size( ref_img, /dim )
      cdims = size( cam_img, /dim )
      if( n_elements(rdims) ne 2 || n_elements(cdims) ne 2 ) then begin
        print, inam, ' : pinhole data is not 2-dimensional: ', ref_fn, cam_fn
        continue
      endif
      
      red_append, alignments, { state1:ref_states_unique[iref] $
                                , state2:cam_states[cam_idx] $
                                , map:rdx_img_align( ref_img, cam_img, nref=nref $
                                                     , h_init=this_init, threshold=threshold $
                                                     , verbose=verbose, max_shift=max_shift ) }

      last_prefilter = this_prefilter
    endfor                      ; iref
    
    okmaps = where( (alignments.state2.camera eq cams[icam]) and (alignments.map[2,2] eq 1) )
    failedmaps = where( (alignments.state2.camera eq cams[icam]) and (alignments.map[2,2] ne 1) )
    if (max(failedmaps) ge 0) and (max(okmaps) ge 0) then begin
      avg_map = alignments[okmaps].map
      if n_elements(okmaps) gt 1 then avg_map = total(avg_map,3)/n_elements(okmaps)
      for ifailed=0, n_elements(failedmaps)-1 do begin
        this_init = avg_map
        ref_img = red_readdata( alignments[failedmaps[ifailed]].state1.filename, /silent )
        img2 = red_readdata( alignments[failedmaps[ifailed]].state2.filename, /silent )
        alignments[failedmaps[ifailed]].map = rdx_img_align( ref_img, img2 $
                                                             , nref=nref, h_init=this_init $
                                                             , threshold=threshold $
                                                             , verbose=verbose )
      endfor                    ; ifailed
    endif

  endfor                        ; icam

  save, file = output_file, alignments
  
end
