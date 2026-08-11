; docformat = 'rst'

;+
;   Calculate align-clip and offset files for an arbitrary
;   subset of cameras and states.
;
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;     Tomas Hillberg, Institute for Solar Physics, 2016
;
; 
; :Params:
; 
;   
; 
; :Keywords:
;    
;    align : out
;
;      Array of structures containing clip, offset filenames and
;      state information.
;
;    prefilters : in, optional, type=strarr
;
;      Indicate the prefilters to calculate the clips/offsets for.
;      Default is to do it for all prefilters.
;
;    extraclip : in, optional, type=intarr(4)
;
;      Exclude a margin around the clip-area
;
;    cams : in, optional, type=strarr
;
;      Indicate the cameras to calculate the clips/offsets for.
;      Default is to do it for all cameras.
;
;    refcam : in, optional, type=integer, default=0
;
;      Select reference camera
;
;    overwrite : in, optional, type=boolean, default=false
;
;      Force re-calculation if files already exist.
;
;    output_dir : in, optional, type=string, default="self.out_dir+'/calib/'"
;
;      The directory in which to write offset files.
; 
; :History:
;
;   2016-09-21 : First version
;
;   2017-01-18 : MGL. New keyword output_dir.
;
;   2017-04-10 : MGL. New keyword makeoffsets.
;
;-
pro red::getalignment, align = align, $
                       prefilters = prefilters, $
                       extraclip = extraclip, $
                       cams = cams, $
                       refcam = refcam, $
                       makeoffsets = makeoffsets, $
                       output_dir = output_dir, $
                       overwrite = overwrite

  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)              
  
  if( n_elements(cams) gt 0 ) then cams = [cams] $
  else if ptr_valid(self.cameras) then cams = [*self.cameras]

  if(n_elements(refcam) eq 0) THEN refcam = self.refcam
  if(n_elements(verbose) eq 0) THEN verbose = 0

  if n_elements(output_dir) eq 0 then output_dir = self.out_dir+'/calib/'
  alignfile = self.out_dir+'/calib/alignments.sav'

  if file_test(alignfile) then begin
    restore, alignfile
  endif else begin
    print, inam + ' : Alignment file missing.'
    print, inam + ' : Did you run a -> pinholecalib ?'
    return
  endelse

  Ncams = n_elements(cams)
  if Ncams lt 2 then begin
    print, inam, ' : Need at least 2 cameras.'
    print, inam, ' : cams: ',cams
    print, inam, ' : self.cameras: ',*self.cameras
    return
  endif
  
  if refcam ge Ncams then begin
    print, inam, ' : index of reference camera out of range: ', refcam, ' >= ', Ncams
    return
  endif
  
  refcam_name = cams[refcam]

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

  ;; filter selection
  Nalign = n_elements(alignments)
  Npref = n_elements(prefilters)
  for ialign=0, Nalign-1 do begin
    skip = (max(where(cams eq alignments[ialign].state1.camera)) lt 0 $
            or max(where(cams eq alignments[ialign].state2.camera)) lt 0)
    if skip eq 0 and Npref gt 0 then begin
      skip = (max(where(prefilters eq alignments[ialign].state1.prefilter)) lt 0 $
              and max(where(prefilters eq alignments[ialign].state2.prefilter)) lt 0)
    endif
    if skip eq 0 then red_append, selected, ialign
  endfor

  alignments = alignments[ selected ]
  Nalign = n_elements(alignments)
  if Nalign eq 0 then begin
    print, inam, ' : No alignments left, was the filtering too strict?'
    return
  endif
  
  allstates = [[alignments.state1],[alignments.state2]]
  ref_idx = where( allstates.camera eq refcam_name )
  if max(ref_idx) lt 0 then begin
    print, inam, ' : Alignment not found for reference camera: ', refcam_name
    return
  endif

  pos = strpos(allstates[ref_idx[0]].filename,'pinhs/')
  head = red_readhead(strmid(allstates[ref_idx[0]].filename, pos), /silent )
  ref_dims = fxpar(head, 'NAXIS*')
  ref_corners = red_make_corners( [ 0, ref_dims[0]-1 , $
                                    0, ref_dims[1]-1 ] )

  ;; TODO invert, or re-map, if refcam is not the same as when pinholecalib was run.
  halfsize = ref_dims/2
  common_fov = ref_corners

  for ialign=0, Nalign-1 do begin

    pos = strpos(alignments[ialign].state2.filename,'pinhs/')
    head = red_readhead(strmid(alignments[ialign].state2.filename, pos), /silent )
    dims = fxpar(head, 'NAXIS*')
    corners = red_make_corners( [ 0, dims[0]-1 , $
                                  0, dims[1]-1 ] )
                                ; transform the corners of the FOV to the non-reference camera
    common_fov = common_fov # alignments[ialign].map
    idx = where(common_fov(*,2) ne 0, COMPLEMENT=idx_c)
    if max(idx) ne -1 then begin
      common_fov(idx,0) /= common_fov(idx,2)
      common_fov(idx,1) /= common_fov(idx,2)
    endif
    if max(idx_c) ne -1 then common_fov(idx_c,*) = 0
    common_fov(*,2) = 1
    
                                ; if a mapped corner falls outside the range, clip it
    idx = where(common_fov(*,0) lt corners(0,0) )
    if max(idx) ne -1 then  common_fov(idx,0) = corners(0,0)
    idx = where(common_fov(*,0) gt corners(1,0) )
    if max(idx) ne -1 then  common_fov(idx,0) = corners(1,0)
    idx = where(common_fov(*,1) lt corners(0,1) )
    if max(idx) ne -1 then  common_fov(idx,1) = corners(0,1)
    idx = where(common_fov(*,1) gt  corners(2,1) )
    if max(idx) ne -1 then  common_fov(idx,1) = corners(2,1)
    
                                ; transform back to reference camera
    common_fov = common_fov # invert(alignments[ialign].map)
    idx = where(common_fov(*,2) ne 0, COMPLEMENT=idx_c)
    if max(idx) ne -1 then begin
      common_fov(idx,0) /= common_fov(idx,2)
      common_fov(idx,1) /= common_fov(idx,2)
    endif
    if max(idx_c) ne -1 then common_fov(idx_c,*) = 0
    common_fov(*,2) = 1
    
                                ; clip again
    idx = where( common_fov(*,0) lt ref_corners(0,0) )
    if max(idx) ne -1 then  common_fov(idx,0) = ref_corners(0,0)
    idx = where( common_fov(*,0) gt ref_corners(1,0) )
    if max(idx) ne -1 then  common_fov(idx,0) = ref_corners(1,0)
    idx = where( common_fov(*,1) lt ref_corners(0,1) )
    if max(idx) ne -1 then  common_fov(idx,1) = ref_corners(0,1)
    idx = where( common_fov(*,1) gt ref_corners(2,1) )
    if max(idx) ne -1 then  common_fov(idx,1) = ref_corners(2,1)
    
  endfor

                                ; round inwards
  idx = where( common_fov[*,0] lt halfsize[0], compl=cidx )
  common_fov[idx,0] = ceil(common_fov[idx,0] + extraclip[0])
  common_fov[cidx,0] = floor(common_fov[cidx,0] - extraclip[1])
  idx = where( common_fov[*,1] lt halfsize[1], compl=cidx )
  common_fov[idx,1] = ceil(common_fov[idx,1] + extraclip[2])
  common_fov[cidx,1] = floor(common_fov[cidx,1] - extraclip[3])
  
  common_fov = fix(common_fov)
  
  ref_clip = intarr(4)
  ref_clip[0] = max(common_fov([0,2],0))
  ref_clip[1] = min(common_fov([1,3],0))
  ref_clip[2] = max(common_fov([0,1],1))
  ref_clip[3] = min(common_fov([2,3],1))

  ref_mid = [ (ref_clip[0]+ref_clip[1])/2.0, (ref_clip[2]+ref_clip[3])/2.0, 1 ]
  ref_sz = [ abs(ref_clip[0]-ref_clip[1])+1, abs(ref_clip[2]-ref_clip[3])+1 ]
  ref_origin = [ min(ref_clip[0:1]), min(ref_clip[2:3]) ]
  
  sx = long(ref_sz[0])
  sy = long(ref_sz[1])

  indices = [ [[dindgen(sx)#replicate(1.d0, sy) + ref_origin(0)]], $
              [[replicate(1.d0, sx)#dindgen(sy) + ref_origin(1)]], $
              [[replicate(1.d0, sx, sy)]]]

  ref_align = { clip:ref_clip, $
                state1:alignments[0].state1, $
                state2:alignments[0].state1, $
                map:identity(3), $
                xoffs_file:'', yoffs_file:'' }
  align = replicate( ref_align, Nalign+1 )

  for ialign=0, Nalign-1 do begin
    
    ;; map center and generate a clip around it
    h = alignments[ialign].map
    mid = ref_mid # alignments[ialign].map
    if abs(mid[2]) gt 1E-6 then begin
      mid /= mid[2]
    endif
    
    align[ialign+1].clip[0] = round(mid[0]-ref_sz[0]/2)
    align[ialign+1].clip[1] = align[ialign+1].clip[0] + ref_sz[0] - 1
    align[ialign+1].clip[2] = round(mid[1]-ref_sz[1]/2)
    align[ialign+1].clip[3] = align[ialign+1].clip[2] + ref_sz[1] - 1
    align[ialign+1].state1 = alignments[ialign].state1
    align[ialign+1].state2 = alignments[ialign].state2
    align[ialign+1].map = h

    if keyword_set(makeoffsets) then begin

      fname = strjoin([ align[ialign+1].state1.detector, $
                        align[ialign+1].state1.fullstate, $
                        align[ialign+1].state2.detector, $
                        align[ialign+1].state2.fullstate ],'_')
      align[ialign+1].xoffs_file = output_dir + fname + '.xoffs'
      align[ialign+1].yoffs_file = output_dir + fname + '.yoffs'
      flipped = transpose([0,0,1]) # h gt mid
      if flipped[0] gt 0 then begin
        align[ialign+1].clip[0:1] = reverse(align[ialign+1].clip[0:1])
      endif
      if flipped[1] gt 0 then begin
        align[ialign+1].clip[2:3] = reverse(align[ialign+1].clip[2:3])
      endif

      if( keyword_set(overwrite) || ~file_test(align[ialign+1].xoffs_file) $
          || ~file_test(align[ialign+1].xoffs_file)) then begin

        red_make_offs, h, xoff, yoff, align[ialign+1].clip+1, ref_clip=ref_clip+1
        
        print, 'Saving files ' + fname +'.(x|y)offs'
        red_writedata, align[ialign+1].xoffs_file, xoff $
                       , filetype='ANA', /overwrite
        red_writedata, align[ialign+1].yoffs_file, yoff $
                       , filetype='ANA', /overwrite
      endif
      
    endif

  endfor                        ; ialign
  
  align.clip += 1               ; return 1-based clips

end
