; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
; 
; 
; 
; 
; :Keywords:
; 
;    date_obs : in, optional, type=string
;   
;      The date of observations in ISO (YYYY-MM-DD) format. If not
;      given, taken from class object.
;   
;    dirs : in, optional, type=strarr
;   
;       The data directories to process. Or just the timestamps.
;   
;    div : 
;
;
;
;    escan :  in, optional, type=integer
;   
;       Scan number.
;   
;    extraclip : in, optional, type="intarr(4)", default="[0,0,0,0]"
; 
;      A margin (in pixels) to apply to the FOV edges of the reference
;      channel when calculating the largest common FOV. The order is
;      [left, right, top, bottom]
;
;    fill_fov : in, optional, type=boolean
; 
;      Fill the illuminated field of view with momfbd subfields.
; 
;    global_keywords : in, optional, type=strarr
;   
;      Any global keywords that you want to add to the momfbd config file.
;   
;    margin : in, optional, type=integer, default=5
; 
;      A margin (in pixels) to disregard from the FOV edges when
;      constructing the grid of MOMFBD subfields. A negative margin
;      adds to the FOV instead, which can be useful when /fill_fov is
;      used (in which case the default is -numpoints/2). 
;
;    modes : in, optional, type=string, default="'2-45,50,52-55,65,66'"
;   
;      The modes to include in the expansions of the wavefront phases.
;   
;    momfbddir :  in, optional, type=string, default='momfbd'
;   
;       Top directory of output tree.
; 
;    nfac : in, optional, type=float
;
;       Noise factor.
;   
;    nmodes : in, optional, type=integer, default=51
;
;      If keyword modes is not given, use the Nmodes most significant
;      KL modes.
;   
;    no_descatter : in, optional, type=boolean
;   
;       Set this if your data is from a near-IR (777 or 854 nm) line
;       and you do not want to do backscatter corrections.
;   
;    no_fitsheaders : in, optional, type=boolean
;   
;       Don't write fitsheader files.
;   
;    no_narrowband : in, optional, type=boolean
;   
;       Set this to include only wideband data, including diversity
;       data if any.
;   
;    no_pd : in, optional, type=boolean
;   
;       Set this to exclude phase diversity data and processing. 
;   
;    nremove : 
;   
;   
;    numpoints : in, optional, type=integer, default=88
;   
;      The size of MOMFBD subfields.
;   
;    old_wb_states : in, optional, type=boolean
;   
;      Generate extra WB objects, one per NB state, using the old
;      mechanism where they are explicitly specified in the config
;      file.
;   
;    oldgains :
;
;
;   
;    pref : 
;   
;   
; 
;    skip : 
;   
;   
;   
;    state : 
;   
;   
;   
;    unpol : 
;   
;   
;   
;    wb_states : in, optional, type=boolean
;   
;      Generate extra WB objects, one per NB state, using the TRACE
;      config file keyword.
;   
;
;
;
;
; :History:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
;
;   2013-06-13 : JdlCR. added support for scan-dependent gains ->
;                using keyword "/newgains".
;
;   2013-06-28 : JdlCR. added NF (object) option 
; 
;   2013-08-27 : MGL. Added support for logging. Let the subprogram
;                find out its own name.
;
;   2013-12-19   PS. Work based on the link directory guess date
;                before asking adapt to changed link directory names
;                NEWGAINS is the default now (removed), use OLDGAINS
;
;   2014-01-10   PS. Remove keyword outformat, use self.filetype. to
;                not be a string.
;
;   2016-02-15 : MGL. Use red_loadbackscatter. Remove keyword descatter,
;                new keyword no_descatter.
;
;   2016-02-15 : MGL. Get just the file names from
;                red_loadbackscatter, do not read the files.
;
;   2016-04-18 : THI. Added margin keyword to allow for user-defined edge trim
;                Changed numpoints keyword to be a number rather than a string.
;
;   2016-04-21 : MGL. Added some documentation. Use n_elements, not
;                keyword_set, to find out if a keyword needs to be set
;                to a default value.
;
;   2016-06-06 : MGL. Default date from class object. Better loop
;                indices. Added dirs keyword. Added keyword mfbddir.
;
;   2016-06-08 : MGL. New keyword no_pd. Renamed keyword nf to nfac so
;                scope_varfetch does not get confused in writelog.
;
;   2016-08-23 : THI. Rename camtag to detector and channel to camera,
;                so the names match those of the corresponding SolarNet
;                keywords.
;
;   2016-09-22 : THI. Re-write to support the new data- and class-organization.
;                Added keyword refcam to allow selecting the annchor-channel.
;                Added keyword redux to add keywords/options specific to redux.
;
;   2016-10-13 : MGL & THI. Implement the nremove mechanism.
;
;   2016-10-14 : MGL. Time-varying gaintables.
;
;   2017-01-18 : MGL. Implemented phase diversity.
;
;   2017-02-10 : JdlCR. Bugfix, maxshift was not used when provided as
;                a keyword.
;
;   2017-04-07 : MGL. New keyword, Nmodes.
;
;   2017-04-10 : MGL. When redux flag is set, use geometry maps
;                instead of offset files and align_clip. Run
;                prepmomfbd_fitsheaders when done.
;
;   2017-05-14 : THI. Calculate the patches within the common FOV, using global
;                coordinates for the redux-code.
;
;   2017-06-19 : MGL. Use gradient_vogel if there is PD data.
;
;   2018-03-29 : MGL. New keyword old_wb_states. With /wb_states,
;                generate extra WB objects by use of TRACE keyword.
;
;   2018-12-18 : MGL. Make the redux keyword obsolete by
;                non-optionally setting it to 1.
;
;   2019-09-20 : MGL. Remove all code for ~keyword_set(redux). Keep
;                the keyword itself for backward compatibility.
;
;   2022-09-11 : MGL. New keyword fill_fov.
;
;   2022-09-27 : MGL. Make FOV mask. New keyword no_fitsheaders.
;
;   2023-11-24 : MGL. New keyword no_narrowband.
;
;-
pro red::prepmomfbd, cams = cams $
                     , date_obs = date_obs $
                     , dirs = dirs $
                     , div = div $
                     , escan = escan $
                     , extraclip = extraclip $
                     , fill_fov = fill_fov $
                     , global_keywords = global_keywords $
                     , margin = margin $
                     , maxshift = maxshift $
                     , modes = modes $
                     , momfbddir = momfbddir $
                     , nfac = nfac $
                     , nmodes = nmodes $
                     , no_descatter = no_descatter $
                     , no_fitsheaders = no_fitsheaders $
                     , no_narrowband = no_narrowband $
                     , no_pd = no_pd $
                     , nremove = nremove $
                     , numpoints = numpoints $                     
                     , oldgains = oldgains $
                     , pref = pref $
                     , redux = redux $
                     , refcam = refcam $
                     , skip = skip $
                     , state = state $
                     , unpol = unpol $
                     , weight = weight $
                     , wb_states = wb_states

  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)
  
  instrument = ((typename(self)).tolower())

  case n_elements(escan) of
     0 : 
     1 : if size(escan, /tname) ne 'INT' then begin
        print, "Error: 'escan' keyword should be an integer number."
        return
     endif
     else : begin
        print, "Error: 'escan' keyword should be an integer number."
        return
     end
  endcase

  ;; Cameras
  if ~keyword_set(cams) then cams = *self.cameras
  iswb = strmatch(cams,'*-W') or strmatch(cams,'*-D')
  ispd = strmatch(cams,'*-D')

  if total(ispd) eq 0 then no_pd = 1 ; If no PD camera, run in no_pd mode
  
  offset_dir = self.out_dir + '/calib/'

  if(n_elements(maxshift) eq 0) then maxshift='30'
  maxshift = strcompress(string(maxshift), /remove_all)
  
  LF = string(10b)

  ;; Get keywords
  if n_elements(momfbddir) eq 0 then begin
    if keyword_set(no_narrowband) then begin
      if keyword_set(no_pd) then begin
        momfbddir = 'mfbd_nopd' 
      endif else begin
        momfbddir = 'mfbd' 
      endelse
    endif else begin
      if keyword_set(no_pd) then begin
        momfbddir = 'momfbd_nopd' 
      endif else begin
        momfbddir = 'momfbd' 
      endelse
    endelse
  endif
  if n_elements(date_obs) eq 0 then date_obs = self.isodate

  ;; If the modes are specified, then use them. Otherwise look for Nmodes.
  if n_elements(modes) eq 0 then begin
    ;; Is the number of modes specified or do we need a default number?
    if n_elements(Nmodes) eq 0 then Nmodes = 51
    ;; KL modes in variance order:
    manymodes = red_expandrange('2-6,9,10,7,8,14,15,11-13,20,21,18,19,16,17,27,28,25,26,22-24,35,36,33,34,44,45,31,32,29,30,54,55,42,43,40,41,37-39,65,66,52,53,50,51,77,78,48,49,46,47,63,64,61,62,90,91,75,76,59,60,56-58,104,105,73,74,88,89,71,72,69,70,67,68,119,120,86,87,102,103,84,85,135,136,82,83,79-81,100,101,117,118,152,153,98,99,96,97,115,116,133,134,94,95,92,93,170,171,113,114,150,151,131,132,111,112,189,190,109,110,106-108,129,130,168,169,148,149,209,210,127,128,125,126,123,124,121,122,146,147,187,188,230,231,166,167,144,145,142,143,207,208,164,165,252,253,140,141,137-139,185,186,162,163,275,276,228,229,160,161,183,184,205,206,158,159,156,157,154,155,181,182,299,300,250,251,203,204,179,180,226,227,177,178,175,176,324,325,172-174,201,202,273,274,224,225,248,249,199,200,350,351,197,198,297,298,222,223,195,196,193,194,191,192,246,247,271,272,220,221,377,378,322,323,218,219,244,245,216,217,295,296,269,270,214,215,211-213,405,406,242,243,348,349,240,241,267,268,320,321,293,294,434,435,238,239,236,237,375,376,265,266,234,235,232,233,291,292,346,347,263,264,464,465,318,319,261,262,403,404,289,290,259,260,257,258,254-256,495,496,373,374,316,317,344,345,287,288,432,433,285,286,314,315,527,528,283,284,401,402,342,343,281,282,371,372,279,280,277,278,462,463,312,313,560,561,310,311,340,341,430,431,308,309,369,370,399,400,493,494,306,307,338,339,304,305,301-303,594,595,367,368,336,337,460,461,397,398,525,526,428,429,334,335,629,630,365,366,332,333,330,331,328,329,326,327,491,492,395,396,363,364,558,559,458,459,426,427,665,666,361,362,393,394,359,360,523,524,357,358,592,593,424,425,355,356,352-354,489,490,391,392,456,457,702,703,389,390,422,423,556,557,627,628,387,388,454,455,521,522,740,741,385,386,487,488,420,421,383,384,381,382,379,380,590,591,418,419,452,453,663,664,779,780,485,486,416,417,554,555,519,520,450,451,414,415,412,413,625,626,410,411,700')
    ;; Use the Nmodes most significant modes:
    if Nmodes le n_elements(manymodes) then modes = red_collapserange(manymodes[0:Nmodes-1], ld = '', rd = '') else stop
  endif

  if n_elements(nremove) eq 0 then nremove=0
  ;;if n_elements(nfac) eq 0 then nfac = 1.
  if n_elements(nfac) eq 1 then nfac = replicate(nfac,3)

  Ndirs = n_elements(dirs)
  if Ndirs gt 0 then begin
    if Ndirs eq 1 then dirs = [dirs] 
    for idir = 0, Ndirs-1 do begin
      if ~file_test(dirs[idir]) then begin
        if file_test(self.out_dir+'data/'+dirs[idir]) then begin
           dirs[idir] = self.out_dir+'data/'+dirs[idir]
        endif else begin
          print,'The directory ',self.out_dir+'data/'+dirs[idir], " doesn't exist."
          print,"Run a->link_data, dir=['",dirs[idir],"'] first."
          return
        endelse
      endif
    endfor                      ; idir
  endif else begin
;    if ~ptr_valid(self.data_dirs) then begin
;      print, inam+' : ERROR : undefined data_dir'
;      return
;    endif
    dirs = file_search(self.out_dir+'data/*')
;    dirs = *self.data_dirs
    Ndirs = n_elements(dirs)
  endelse

  if Ndirs eq 0 then begin
    print, inam+' : ERROR : no directories defined'
    return
  endif else begin
    if Ndirs gt 1 then dirstr = '['+ strjoin(dirs,';') + ']' $
    else dirstr = dirs[0]
  endelse

  ;; Get states from the data folder
  ;;  d_dirs = file_search(self.out_dir+'/data/*', /TEST_DIR, count = Ndirs)
  IF Ndirs EQ 0 THEN BEGIN
    print, inam + ' : ERROR -> no frames found in '+self.out_dir+'/data'
    print, inam + '   Did you run link_data?'
    return
  ENDIF

  ismos = bytarr(Ndirs)         ; Should be TRUE for directories with automatic mosaic data. 
  
  if keyword_set(no_pd) then begin
    ;; Remove PD camera if any
    indx = where(~ispd)
    cams = cams[indx]
    iswb = iswb[indx]
    ispd = ispd[indx]
    Ncams = n_elements(cams)    ; Number of cameras
  endif else begin
    ;; Establish PD camera
    indx = where(iswb and ispd)
    if max(indx) ge 0 then pdcam = indx[0] $
    else pdcam = self.pdcam 
    Ncams = n_elements(cams)    ; Number of cameras
    if pdcam ge Ncams then begin
      print, inam, ' : index of PD camera out of range: ', pdcam, ' >= ', Ncams
      return
    endif
    pdcam_name = cams[pdcam]
 
    ;; Get amount of diversity
    if file_test('calib/diversity.txt') then begin
      ;; Read calibrated amount of diversity
      spawn, 'cat calib/diversity.txt', diversity_string
    endif else begin
      if self.diversity ne '' then begin 
        ;; Use the nominal amount specified in the config file
        diversity_string = strtrim(string(float(self.diversity)*1e3,format='(f5.2)')+' mm', 2)
      endif else begin
        print, inam+' : Diversity not specified.'
        stop
      endelse
    endelse
  endelse

  ;; Establish reference camera
  if(n_elements(refcam) eq 0) then begin
    indx = where(iswb and ~ispd)
    if max(indx) ge 0 then refcam = indx[0] $
    else refcam = self.refcam
  endif
  if refcam ge Ncams then begin
    print, inam, ' : index of reference camera out of range: ', refcam, ' >= ', Ncams
    return
  endif
  refcam_name = cams[refcam]

  ;; NB: this will overwrite existing offset files !!
  self -> getalignment, align=align, cams=cams, refcam=refcam, prefilters=pref $
                        , output_dir = offset_dir $
                        , extraclip=extraclip, /overwrite $
                        , makeoffsets = 0

  ref_idx = where( align.state1.camera eq refcam_name, Nref)
  if Nref eq 0 then begin
    print, inam, ' : Failed to get alignment for refererence camera: ', refcam_name
    return
  endif

  detectors = strarr(Ncams)
  for icam = 0, Ncams-1 do detectors[icam] = self -> getdetector(cams[icam])
  ;;self -> getdetectors, dir = self.data_dir

  ;; Print cams
  print, inam + ' : cameras found:'
  for icam = 0, Ncams-1 do begin
    outstr = '    ' + cams[icam] + ' ' + detectors[icam] + ' '
    if iswb[icam] then outstr += 'WB ' else outstr += 'NB '
    if ispd[icam] then outstr += 'PD '
    print, outstr
  endfor                        ; icam

  ;; Use a narrowband camera when searching for files, so we are sure
  ;; to get the states information.
  pos = where(~iswb and ~ispd)
  searchcam = cams[pos[0]]
  searchdet = detectors[pos[0]]

  if n_elements(numpoints) eq 0 then begin
    ;; About the same subfield size in arcsec as CRISP:
    numpoints = strtrim(round(88*0.0590/self.image_scale/2)*2, 2)
  endif else begin
    ;; Convert strings, just to avoid breaking existing codes.
    if( size(numpoints, /type) eq 7 ) then numpoints = fix(numpoints) 
  endelse

  if n_elements(margin) eq 0 then margin = 5

  if keyword_set(fill_fov) then begin
    ;; Use masks made from gaintable to find the part of the FOV that
    ;; has light. To be used for CRISP with new cameras.
    
    ;; Read gains for all cameras, remap them to the WB orientation
    ;; and position, AND them gt 0, to make a mask.

    
    maps = fltarr(3, 3, Ncams)
    for icam = 0, Ncams-1 do begin
      gfiles = file_search('gaintables/'+(*self.detectors)[icam]+'_*.fits')
      g = red_readdata(gfiles[0])
      indx = where(align.state2.camera eq cams[icam])
      maps[*, *, icam] = align[indx[0]].map
      if icam eq 0 then begin
        dims = size(g, /dim)
        masks = fltarr([dims, Ncams])
        mask = rdx_img_transform(invert(maps[*, *, 0]), g, /preserve) gt 0
      endif else begin
        mask AND= rdx_img_transform(invert(maps[*, *, icam]), g, /preserve) gt 0
      endelse
      masks[*, *, icam] = g gt 0 ; Masks for individual cameras
    endfor                       ; icam

    ;; Close small holes in the mask
    ste=[[0,1,0],[1,1,1],[0,1,0]]
    Nclose = 10
    for iclose = 0, Nclose do mask = dilate(mask, ste)     
    for iclose = 0, Nclose do mask = erode(mask, ste)     


    ;; Save this mask later in the cfg directories, it will be used when making cubes
    ;;  dims = size(mask, /dim)    
    fov_mask = mask
    fov_indices = Where(fov_mask EQ 1)
    fov_boundaryPts = Find_Boundary(fov_indices, XSize=dims[0], YSize=dims[1])
    fov_roiObj = Obj_New('IDLanROI', fov_boundaryPts[0,*], fov_boundaryPts[1,*])

    ;; Apply margin keyword
    cmargin = -numpoints/2
    if n_elements(cmargin) ne 0 then begin
      if cmargin gt 0 then begin
        ;; Make mask smaller
        for ierode = 0, cmargin-1 do mask = erode(mask, ste)
      endif else begin
        ;; Make mask larger
        for idilate = 0, abs(cmargin)-1 do mask = dilate(mask, ste)
      endelse 
    endif

    ;; Zero MAXSHIFT outermost rows and columns in mask, in all cameras.
    tmp = bytarr(dims[0]-2*maxshift, dims[1]-2*maxshift) + 1
    tmp = red_centerpic(tmp, xs = dims[0], ys = dims[1], z = 0)
    for icam = 0, Ncams-1 do mask AND= rdx_img_transform(invert(maps[*, *, icam]), tmp, /preserve)
;    mask AND= tmp
    
    ;; Fill the masked FOV with subfield coordinates using sim_xy
    
    ;; Define ROI based on mask
    indices = Where(mask EQ 1)
    boundaryPts = Find_Boundary(indices, XSize=dims[0], YSize=dims[1])

    ;; Apply margin
    boundaryPts[0, *] = boundaryPts[0, *] >margin <(dims[0]-margin-1)    
    boundaryPts[1, *] = boundaryPts[1, *] >margin <(dims[1]-margin-1)    

    roiObj = Obj_New('IDLanROI', boundaryPts[0,*], boundaryPts[1,*])

    ;; sim_x and sim_y based on mask ROI and margin
    sim_roi = round([roiobj.roi_xrange, roiobj.roi_yrange])
    ;;sim_roi[[0,2]] += margin    ; shrink the common FOV by margin along all edges.
    ;;sim_roi[[1,3]] -= margin
    sim_x = rdx_segment( sim_roi[0], sim_roi[1], numpoints, /momfbd )
    sim_y = rdx_segment( sim_roi[2], sim_roi[3], numpoints, /momfbd )

    scrollwindow, xs = 2000, ys = 2000
    cgplot, (*roiobj.data)[0,*], (*roiobj.data)[1,*], psym=3, /aspect $
            , xrange = [0, dims[0]], yrange = [0, dims[1]]
    cgplot, (*fov_roiobj.data)[0,*], (*fov_roiobj.data)[1,*], psym=3, /over, color = 'red'      

    
    ;; Put coordinate pairs that are inside the mask in sim_xy.
    for ix = 0, n_elements(sim_x)-1 do begin    
      for iy = 0, n_elements(sim_y)-1 do begin    

        if roiObj -> ContainsPoints(sim_x[ix], sim_y[iy]) then begin

          for icam = 0, Ncams-1 do begin
            ;; Are the coordinates within the NB frame?
            coords = rdx_point_transform( maps[*, *, icam], [sim_x[ix], sim_y[iy]])
            coords_limited = coords
            coords_limited[0] = coords[0] >(numpoints/2+margin) <(dims[0]-1-numpoints/2-margin)            
            coords_limited[1] = coords[1] >(numpoints/2+margin) <(dims[1]-1-numpoints/2-margin)            
            if coords_limited[0] ne coords[0] or coords_limited[1] ne coords[1] then begin
              ;; If not, nudge them inside and transform back.
              coords = rdx_point_transform( invert(maps[*, *, icam]), coords_limited)
              print, 'Changed!', ix, iy, icam
              sim_x[ix] = coords[0]
              sim_y[iy] = coords[1]              
            endif
          endfor                ; icam
          
          red_append, sim_xy, [sim_x[ix], sim_y[iy]]
          cgplot, /over, sim_x[ix], sim_y[iy], psym = 16
          cgplot, /over, sim_x[ix]+numpoints/2*[1, 1, -1, -1, 1], sim_y[iy]+numpoints/2*[1, -1, -1, 1, 1]
        endif
        
      endfor                    ; iy
    endfor                      ; ix
 
    sim_xy_string = strjoin(strtrim(sim_xy,2), ',')
    
  endif else begin

    ref_clip = align[0].clip
    sim_roi = ref_clip
    sim_roi[[0,2]] += margin    ; shrink the common FOV by margin along all edges.
    sim_roi[[1,3]] -= margin
    if sim_roi[0] gt sim_roi[1] || sim_roi[2] gt sim_roi[3] then begin
      print, inam + ' : Error: The region of interest looks weird. sim_roi = [' + strjoin(strtrim(sim_roi,2),',') + ']'
      print, inam + '                                               margin = ' + strtrim(margin,2)
      return
    endif
    sim_x = rdx_segment( sim_roi[0], sim_roi[1], numpoints, /momfbd )
    sim_y = rdx_segment( sim_roi[2], sim_roi[3], numpoints, /momfbd )
    sim_x_string = strjoin(strtrim(sim_x,2), ',')
    sim_y_string = strjoin(strtrim(sim_y,2), ',')

  endelse

  
  for idir=0L, Ndirs-1 do begin
    
    dir = dirs[idir]+'/'
    folder_tag = file_basename(dir)
    
    if file_test(dir + refcam_name + '_nostate/',/directory) then subdir = refcam_name + '_nostate/'
    
    print, inam + ' : Search for reference files in ' + dir
    self->selectfiles, cam=refcam_name, dirs=dir, prefilter=pref, subdir=subdir, $
                       files=ref_files, states=ref_states, /force, scan = escan

    if n_elements(ref_states) eq 0 then begin
      print, inam, ' : Failed to find files/states for the reference channel in ', dir
      return
    endif
    
    ref_img_dir = file_dirname(file_expand_path(ref_states[0].filename),/mark)
    ref_caminfo = red_camerainfo(detectors[refcam])
    ref_head = red_readhead(ref_states[0].filename)

    ;; Are these data automatic mosaics?
    ismos[idir] = strmatch(fxpar(ref_head, 'FILENAME'), '*_mos*')

    ;; Pixel size and binning
    pixelsize = ref_caminfo.pixelsize
    nbin = fxpar(ref_head, 'NBIN*', count = Nnbin)
    if Nnbin ne 0 then begin
      pixelsize *= nbin[0]      ; Assume binning same along both axes
    endif
    pixelsize = strtrim(pixelsize, 2)

    ;; unique prefilters
    upref = ref_states[uniq(ref_states.prefilter, sort(ref_states.prefilter))].prefilter
    Nprefs = n_elements(upref)

    ;; unique scan numbers
    uscan = ref_states[uniq(ref_states.scannumber, sort(ref_states.scannumber))].scannumber
    Nscans = n_elements(uscan)
    
    ;; base output location
    cfg_base_dir = self.out_dir + PATH_SEP() + momfbddir + PATH_SEP() + folder_tag


    if ~keyword_set(no_pd) then begin
      print, inam + ' : Search for PD files in ' + dir
      if file_test(dir + pdcam_name + '_nostate/',/directory) then subdir = pdcam_name + '_nostate/'
      self->selectfiles, cam=pdcam_name, dirs=dir, prefilter=pref, subdir=subdir, $
                         files=pd_files, states=pd_states, /force, scan = escan

      if n_elements(pd_states) eq 0 then begin
        print, inam, ' : Failed to find files/states for the pd channel in ', dir
        stop
      endif
      
      pd_img_dir = file_dirname(file_expand_path(pd_states[0].filename),/mark)
;      pd_caminfo = red_camerainfo(detectors[pdcam])

    endif


    for iscan=0L, Nscans-1 do begin

      red_progressbar, iscan, Nscans, 'Config info for WB', /predict

      if n_elements(escan) ne 0 then scannumber = escan else scannumber = uscan[iscan]

      scanstring = string(scannumber,format='(I05)')

      for ipref=0L, Nprefs-1 do begin
        
        self->selectfiles, prefilter=upref[ipref], scan=scannumber, $
                           files=ref_files, states=ref_states, selected=ref_sel
        
        if ref_states[ref_sel[0]].nframes eq 1 then begin
          if nremove lt n_elements(ref_sel) then ref_sel = ref_sel[nremove:*] else continue
        endif

        if( max(ref_sel) lt 0 ) then continue
        
        filename = file_basename(ref_states[ref_sel[0]].filename)
        pos = STREGEX(filename, '[0-9]{7}', length=len)
        ref_fn_template = strmid(filename, 0, pos) + '%07d' + strmid(filename, pos+len)
        
        self -> get_calib, ref_states[ref_sel[0]] $
                           , gainname = ref_gainname, darkname = ref_darkname, status = status
        if( status lt 0 ) then continue

        if ~keyword_set(no_pd) then begin

          filename = file_basename(pd_states[ref_sel[0]].filename)
          pos = STREGEX(filename, '[0-9]{7}', length=len)
          pd_fn_template = strmid(filename, 0, pos) + '%07d' + strmid(filename, pos+len)
          
          self -> get_calib, pd_states[ref_sel[0]] $
                             , gainname = pd_gainname, darkname = pd_darkname, status = status
          if( status lt 0 ) then continue

        endif

        cfg_dir = cfg_base_dir + '/'+upref[ipref]+'/cfg/'
        rdir = cfg_dir + 'results/'
        ddir = cfg_dir + 'data/'
        cfg = { dir:cfg_dir, $
                file:cfg_dir+'momfbd_reduc_'+upref[ipref]+'_'+scanstring+'.cfg', $
                globals:'', $
                objects:'', $
                framenumbers:ptr_new(ref_states[ref_sel].framenumber) $
              }
        
        ;; Global keywords
        cfg.globals += 'DATE_OBS=' + date_obs + LF
        cfg.globals += 'PROG_DATA_DIR=./data/' + LF
        cfg.globals += 'NEW_CONSTRAINTS' + LF
        cfg.globals += 'FAST_QR' + LF
        cfg.globals += 'FPMETHOD=horint' + LF
        cfg.globals += 'BASIS=Karhunen-Loeve' + LF
        cfg.globals += 'GETSTEP=getstep_conjugate_gradient' + LF
        if keyword_set(no_pd) then begin
          cfg.globals += 'GRADIENT=gradient_diff' + LF
        endif else begin
          cfg.globals += 'GRADIENT=gradient_Vogel' + LF
        endelse
        cfg.globals += 'MODES=' + modes + LF
        cfg.globals += 'TELESCOPE_D=0.97' + LF
        cfg.globals += 'MAX_LOCAL_SHIFT='+string(maxshift,format='(I0)') + LF
        cfg.globals += 'NUM_POINTS=' + strtrim(numpoints,2) + LF
        cfg.globals += 'ARCSECPERPIX=' + self.image_scale + LF
        cfg.globals += 'PIXELSIZE=' + pixelsize + LF
        cfg.globals += 'FILE_TYPE=' + self.filetype + LF
        case self.filetype of
          'ANA' : begin
            cfg.globals += 'DATA_TYPE=FLOAT' + LF
          end
          'MOMFBD' : begin
            cfg.globals += 'GET_PSF' + LF
            cfg.globals += 'GET_PSF_AVG' + LF
          end
        endcase
        if keyword_set(fill_fov) then begin
          cfg.globals += 'SIM_XY=' + sim_xy_string + LF
        endif else begin
          cfg.globals += 'SIM_X=' + sim_x_string + LF
          cfg.globals += 'SIM_Y=' + sim_y_string + LF
        endelse
        if keyword_set(wb_states) then cfg.globals += 'TRACE' + LF

        ;; External keywords?
        if(keyword_set(global_keywords)) then begin
          nk = n_elements(global_keywords)
          for ki=0L, nk-1 do cfg.globals += global_keywords[ki] + LF
        endif
        
        ;; Reference object
        cfg.objects += 'object{' + LF
        cfg.objects += '    WAVELENGTH=' + strtrim(ref_states[ref_sel[0]].pf_wavelength,2) + LF
        cfg.objects += '    OUTPUT_FILE=results/' + detectors[refcam] $
                       + '_' + date_obs+'T'+folder_tag $
                       + '_' + scanstring + '_' + upref[ipref] + LF
        if(n_elements(weight) gt 0 ) then $
           cfg.objects += '    WEIGHT=' + strtrim(weight[0],2) + LF
        cfg.objects += '    channel{' + LF
        cfg.objects += '        IMAGE_DATA_DIR=' + ref_img_dir + LF
        cfg.objects += '        FILENAME_TEMPLATE=' + ref_fn_template + LF
        cfg.objects += '        GAIN_FILE=' + ref_gainname + LF
        cfg.objects += '        DARK_TEMPLATE=' + ref_darkname + LF
        cfg.objects += '        DARK_NUM=0000001' + LF
        cfg.objects += '        ALIGN_MAP='+strjoin(strtrim(reform(align[0].map, 9), 2), ',') + LF

        if( upref[ipref] EQ '8542' OR upref[ipref] EQ '7772' ) AND $
           ~keyword_set(no_descatter) AND self.dodescatter then begin
          self -> loadbackscatter, detectors[refcam], upref[ipref] $
                                   , bgfile = bgf, bpfile = psff
          if(file_test(psff) AND file_test(bgf)) then begin
            cfg.objects += '        PSF=' + psff + LF
            cfg.objects += '        BACK_GAIN=' + bgf + LF
          endif else begin
            print, inam, ' : No backscatter files found for prefilter: ', upref[ipref]
          endelse
        endif
        if(n_elements(nfac) gt 0) then $
           cfg.objects += '        NF=' + red_stri(nfac[0]) + LF
        cfg.objects += '        INCOMPLETE' + LF
        cfg.objects += '    }' + LF

        if ~keyword_set(no_pd) then begin ; PD channel

          align_idx = where( align.state2.camera eq cams[pdcam], count)
          if count eq 0 then begin
            ;;print, inam, ' : Failed to get ANY alignment for camera/state ', cams[icam] + ':' + thisstate
            ;;stop
            continue
          endif
          
          if count gt 1 then begin
            pd_map = median(align[align_idx].map, dim = 3)
          endif else begin
            pd_map = align[align_idx].map
          endelse

          if n_elements(align_idx) gt 1 then align_idx = align_idx[0] ; just pick the first one for now
          state_align = align[align_idx]

          cfg.objects += '    channel{' + LF
          cfg.objects += '        IMAGE_DATA_DIR=' + pd_img_dir + LF
          cfg.objects += '        FILENAME_TEMPLATE=' + pd_fn_template + LF
          cfg.objects += '        GAIN_FILE=' + pd_gainname + LF
          cfg.objects += '        DARK_TEMPLATE=' + pd_darkname + LF
          cfg.objects += '        DARK_NUM=0000001' + LF
          cfg.objects += '        ALIGN_MAP='+strjoin(strtrim(reform(pd_map, 9), 2), ',') + LF
          if( upref[ipref] EQ '8542' OR upref[ipref] EQ '7772' ) AND $
             ~keyword_set(no_descatter) AND self.dodescatter then begin
            self -> loadbackscatter, detectors[refcam], upref[ipref] $
                                     , bgfile = bgf, bpfile = psff
            if(file_test(psff) AND file_test(bgf)) then begin
              cfg.objects += '        PSF=' + psff + LF
              cfg.objects += '        BACK_GAIN=' + bgf + LF
            endif else begin
              print, inam, ' : No backscatter files found for prefilter: ', upref[ipref]
            endelse
          endif
          if(n_elements(nfac) gt 0) then $
             cfg.objects += '        NF=' + red_stri(nfac[0]) + LF
          cfg.objects += '        INCOMPLETE' + LF
          cfg.objects += '        DIVERSITY=' + diversity_string + LF
          cfg.objects += '    }' + LF
        endif                   ; PD channel
        cfg.objects += '}' + LF

        red_append, cfg_list, cfg

      endfor                    ; ipref

    endfor                      ; iscan 

    if keyword_set(no_narrowband) then continue ; Don't add NB objects
    
    for icam=0L, Ncams-1 do begin

      ;;if icam ne refcam then begin
      if ~iswb[icam] then begin ; This excludes also a WB PD camera
        
        ;; Get a list of all states for this camera
        self->selectfiles, cam=cams[icam], dirs=dir, files=files, $
                           states=states, /force, scan = escan ; , nremove=remove

        for iscan=0L, Nscans-1 do begin
          
          red_progressbar, iscan, Nscans, 'Config info for NB '+cams[icam], /predict
      
          if n_elements(escan) ne 0 then scannumber = escan else scannumber = uscan[iscan]

          scanstring = string(scannumber,format='(I05)')

          for ipref=0L, Nprefs-1 do begin        
            ;; select a subset of the states which matches prefilter & scan number.
            self->selectfiles, cam=cams[icam], dirs=dir, prefilter=upref[ipref], scan=scannumber, $
                               files=files, states=states, selected=sel ; , nremove=remove

            if( max(sel) lt 0 ) then continue
            
            cfg_dir = cfg_base_dir + '/'+upref[ipref]+'/cfg/'
            cfg_file = cfg_dir+'momfbd_reduc_'+upref[ipref]+'_'+scanstring+'.cfg'
            cfg_idx = where( cfg_list.file eq cfg_file )
            
            if( max(cfg_idx) lt 0 ) then continue
            
            state_list = states[sel]
            ustates = state_list[uniq(state_list.fullstate, sort(state_list.fullstate))]
            Nstates = n_elements(ustates)
            
            ;; Loop over states and add object to cfg_list
            for istate=0L, Nstates-1 do begin

              if instrument eq 'chromis' then begin
                thisstate = ustates[istate].fpi_state
                state_idx = where(state_list.fpi_state eq thisstate)
              endif else begin
                ;; Or does this work for chromis as well? (Se also
                ;; align further down!)
                thisstate = ustates[istate].fullstate
                state_idx = where(state_list.fullstate eq thisstate)
              endelse
              if max(state_idx) lt 0 then continue
              
;              self -> get_calib, state_list[state_idx[0]] $
;                                 , gainname = gainname, darkname = darkname, status = status
;              if( status lt 0 ) then stop ;continue
              
;              darkname = self -> filenames('dark', state_list[state_idx[0]], /no_fits)
              darkname = self -> filenames('dark', state_list[state_idx[0]])
              if(~keyword_set(unpol)) then begin
                if(keyword_set(oldgains)) then begin
;                  search = self.out_dir+'/gaintables/'+self.camttag + '.' + ustat1[ii] + '*.gain'
                  stop
                  ;; Not implemented in chromis::filenames yet.
;                  gainname = a -> filenames('oldgain', state_list[state_idx[0]], /no_fits)
                  gainname = self -> filenames('oldgain', state_list[state_idx[0]])
                endif else begin
                  ;;search =
                  ;;self.out_dir+'/gaintables/'+folder_tag+'/'+self.camttag
                  ;;+ '.' + istate+'.gain'
                  gainname = self -> filenames('scangain', state_list[state_idx[0]] $
                                               , timestamp = stregex(state_list[state_idx[0]].filename $
                                                                     ,'[0-9][0-9]:[0-9][0-9]:[0-9][0-9]' $
                                                                     ,/extr))
;                                                                     ,/extr), /no_fits)
                endelse
              endif else begin
;                  gainname = self -> filenames('cavityfree_gain', state_list[state_idx[0]], /no_fits)
                  gainname = self -> filenames('cavityfree_gain', state_list[state_idx[0]])
;
;                 search = self.out_dir+'/gaintables/'+self.camttag + $
;                          '.' + strmid(ustat1[ii], idx[0], $
;                                       idx[nidx-1])+ '*unpol.gain'
              endelse
  
              if state_list[state_idx[0]].nframes eq 1 then begin
                if nremove ge n_elements(state_idx) then continue
                if nremove ne 0 then begin
                  *(cfg_list[cfg_idx].framenumbers) = red_strip(*(cfg_list[cfg_idx].framenumbers) $
                                                                , state_list[state_idx[0:nremove-1]].framenumber)
                  state_idx = state_idx[nremove:*] 
                endif
              endif
              red_append, *(cfg_list[cfg_idx].framenumbers), state_list[state_idx].framenumber
              
              img_dir = file_dirname(file_expand_path(state_list[state_idx[0]].filename),/mark)
              filename = file_basename(state_list[state_idx[0]].filename)
              pos = STREGEX(filename, '[0-9]{7}', length=len)
              fn_template = strmid(filename, 0, pos) + '%07d' + strmid(filename, pos+len)
              
              ;; Use median of align info for camera+prefilter combo.
              ;; (There does not seem to be any consistency in
              ;; parameters that motivate wavelength interpolation.)
              align_idx = where( align.state2.camera eq cams[icam] and $
                                 align.state2.prefilter eq ustates[istate].prefilter, Nalign)
              case Nalign of
                 0 : begin
                    print, inam, ' : Failed to get ANY alignment for camera/state ' $
                           , cams[icam] + ':' + thisstate
                    stop
                 end
                 1    : align_map =        align[align_idx].map
                 else : align_map = median(align[align_idx].map,dim=3)
              endcase
               
              ;; Create cfg object
              cfg_list[cfg_idx].objects += 'object{' + LF
              cfg_list[cfg_idx].objects += '    WAVELENGTH=' + strtrim(ustates[istate].pf_wavelength,2) + LF
              cfg_list[cfg_idx].objects += '    OUTPUT_FILE=results/' + detectors[icam] $
                                           + '_' + date_obs+'T'+folder_tag $
                                           + '_' + scanstring + '_'+ustates[istate].fullstate + LF
              if(n_elements(weight) gt 1) then $
                 cfg_list[cfg_idx].objects += '    WEIGHT=' + strtrim(weight[1],2) + LF
              cfg_list[cfg_idx].objects += '    channel{' + LF
              cfg_list[cfg_idx].objects += '        IMAGE_DATA_DIR=' + img_dir + LF
              cfg_list[cfg_idx].objects += '        FILENAME_TEMPLATE=' + fn_template + LF ; ******
              cfg_list[cfg_idx].objects += '        GAIN_FILE=' + gainname + LF            ; ******
              cfg_list[cfg_idx].objects += '        DARK_TEMPLATE=' + darkname + LF
              cfg_list[cfg_idx].objects += '        DARK_NUM=0000001' + LF
              cfg_list[cfg_idx].objects += '        ALIGN_MAP='+strjoin(strtrim(reform(align_map, 9), 2), ',') + LF
              if( ustates[istate].prefilter EQ '8542' OR ustates[istate].prefilter EQ '7772' ) AND $
                 ~keyword_set(no_descatter) AND self.dodescatter then begin
                self -> loadbackscatter, detectors[icam], ustates[istate].prefilter, bgfile = bgf, bpfile = psff
                if(file_test(psff) AND file_test(bgf)) then begin
                  cfg_list[cfg_idx].objects += '        PSF=' + psff + LF
                  cfg_list[cfg_idx].objects += '        BACK_GAIN=' + bgf + LF
                endif else begin
                  print, inam, ' : No backscatter files found for prefilter: ' $
                         , ustates[istate].prefilter
                endelse
              endif
              if(n_elements(nfac) gt 1) then $
                 cfg_list[cfg_idx].objects += '        NF=' + red_stri(nfac[1]) + LF
              if max(nremove) gt 0 then $
                 cfg_list[cfg_idx].objects += '        DISCARD=' $
                                              + strjoin(strtrim(nremove,2),',') + LF
              cfg_list[cfg_idx].objects += '        INCOMPLETE' + LF
              cfg_list[cfg_idx].objects += '    }' + LF
              cfg_list[cfg_idx].objects += '}' + LF
              
              if keyword_set(old_wb_states) then begin

                ;; select WB files with the same framenumbers
                self->selectfiles, prefilter=upref[ipref], scan=scannumber, $
                                   files=ref_files, states=ref_states, selected=ref_sel $
                                   , framenumbers = state_list[state_idx].framenumber
                
                if( max(ref_sel) lt 0 ) then continue
                
                self -> get_calib, ref_states[ref_sel[0]] $
                                   , gainname = gainname, darkname = darkname, status = status
                if( status lt 0 ) then continue
                
                ;; Use the real file name, not the link, because the
                ;; links (even in the dir w/o _nostate) do not have
                ;; the NB state info.
                fullname = ref_states[ref_sel[0]].filename
                if file_test(fullname,/symlink) gt 0 then fullname = file_readlink(fullname)
                
                img_dir = file_dirname(fullname,/mark)
                filename = file_basename(fullname)
                pos = STREGEX(filename, '[0-9]{7}', length=len)
                fn_template = strmid(filename, 0, pos) + '%07d' + strmid(filename, pos+len)
                
                cfg_list[cfg_idx].objects += 'object{' + LF
                cfg_list[cfg_idx].objects += '    WAVELENGTH=' $
                                             + strtrim(ustates[istate].pf_wavelength,2) + LF
                if(n_elements(weight) gt 2) then $
                   cfg.objects += '    WEIGHT=' + strtrim(weight[3],2) + LF $
                else cfg_list[cfg_idx].objects += '    WEIGHT=0.00' + LF
                cfg_list[cfg_idx].objects += '    OUTPUT_FILE=results/'+detectors[refcam] + '_' $
                                             + date_obs+'T'+folder_tag $
                                             + '_' + scanstring + '_'+ustates[istate].fullstate + LF
                cfg_list[cfg_idx].objects += '    channel{' + LF
                cfg_list[cfg_idx].objects += '        IMAGE_DATA_DIR=' + img_dir + LF
                cfg_list[cfg_idx].objects += '        FILENAME_TEMPLATE=' + fn_template + LF
                cfg_list[cfg_idx].objects += '        GAIN_FILE=' + gainname + LF
                cfg_list[cfg_idx].objects += '        DARK_TEMPLATE=' + darkname + LF
                cfg_list[cfg_idx].objects += '        DARK_NUM=0000001' + LF
                cfg_list[cfg_idx].objects += '        ALIGN_MAP='+strjoin(strtrim(reform(align[0].map, 9), 2), ',') + LF
                if( ustates[istate].prefilter EQ '8542' OR $
                    ustates[istate].prefilter EQ '7772' ) AND $
                   ~keyword_set(no_descatter) AND self.dodescatter then begin
                  self -> loadbackscatter, detectors[refcam], ustates[istate].prefilter $
                                           , bgfile = bgf, bpfile = psff
                  if(file_test(psff) AND file_test(bgf)) then begin
                    cfg_list[cfg_idx].objects += '        PSF=' + psff + LF
                    cfg_list[cfg_idx].objects += '        BACK_GAIN=' + bgf + LF
                  endif else begin
                    print, inam, ' : No backscatter files found for prefilter: ' $
                           , ustates[istate].prefilter
                  endelse
                endif
                if(n_elements(nfac) gt 2) then cfg_list[cfg_idx].objects $
                   += '        NF=' + red_stri(nfac[2]) + LF
                if max(nremove) gt 0 then $
                   cfg_list[cfg_idx].objects += '        DISCARD=' $
                                                + strjoin(strtrim(nremove,2),',') + LF
                cfg_list[cfg_idx].objects += '        INCOMPLETE' + LF
                cfg_list[cfg_idx].objects += '    }' + LF
                cfg_list[cfg_idx].objects += '}' + LF
              endif
              
            endfor              ; istate
            
          endfor                ; ipref
          
        endfor                  ; iscan
        
      endif                     ; NB camera?
      
    endfor                      ; icam 
    
  endfor                        ; idir 
  
  for icfg=0, n_elements(cfg_list)-1 do begin
    
    red_progressbar, icfg, n_elements(cfg_list) $
                     , 'Write config file ' + red_strreplace(cfg_list[icfg].file, self.out_dir, '')
    
    if( ~file_test(cfg_list[icfg].dir+'/data/', /directory) ) then begin
      file_mkdir, cfg_list[icfg].dir+'/data/'
      file_mkdir, cfg_list[icfg].dir+'/results/'
    endif
    
    if n_elements(fov_mask) gt 0 then writefits, cfg_list[icfg].dir+'/fov_mask.fits', fov_mask
    
;    number_str = string(cfg_list[icfg].first_file, format='(I07)') $
;                 + '-' + string(cfg_list[icfg].last_file,format='(I07)')
;        cfg_list[icfg].globals += 'IMAGE_NUMS=' + number_str +
;        LF + LF
    frmnums = *(cfg_list[icfg].framenumbers)
    frmnums = frmnums[uniq(frmnums, sort(frmnums))]
    cfg_list[icfg].globals += 'IMAGE_NUMS=' + red_collapserange(frmnums, ld='', rd='') + LF

;    print,'Writing: ', cfg_list[icfg].file
    openw, lun, cfg_list[icfg].file, /get_lun, width=2500
    printf, lun, cfg_list[icfg].objects + cfg_list[icfg].globals
    free_lun, lun
    
    ptr_free, cfg_list[icfg].framenumbers


  endfor                        ; icfg


;  ;; Make header-only fits files to be read post-momfbd.
;  if ~keyword_set(no_fitsheaders) && ~ismos[idir] then $
;     self -> prepmomfbd_fitsheaders, dirs=dirs, momfbddir=momfbddir, pref = pref, scanno = escan
  
  ;; Call prepmomfbd_mosaic on dirs that are automatic mosaics.
  for idir = 0, Ndirs-1 do begin
    if ismos[idir] then begin
      ;; Replicate config files for each mosaic tile
      self -> prepmomfbd_mosaic, dirs=dirs[idir], momfbddir=momfbddir, pref = pref, no_fitsheaders = no_fitsheaders
    endif else begin
      ;; Make header-only fits files to be read post-momfbd.
      if ~keyword_set(no_fitsheaders) then $
         self -> prepmomfbd_fitsheaders, dirs=dirs, momfbddir = momfbddir, pref = pref $
                                         , scanno = escan, no_narrowband = no_narrowband
    endelse
  endfor                        ; idir

  
  return
  
  
  
;;;;     
;;;;     for idir = 0L, Ndirs - 1 do begin
;;;;   
;;;;        data_dir = dirs[idir]
;;;;        folder_tag = file_basename(data_dir)
;;;;        searchcamdir = data_dir + '/' + searchcam + '/'
;;;;   
;;;;   ;     search = self.out_dir+'/data/'+folder_tag+'/'+self.camt
;;;;        search = searchcamdir + '*' + searchdet + '*'
;;;;        files = file_search(search, count = Nfiles) 
;;;;        
;;;;        IF Nfiles EQ 0 THEN BEGIN
;;;;            print, inam + ' : ERROR -> no frames found : '+search
;;;;            print, inam + '   Did you run link_data?'
;;;;            return
;;;;        ENDIF 
;;;;   
;;;;        files = red_sortfiles(temporary(files))
;;;;        
;;;;        ;; Get image unique states
;;;;        self -> extractstates, files, states
;;;;   ;     stat = red_getstates(files, /LINKS)
;;;;        
;;;;        ;; skip leading frames?
;;;;        IF nremove GT 0 THEN red_flagtuning, stat, nremove
;;;;   
;;;;        ;; Get unique prefilters
;;;;        upref = states[uniq(states.prefilter, sort(states.prefilter))].prefilter
;;;;        Nprefs = n_elements(upref)
;;;;   
;;;;        ;; Get scan numbers
;;;;        uscan = states[uniq(states.scannumber, sort(states.scannumber))].scannumber
;;;;        Nscans = n_elements(uscan)
;;;;   
;;;;        ;;states = stat.hscan+'.'+stat.state
;;;;        ;;pos = uniq(states, sort(states))
;;;;        ;;ustat = stat.state[pos]
;;;;        ;;ustatp = stat.pref[pos]
;;;;        ;;                           ;ustats = stat.scan[pos]
;;;;   
;;;;        ;;ntt = n_elements(ustat)
;;;;        ;;hscans = stat.hscan[pos]
;;;;   
;;;;        ;; Create a reduc file per prefilter and scan number?
;;;;        outdir0 = self.out_dir + '/' + momfbddir + '/' + folder_tag
;;;;   
;;;;        ;; Choose offset state
;;;;        for iscan = 0L, Nscans-1 do begin
;;;;   
;;;;           IF n_elements(escan) NE 0 THEN IF iscan NE escan THEN CONTINUE 
;;;;   
;;;;           scannumber = uscan[iscan]
;;;;           scanstring = string(scannumber,format='(I05)')
;;;;   
;;;;           for ipref = 0L, Nprefs-1 do begin
;;;;   
;;;;              if(keyword_set(pref)) then begin
;;;;                 if(upref[ipref] NE pref) then begin
;;;;                    print, inam + ' : Skipping prefilter -> ' + upref[ipref]
;;;;                    continue
;;;;                 endif
;;;;              endif
;;;;   
;;;;              ;; Load align clips
;;;;   ;            clipfile = self.out_dir + '/calib/align_clips.'+upref[ipref]+'.sav'
;;;;   ;            IF(~file_test(clipfile)) THEN BEGIN
;;;;   ;               print, inam + ' : ERROR -> align_clip file not found'
;;;;   ;               print, inam + ' : -> you must run red::getalignclips first!'
;;;;   ;               continue
;;;;   ;            endif
;;;;   ;            restore, clipfile
;;;;   ;            wclip = acl[0]
;;;;   ;            tclip = acl[1]
;;;;   ;            rclip = acl[2]
;;;;   
;;;;              wclip = align[0].clip
;;;;              xsz = abs(wclip[0,0]-wclip[1,0])+1
;;;;              ysz = abs(wclip[2,0]-wclip[3,0])+1
;;;;              this_margin = max([0, min([xsz/3, ysz/3, margin])])  ; prevent silly margin values
;;;;              ; generate patch positions with margin
;;;;              sim_x = rdx_segment( this_margin, xsz-this_margin, numpoints, /momfbd )
;;;;              sim_y = rdx_segment( this_margin, ysz-this_margin, numpoints, /momfbd )
;;;;              sim_x_string = strjoin(strtrim(sim_x,2), ',')
;;;;              sim_y_string = strjoin(strtrim(sim_y,2), ',')
;;;;   
;;;;              lam = strmid(string(float(upref[ipref]) * 1.e-10), 2)
;;;;   
;;;;              cfg_file = 'momfbd.reduc.'+upref[ipref]+'.'+scanstring+'.cfg'
;;;;              outdir = outdir0 + '/'+upref[ipref]+'/cfg/'
;;;;              rdir = outdir + 'results/'
;;;;              ddir = outdir + 'data/'
;;;;              if( ~file_test(rdir, /directory) ) then begin
;;;;                  file_mkdir, rdir
;;;;                  file_mkdir, ddir
;;;;              endif
;;;;              if(n_elements(lun) gt 0) then free_lun, lun
;;;;              openw, lun, outdir + cfg_file, /get_lun, width=2500
;;;;   
;;;;              ;; Image numbers
;;;;              numpos = where((states.scannumber eq uscan[iscan]) AND (states.skip eq 0B) AND (states.prefilter eq upref[ipref]), ncount)
;;;;              if( ncount eq 0 ) then continue
;;;;              n0 = min(states[numpos].framenumber)
;;;;              n1 = max(states[numpos].framenumber)
;;;;              numbers = states[numpos].framenumber
;;;;              number_range = [ min(states[numpos].framenumber), max(states[numpos].framenumber)]
;;;;              number_str = string(number_range[0],format='(I07)') + '-' + string(number_range[1],format='(I07)')
;;;;              ;nall = strjoin(strtrim(states[numpos].framenumber,2),',')
;;;;              print, inam+' : Prefilter = '+upref[ipref]+' -> scan = '+scanstring+' -> image range = ['+number_str+']'
;;;;   
;;;;              self -> get_calib, states[numpos[0]], gainname = gainname, darkname = darkname, status = status
;;;;              if status ne 0 then begin
;;;;                  print, inam+' : no dark/gain found for camera: ', refcam_name
;;;;                  continue
;;;;              endif
;;;;              
;;;;              ;; WB anchor channel
;;;;              printf, lun, 'object{'
;;;;              printf, lun, '  WAVELENGTH=' + lam
;;;;              printf, lun, '  OUTPUT_FILE=results/'+align[0].state1.detector+'.'+scanstring+'.'+upref[ipref]
;;;;              if(n_elements(weight) eq 3) then printf, lun, '  WEIGHT='+string(weight[0])
;;;;              printf, lun, '  channel{'
;;;;              printf, lun, '    IMAGE_DATA_DIR='+self.out_dir+'/data/'+folder_tag $
;;;;                      + '/' +align[0].state1.camera+'_nostate/'
;;;;              printf, lun, '    FILENAME_TEMPLATE='+align[0].state1.detector+'.'+scanstring $
;;;;                      +'.'+upref[ipref]+'.%07d'
;;;;              printf, lun, '    GAIN_FILE=' + gainname
;;;;              printf, lun, '    DARK_TEMPLATE=' + darkname
;;;;              printf, lun, '    DARK_NUM=0000001'
;;;;              printf, lun, '    ALIGN_CLIP=' + strjoin(strtrim(wclip,2),',')
;;;;              if (upref[ipref] EQ '8542' OR upref[ipref] EQ '7772' ) AND ~keyword_set(no_descatter) then begin
;;;;                 self -> loadbackscatter, align[0].state1.detector, upref[ipref], bgfile = bgf, bpfile = psff
;;;;   ;              psff = self.descatter_dir+'/'+align[0].state1.detector+'.psf.f0'
;;;;   ;              bgf = self.descatter_dir+'/'+align[0].state1.detector+'.backgain.f0'
;;;;   ;              if(file_test(psff) AND file_test(bgf)) then begin
;;;;                 printf, lun, '    PSF='+psff
;;;;                 printf, lun, '    BACK_GAIN='+bgf
;;;;   ;              endif
;;;;              endif 
;;;;   
;;;;              if(keyword_set(div)) then begin
;;;;                 printf, lun, '    DIVERSITY='+string(div[0])+' mm'
;;;;              endif
;;;;              
;;;;              if( align[0].xoffs_file ne '' && file_test(align[0].xoffs_file)) then $
;;;;                 printf, lun, '    XOFFSET='+align[0].xoffs_file
;;;;              if( align[0].yoffs_file ne '' && file_test(align[0].yoffs_file)) then $
;;;;                 printf, lun, '    YOFFSET='+align[0].yoffs_file
;;;;   
;;;;              if(n_elements(nfac) gt 0) then printf,lun,'    NF=',red_stri(nfac[0])
;;;;              printf, lun, '  }'
;;;;              printf, lun, '}'
;;;;   
;;;;              ;; Loop all wavelengths
;;;;              pos1 = where((ustatp eq upref[ipref]), count)
;;;;              if(count eq 0) then continue
;;;;              ustat1 = ustat[pos1]
;;;;   
;;;;              for ii = 0L, count - 1 do BEGIN
;;;;                 
;;;;                 ;; External states?
;;;;                 if(keyword_set(state)) then begin
;;;;                    dum = where(state eq ustat1[ii], cstate)
;;;;                    if(cstate eq 0) then continue
;;;;                    print, inam+' : found '+state+' -> scan = '+scanstring
;;;;                 endif
;;;;   
;;;;                 self -> whichoffset, ustat1[ii], xoff = xoff, yoff = yoff
;;;;   
;;;;                 ;; Trans. camera
;;;;                 istate = red_encode_scan(hscans[pos1[ii]], scan)+'.'+ustat1[ii]
;;;;   
;;;;                 ;; lc4?
;;;;                 tmp = strsplit(istate,'.', /extract)
;;;;                 ntmp = n_elements(tmp)
;;;;   
;;;;                 idx = strsplit(ustat1[ii],'.')
;;;;                 nidx = n_elements(idx)
;;;;                 iwavt = strmid(ustat1[ii], idx[0], idx[nidx-1]-1)
;;;;   
;;;;                 if(keyword_set(skip)) then begin
;;;;                    dum = where(iwavt eq skip, ccout)
;;;;                    if ccout ne 0 then begin
;;;;                       print, inam+' : skipping state -> '+ustat1[ii]
;;;;                       continue
;;;;                    endif
;;;;                 endif
;;;;   
;;;;                 printf, lun, 'object{'
;;;;                 printf, lun, '  WAVELENGTH=' + lam
;;;;                 printf, lun, '  OUTPUT_FILE=results/'+self.camttag+'.'+istate 
;;;;                 if(n_elements(weight) eq 3) then printf, lun, '  WEIGHT='+string(weight[1])
;;;;                 printf, lun, '  channel{'
;;;;                 printf, lun, '    IMAGE_DATA_DIR='+self.out_dir+'/data/'+folder_tag+ '/' +self.camt+'/'
;;;;                 printf, lun, '    FILENAME_TEMPLATE='+self.camttag+'.'+istate+'.%07d'
;;;;   
;;;;                 if(~keyword_set(unpol)) then begin
;;;;                    if(keyword_set(oldgains)) then begin
;;;;                       search = self.out_dir+'/gaintables/'+self.camttag + '.' + ustat1[ii] + '*.gain'
;;;;                    endif else begin
;;;;                       search = self.out_dir+'/gaintables/'+folder_tag+'/'+self.camttag + '.' + istate+'.gain'
;;;;                    endelse
;;;;                 endif Else begin
;;;;   
;;;;                    search = self.out_dir+'/gaintables/'+self.camttag + $
;;;;                             '.' + strmid(ustat1[ii], idx[0], $
;;;;                                          idx[nidx-1])+ '*unpol.gain'
;;;;                 endelse
;;;;                 printf, lun, '    GAIN_FILE=' + file_search(search)
;;;;                 printf, lun, '    DARK_TEMPLATE='+self.out_dir+'/darks/'+self.camttag+'.summed.0000001'
;;;;                 printf, lun, '    DARK_NUM=0000001'
;;;;                 printf, lun, '    ' + tclip
;;;;   
;;;;                 xofile = offset_dir+self.camttag+'.'+xoff
;;;;                 yofile = offset_dir+self.camttag+'.'+yoff
;;;;                 if(file_test(xofile)) then printf, lun, '    XOFFSET='+xofile
;;;;                 if(file_test(yofile)) then printf, lun, '    YOFFSET='+yofile
;;;;   
;;;;                 if (upref[ipref] EQ '8542' OR upref[ipref] EQ '7772' ) AND ~keyword_set(no_descatter) then begin
;;;;                    self -> loadbackscatter, self.camttag, upref[ipref], bgfile = bgf, bpfile = psff
;;;;   ;                 psff = self.descatter_dir+'/'+self.camttag+'.psf.f0'
;;;;   ;                 bgf = self.descatter_dir+'/'+self.camttag+'.backgain.f0'
;;;;   ;                 if(file_test(psff) AND file_test(bgf)) then begin
;;;;                    printf, lun, '    PSF='+psff
;;;;                    printf, lun, '    BACK_GAIN='+bgf
;;;;   ;              endif
;;;;                 endif 
;;;;   
;;;;                 if(keyword_set(div)) then begin
;;;;                    printf, lun, '    DIVERSITY='+string(div[1])+' mm'
;;;;                 endif
;;;;                 if(n_elements(nfac) gt 0) then printf,lun,'    NF=',red_stri(nfac[1])
;;;;              
;;;;                 printf, lun, '    INCOMPLETE'
;;;;                 printf, lun, '  }'
;;;;                 printf, lun, '}'  
;;;;   
;;;;                 ;; Reflected camera
;;;;                 printf, lun, 'object{'
;;;;                 printf, lun, '  WAVELENGTH=' + lam
;;;;                 printf, lun, '  OUTPUT_FILE=results/'+self.camrtag+'.'+istate 
;;;;                 if(n_elements(weight) eq 3) then printf, lun, '  WEIGHT='+string(weight[2])
;;;;                 printf, lun, '  channel{'
;;;;                 printf, lun, '    IMAGE_DATA_DIR='+self.out_dir+'/data/'+folder_tag+ '/' +self.camr+'/'
;;;;                 printf, lun, '    FILENAME_TEMPLATE='+self.camrtag+'.'+istate+'.%07d'
;;;;                                   ;   printf, lun, '    DIVERSITY=0.0 mm' 
;;;;                 if(~keyword_set(unpol)) then begin
;;;;                    if(keyword_set(oldgains)) then begin
;;;;                       search = self.out_dir+'/gaintables/'+self.camrtag + '.' + ustat1[ii] + '*.gain'
;;;;                    endif else begin
;;;;                       search = self.out_dir+'/gaintables/'+folder_tag+'/'+self.camrtag + '.' + istate+'.gain'
;;;;                    endelse
;;;;                 endif Else begin
;;;;                    idx = strsplit(ustat1[ii],'.')
;;;;                    nidx = n_elements(idx)
;;;;                    search = file_search(self.out_dir+'/gaintables/'+self.camrtag + $
;;;;                                         '.' + strmid(ustat1[ii], idx[0], $
;;;;                                                      idx[nidx-1])+ '*unpol.gain')
;;;;                                   ;if tmp[ntmp-1] eq 'lc4' then search = self.out_dir+'/gaintables/'+$
;;;;                                   ;                                      self.camrtag + '.' + ustat[pos[ii]] + $
;;;;                                   ;                                      '*.gain'
;;;;                 endelse
;;;;                 printf, lun, '    GAIN_FILE=' + file_search(search)
;;;;                 printf, lun, '    DARK_TEMPLATE='+self.out_dir+'/darks/'+self.camrtag+'.summed.0000001'
;;;;                 printf, lun, '    DARK_NUM=0000001'
;;;;                 printf, lun, '    ' + rclip
;;;;                 xofile = offset_dir+self.camrtag+'.'+xoff
;;;;                 yofile = offset_dir+self.camrtag+'.'+yoff
;;;;                 if(file_test(xofile)) then printf, lun, '    XOFFSET='+xofile
;;;;                 if(file_test(yofile)) then printf, lun, '    YOFFSET='+yofile
;;;;                                   ;
;;;;                 if (upref[ipref] EQ '8542' OR upref[ipref] EQ '7772' ) AND ~keyword_set(no_descatter) then begin
;;;;                    self -> loadbackscatter, self.camrtag, upref[ipref], bgfile = bgf, bpfile = psff
;;;;   ;                 psff = self.descatter_dir+'/'+self.camrtag+'.psf.f0'
;;;;   ;                 bgf = self.descatter_dir+'/'+self.camrtag+'.backgain.f0'
;;;;   ;                 if(file_test(psff) AND file_test(bgf)) then begin
;;;;                    printf, lun, '    PSF='+psff
;;;;                    printf, lun, '    BACK_GAIN='+bgf
;;;;   ;                 endif
;;;;                 endif 
;;;;   
;;;;                 if(keyword_set(div)) then begin
;;;;                    printf, lun, '    DIVERSITY='+string(div[2])+' mm'
;;;;                 endif
;;;;                 if(n_elements(nfac) gt 0) then printf,lun,'    NF=',red_stri(nfac[2])
;;;;   
;;;;                 printf, lun, '    INCOMPLETE'
;;;;                 printf, lun, '  }'
;;;;                 printf, lun, '}'
;;;;   
;;;;                 ;; WB with states (for de-warping to the anchor, only to
;;;;                 ;; remove rubbersheet when differential seeing is
;;;;                 ;; strong)
;;;;                 if(keyword_set(wb_states)) then begin
;;;;                    printf, lun, 'object{'
;;;;                    printf, lun, '  WAVELENGTH=' + lam
;;;;                    printf, lun, '  WEIGHT=0.00'
;;;;                    printf, lun, '  OUTPUT_FILE=results/'+align[0].state1.detector+'.'+istate 
;;;;                    printf, lun, '  channel{'
;;;;                    printf, lun, '    IMAGE_DATA_DIR='+self.out_dir+'/data/'+folder_tag+ '/' +align[0].state1.camera+'/'
;;;;                    printf, lun, '    FILENAME_TEMPLATE='+align[0].state1.detector+'.'+istate+'.%07d'
;;;;                                   ; printf, lun, '    DIVERSITY=0.0 mm'
;;;;                    printf, lun, '    GAIN_FILE=' + file_search(self.out_dir+'/gaintables/'+align[0].state1.detector + $
;;;;                                                                '.' + upref[ipref] + '*.gain')
;;;;                    printf, lun, '    DARK_TEMPLATE='+self.out_dir+'/darks/'+align[0].state1.detector+'.summed.0000001'
;;;;                    printf, lun, '    DARK_NUM=0000001'
;;;;                    printf, lun, '    ' + wclip
;;;;                    
;;;;                    if (upref[ipref] EQ '8542' OR upref[ipref] EQ '7772' ) AND ~keyword_set(no_descatter) then begin
;;;;                       self -> loadbackscatter, align[0].state1.detector, upref[ipref], bgfile = bgf, bpfile = psff
;;;;   ;                    psff = self.descatter_dir+'/'+align[0].state1.detector+'.psf.f0'
;;;;   ;                    bgf = self.descatter_dir+'/'+align[0].state1.detector+'.backgain.f0'
;;;;   ;                    if(file_test(psff) AND file_test(bgf)) then begin
;;;;                       printf, lun, '    PSF='+psff
;;;;                       printf, lun, '    BACK_GAIN='+bgf
;;;;   ;                    endif
;;;;                    endif 
;;;;   
;;;;                    if(keyword_set(div)) then begin
;;;;                       printf, lun, '    DIVERSITY='+string(div[0])+' mm'
;;;;                    endif
;;;;                    xofile = offset_dir+align[0].state1.detector+'.'+xoff
;;;;                    yofile = offset_dir+align[0].state1.detector+'.'+yoff
;;;;                    if(file_test(xofile)) then printf, lun, '    XOFFSET='+xofile
;;;;                    if(file_test(yofile)) then printf, lun, '    YOFFSET='+yofile
;;;;                    if(n_elements(nfac) gt 0) then printf,lun,'    NF=',red_stri(nfac[0])
;;;;   
;;;;                    printf, lun, '    INCOMPLETE'
;;;;                    printf, lun, '  }'
;;;;                    printf, lun, '}'
;;;;                 endif          
;;;;               endfor              ; ii
;;;;   
;;;;              ;; Global keywords
;;;;              printf, lun, 'PROG_DATA_DIR=./data/'
;;;;              printf, lun, 'DATE_OBS='+date_obs
;;;;              printf, lun, 'IMAGE_NUMS='+nall       ;;  n0+'-'+n1
;;;;              printf, lun, 'BASIS=Karhunen-Loeve'
;;;;              printf, lun, 'MODES='+modes
;;;;              printf, lun, 'NUM_POINTS='+strtrim(numpoints,2)
;;;;              printf, lun, 'TELESCOPE_D=0.97'
;;;;              printf, lun, 'ARCSECPERPIX='+self.image_scale
;;;;              printf, lun, 'PIXELSIZE=16.0E-6'
;;;;              printf, lun, 'GETSTEP=getstep_conjugate_gradient'
;;;;              printf, lun, 'GRADIENT=gradient_diff'
;;;;              printf, lun, 'MAX_LOCAL_SHIFT='+string(maxshift,format='(I0)')
;;;;              printf, lun, 'NEW_CONSTRAINTS'
;;;;              printf, lun, 'FILE_TYPE='+self.filetype
;;;;              if self.filetype eq 'ANA' then begin
;;;;                  printf, lun, 'DATA_TYPE=FLOAT'
;;;;              endif 
;;;;              printf, lun, 'FAST_QR'
;;;;              IF self.filetype EQ 'MOMFBD' THEN BEGIN
;;;;                  printf, lun, 'GET_PSF'
;;;;                  printf, lun, 'GET_PSF_AVG'
;;;;              ENDIF
;;;;              printf, lun, 'FPMETHOD=horint'
;;;;              printf, lun, 'SIM_X='+sim_x_string
;;;;              printf, lun, 'SIM_Y='+sim_y_string
;;;;   
;;;;              ;; External keywords?
;;;;              if(keyword_set(global_keywords)) then begin
;;;;                 nk = n_elements(global_keywords)
;;;;                 for ki = 0L, nk -1 do printf, lun, global_keywords[ki]
;;;;              endif
;;;;   
;;;;              free_lun, lun
;;;;              
;;;;           endfor                  ; ipref
;;;;        endfor                     ; iscan
;;;;     endfor                        ; idir
;;;;   
;;;;     print, inam+' : done!'
  
end
