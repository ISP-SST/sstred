; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CRISP pipeline
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
;    wb_states : 
;   
;   
;   
;    numpoints : in, optional, type=integer, default=88
;   
;      The size of MOMFBD subfields.
;   
;    modes : in, optional, type=string, default= '2-45,50,52-55,65,66'
;   
;      The modes to include in the expansions of the wavefront phases.
;   
;    date_obs : in, optional, type=string
;   
;      The date of observations in ISO (YYYY-MM-DD) format. If not
;      given, taken from class object.
;   
;    state : 
;   
;   
;   
;    no_descatter : in, optional, type=boolean
;   
;       Set this if your data is from a near-IR (777 or 854 nm) line
;       and you do not want to do backscatter corrections.
;   
;    global_keywords : in, optional, type=strarr
;   
;      Any global keywords that you want to add to the momfbd config file.
;   
;    unpol : 
;   
;   
;   
;    skip : 
;   
;   
;   
;    pref : 
;   
;   
;   
;    escan : 
;   
;   
;   
;    div : 
;
;    no_pd : in, optional, type=boolean
;   
;       Set this to exclude phase diversity data and processing. 
;   
;    nremove : 
;   
;   
;   
;    oldgains :
;   
;    momfbddir :  in, optional, type=string, default='momfbd'
;   
;       Top directory of output tree.
; 
;
;
;    margin : in, optional, type=integer, default=5
; 
;      A margin (in pixels) to disregard from the FOV edges when
;      constructing the grid of MOMFBD subfields.
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
;-
pro red::prepmomfbd, wb_states = wb_states $
                     , numpoints = numpoints $
                     , modes = modes $
                     , date_obs = date_obs $
                     , dirs = dirs $
                     , state = state $
                     , no_descatter = no_descatter $
                     , global_keywords = global_keywords $
                     , unpol = unpol $
                     , skip = skip $
                     , pref = pref $
                     , escan = escan $
                     , div = div $
                     , nremove = nremove $
                     , oldgains = oldgains $
                     , nfac = nfac $
                     , weight = weight $
                     , maxshift = maxshift $
                     , momfbddir = momfbddir $
                     , margin = margin $
                     , no_pd = no_pd $
                     , refcam = refcam $
                     , extraclip = extraclip $
                     , redux = redux

  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  ;; Logging
  help, /obj, self, output = selfinfo 
  red_writelog, selfinfo = selfinfo


  if(n_elements(msh) eq 0) then maxshift='30'

  ;; Get keywords
  if n_elements(momfbddir) eq 0 then momfbddir = 'momfbd' 
  if n_elements(date_obs) eq 0 then date_obs = self.isodate
  if n_elements(modes) eq 0 then modes = '2-45,50,52-55,65,66'
  if n_elements(nremove) eq 0 then nremove=0
                                ;if n_elements(nfac) eq 0 then nfac = 1.
  if n_elements(nfac) eq 1 then nfac = replicate(nfac,3)

  if n_elements(margin) eq 0 then margin = 5

  if n_elements(dirs) gt 0 then begin
    dirs = [dirs] 
  endif else begin
    if ~ptr_valid(self.data_dirs) then begin
      print, inam+' : ERROR : undefined data_dir'
      return
    endif
    dirs = *self.data_dirs
  endelse

  Ndirs = n_elements(dirs)
  if( Ndirs eq 0) then begin
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

  ;; Cameras
  cams = *self.cameras
  iswb = strmatch(cams,'*-W') or strmatch(cams,'*-D')
  ispd = strmatch(cams,'*-D')

  if keyword_set(no_pd) then begin
    indx = where(~ispd)
    cams = cams[indx]
    iswb = iswb[indx]
    ispd = ispd[indx]
  endif

  Ncams = n_elements(cams)      ; Number of cameras

  if(n_elements(refcam) eq 0) then begin
    indx = where(iswb)
    if max(indx) ge 0 then refcam = indx[0] $
    else refcam = self.refcam
  endif

  if refcam ge Ncams then begin
    print, inam, ' : index of reference camera out of range: ', refcam, ' >= ', Ncams
    return
  endif

  refcam_name = cams[refcam]
                                ; NB: this will overwrite exising offset files !!
  self->getalignment, align=align, cams=cams, refcam=refcam, prefilters=pref $
                      , extraclip=extraclip, /overwrite
  ref_idx = where( align.state1.camera eq refcam_name )
  if max(ref_idx) lt 0 then begin
    print, inam, ' : Failed to get alignment for refererence camera: ', refcam_name
    return
  endif

  detectors = strarr(Ncams)
  for icam = 0, Ncams-1 do detectors[icam] = self -> getdetector(cams[icam])
  ;;self -> getdetectors, dir = self.data_dir

  ;; Print cams
  print, inam + ' : cameras found:'
  for icam = 0, Ncams-1 do begin
    outstr = cams[icam] + ' ' + detectors[icam] + ' '
    if iswb[icam] then outstr += 'WB 'else outstr += 'NB '
    if ispd[icam] then outstr += 'PD '
    print, outstr
  endfor

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


  ref_clip = align[0].clip
  xsz = abs(ref_clip[0,0]-ref_clip[1,0])+1
  ysz = abs(ref_clip[2,0]-ref_clip[3,0])+1
  this_margin = max([0, min([xsz/3, ysz/3, margin])]) ; prevent silly margin values
                                ; generate patch positions with margin
  sim_x = rdx_segment( this_margin, xsz-this_margin, numpoints, /momfbd )
  sim_y = rdx_segment( this_margin, ysz-this_margin, numpoints, /momfbd )
  sim_x_string = strjoin(strtrim(sim_x,2), ',')
  sim_y_string = strjoin(strtrim(sim_y,2), ',')

  for idir=0L, Ndirs-1 do begin
    
    dir = dirs[idir]
    folder_tag = file_basename(dir)
    
    if file_test(dir + refcam_name + '_nostate/',/directory) then subdir = refcam_name + '_nostate/'

    print, inam + ' : Search for reference files in ' + dir
    self->selectfiles, cam=refcam_name, dirs=dir, prefilter=pref, subdir=subdir, $
                       files=ref_files, states=ref_states, nremove=remove, /force ;, /strip_wb

    if n_elements(ref_states) eq 0 then begin
      print, inam, ' : Failed to find files/states for the reference channel in ', dir
      return
    endif
    
    ref_img_dir = file_dirname(file_expand_path(ref_states[0].filename),/mark)
    ref_caminfo = red_camerainfo(detectors[refcam])

    ;; unique prefilters
    upref = ref_states[uniq(ref_states.prefilter, sort(ref_states.prefilter))].prefilter
    Nprefs = n_elements(upref)

    ;; unique scan numbers
    uscan = ref_states[uniq(ref_states.scannumber, sort(ref_states.scannumber))].scannumber
    Nscans = n_elements(uscan)
    
    ;; base output location
    cfg_base_dir = self.out_dir + PATH_SEP() + momfbddir + PATH_SEP() + folder_tag
    
    for iscan=0L, Nscans-1 do begin
      
      if n_elements(escan) ne 0 then if iscan ne escan then continue 

      scannumber = uscan[iscan]
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
        fn_template = strmid(filename, 0, pos) + '%07d' + strmid(filename, pos+len)
        
        self -> get_calib, ref_states[ref_sel[0]] $
                           , gainname = gainname, darkname = darkname, status = status
        if( status lt 0 ) then continue

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
        cfg.globals += 'DATE_OBS=' + date_obs + string(10b)
        cfg.globals += 'PROG_DATA_DIR=./data/' + string(10b)
        cfg.globals += 'NEW_CONSTRAINTS' + string(10b)
        cfg.globals += 'FAST_QR' + string(10b)
        cfg.globals += 'FPMETHOD=horint' + string(10b)
        cfg.globals += 'BASIS=Karhunen-Loeve' + string(10b)
        cfg.globals += 'GETSTEP=getstep_conjugate_gradient' + string(10b)
        cfg.globals += 'GRADIENT=gradient_diff' + string(10b)
        cfg.globals += 'MODES=' + modes + string(10b)
        cfg.globals += 'TELESCOPE_D=0.97' + string(10b)
        cfg.globals += 'MAX_LOCAL_SHIFT='+string(maxshift,format='(I0)') + string(10b)
        cfg.globals += 'NUM_POINTS=' + strtrim(numpoints,2) + string(10b)
        cfg.globals += 'ARCSECPERPIX=' + self.image_scale + string(10b)
        cfg.globals += 'PIXELSIZE=' + strtrim(ref_caminfo.pixelsize, 2) + string(10b)
        cfg.globals += 'FILE_TYPE=' + self.filetype + string(10b)
        if self.filetype eq 'ANA' then begin
          cfg.globals += 'DATA_TYPE=FLOAT' + string(10b)
        endif 
        if self.filetype eq 'MOMFBD' then begin
          cfg.globals += 'GET_PSF' + string(10b)
          cfg.globals += 'GET_PSF_AVG' + string(10b)
        endif
        cfg.globals += 'SIM_X=' + sim_x_string + string(10b)
        cfg.globals += 'SIM_Y=' + sim_y_string + string(10b)

        ;; External keywords?
        if(keyword_set(global_keywords)) then begin
          nk = n_elements(global_keywords)
          for ki=0L, nk-1 do cfg.globals += global_keywords[ki] + string(10b)
        endif
        
        ;; Reference object
        cfg.objects += 'object{' + string(10b)
        cfg.objects += '    WAVELENGTH=' + strtrim(ref_states[ref_sel[0]].pf_wavelength,2) + string(10b)
        cfg.objects += '    OUTPUT_FILE=results/' + detectors[refcam] $
                       + '_' + date_obs+'T'+folder_tag $
                       + '_' + scanstring + '_' + upref[ipref] + string(10b)
        if(n_elements(weight) gt 0 ) then $
           cfg.objects += '    WEIGHT=' + strtrim(weight[0],2) + string(10b)
        cfg.objects += '    channel{' + string(10b)
        cfg.objects += '        IMAGE_DATA_DIR=' + ref_img_dir + string(10b)
        cfg.objects += '        FILENAME_TEMPLATE=' + fn_template + string(10b)
        cfg.objects += '        GAIN_FILE=' + gainname + string(10b)
        cfg.objects += '        DARK_TEMPLATE=' + darkname + string(10b)
        cfg.objects += '        DARK_NUM=0000001' + string(10b)
        cfg.objects += '        ALIGN_CLIP=' $
                       + strjoin(strtrim(ref_clip,2),',') + string(10b)
        if( align[0].xoffs_file ne '' && file_test(align[0].xoffs_file)) then $
           cfg.objects += '        XOFFSET='+align[0].xoffs_file + string(10b)
        if( align[0].yoffs_file ne '' && file_test(align[0].yoffs_file)) then $
           cfg.objects += '        YOFFSET='+align[0].yoffs_file + string(10b)
        if( upref[ipref] EQ '8542' OR upref[ipref] EQ '7772' ) AND $
           ~keyword_set(no_descatter) then begin
          self -> loadbackscatter, detectors[refcam], upref[ipref] $
                                   , bgfile = bgf, bpfile = psff
          if(file_test(psff) AND file_test(bgf)) then begin
            cfg.objects += '        PSF=' + psff + string(10b)
            cfg.objects += '        BACK_GAIN=' + bgf + string(10b)
          endif else begin
            print, inam, ' : No backscatter files found for prefilter: ', upref[ipref]
          endelse
        endif
        if(n_elements(nfac) gt 0) then $
           cfg.objects += '        NF=' + red_stri(nfac[0]) + string(10b)
        cfg.objects += '        INCOMPLETE' + string(10b)
        cfg.objects += '    }' + string(10b)
        cfg.objects += '}' + string(10b)

        red_append, cfg_list, cfg
        
      endfor                    ; prefilter loop
      
    endfor                      ; scan loop

    for icam=0L, Ncams-1 do begin
      if icam ne refcam then begin

        ;; get a list of all states for this camera
        self->selectfiles, cam=cams[icam], dirs=dir, files=files, $
                           states=states, nremove=remove, /force

        for iscan=0L, Nscans-1 do begin
          
          if n_elements(escan) ne 0 then if iscan ne escan then continue 

          scannumber = uscan[iscan]
          scanstring = string(scannumber,format='(I05)')

          for ipref=0L, Nprefs-1 do begin
            
            ;; select a subset of the states which matches prefilter & scan number.
            self->selectfiles, cam=cams[icam], dirs=dir, prefilter=upref[ipref], scan=scannumber, $
                               files=files, states=states, nremove=remove, selected=sel

            if( max(sel) lt 0 ) then continue
            
            cfg_dir = cfg_base_dir + '/'+upref[ipref]+'/cfg/'
            cfg_file = cfg_dir+'momfbd_reduc_'+upref[ipref]+'_'+scanstring+'.cfg'
            cfg_idx = where( cfg_list.file eq cfg_file )
            
            if( max(cfg_idx) lt 0 ) then continue
            
            state_list = states[sel]
            ustates = state_list[uniq(state_list.fullstate, sort(state_list.fullstate))]
            Nstates = n_elements(ustates)
            
            ;; Loop over states and add object to cfg_list
            for is=0L, Nstates-1 do begin
              
              thisstate = ustates[is].fpi_state
              state_idx = where(state_list.fpi_state eq thisstate)
              if max(state_idx) lt 0 then continue
              
;              self -> get_calib, state_list[state_idx[0]] $
;                                 , gainname = gainname, darkname = darkname, status = status
;              if( status lt 0 ) then stop ;continue
              
              darkname = self -> filenames('dark', state_list[state_idx[0]], /no_fits)
              if(~keyword_set(unpol)) then begin
                if(keyword_set(oldgains)) then begin
;                  search = self.out_dir+'/gaintables/'+self.camttag + '.' + ustat1[ii] + '*.gain'
                  stop
                  ;; Not implemented in chromis::filenames yet.
                  gainname = a -> filenames('oldgain', state_list[state_idx[0]], /no_fits)
                endif else begin
                  ;;search =
                  ;;self.out_dir+'/gaintables/'+folder_tag+'/'+self.camttag
                  ;;+ '.' + istate+'.gain'
                  gainname = self -> filenames('scangain', state_list[state_idx[0]] $
                                               , timestamp = stregex(state_list[state_idx[0]].filename $
                                                                     ,'[0-9][0-9]:[0-9][0-9]:[0-9][0-9]' $
                                                                     ,/extr), /no_fits)
                endelse
              endif else begin
                  gainname = self -> filenames('gain', state_list[state_idx[0]], /no_fits)
;
;                 search = self.out_dir+'/gaintables/'+self.camttag + $
;                          '.' + strmid(ustat1[ii], idx[0], $
;                                       idx[nidx-1])+ '*unpol.gain'
              endelse
 

              if state_list[state_idx[0]].nframes eq 1 then begin
                if nremove ge n_elements(state_idx) then continue
                if nremove ne 0 then begin
                  *(cfg_list[cfg_idx].framenumbers) = red_strip(*(cfg_list[cfg_idx].framenumbers), state_list[state_idx[0:nremove-1]].framenumber)
                  state_idx = state_idx[nremove:*] 
                endif
              endif
              red_append, *(cfg_list[cfg_idx].framenumbers), state_list[state_idx].framenumber
              
              img_dir = file_dirname(file_expand_path(state_list[state_idx[0]].filename),/mark)
              filename = file_basename(state_list[state_idx[0]].filename)
              pos = STREGEX(filename, '[0-9]{7}', length=len)
              fn_template = strmid(filename, 0, pos) + '%07d' + strmid(filename, pos+len)

              align_idx = where( align.state2.camera eq cams[icam] and $
                                 align.state2.fpi_state eq thisstate)
              if max(align_idx) lt 0 then begin ; no match for state, try only prefilter
                align_idx = where( align.state2.camera eq cams[icam] and $
                                   align.state2.prefilter eq ustates[is].prefilter)
                if max(align_idx) lt 0 then begin
                  ;;print, inam, ' : Failed to get ANY alignment for camera/state ', cams[icam] + ':' + thisstate
                  ;;stop
                  continue
                endif
              endif

              if n_elements(align_idx) gt 1 then align_idx = align_idx[0] ; just pick the first one for now
              state_align = align[align_idx]
              
              ;; create cfg object
              cfg_list[cfg_idx].objects += 'object{' + string(10b)
              cfg_list[cfg_idx].objects += '    WAVELENGTH=' + strtrim(ustates[is].pf_wavelength,2) + string(10b)
              cfg_list[cfg_idx].objects += '    OUTPUT_FILE=results/' + detectors[icam] $
                                           + '_' + date_obs+'T'+folder_tag $
                                           + '_' + scanstring + '_'+ustates[is].fullstate + string(10b)
              if(n_elements(weight) gt 1) then $
                 cfg_list[cfg_idx].objects += '    WEIGHT=' + strtrim(weight[1],2) + string(10b)
              cfg_list[cfg_idx].objects += '    channel{' + string(10b)
              cfg_list[cfg_idx].objects += '        IMAGE_DATA_DIR=' + img_dir + string(10b)
              cfg_list[cfg_idx].objects += '        FILENAME_TEMPLATE=' + fn_template + string(10b)
              cfg_list[cfg_idx].objects += '        GAIN_FILE=' + gainname + string(10b)
              cfg_list[cfg_idx].objects += '        DARK_TEMPLATE=' + darkname + string(10b)
              cfg_list[cfg_idx].objects += '        DARK_NUM=0000001' + string(10b)
              cfg_list[cfg_idx].objects += '        ALIGN_CLIP=' $
                                           + strjoin(strtrim(state_align.clip,2),',') + string(10b)
              if( state_align.xoffs_file ne '' && file_test(state_align.xoffs_file)) then $
                 cfg_list[cfg_idx].objects += '        XOFFSET='+state_align.xoffs_file + string(10b)
              if( state_align.yoffs_file ne '' && file_test(state_align.yoffs_file)) then $
                 cfg_list[cfg_idx].objects += '        YOFFSET='+state_align.yoffs_file + string(10b)
              if( ustates[is].prefilter EQ '8542' OR ustates[is].prefilter EQ '7772' ) AND $
                 ~keyword_set(no_descatter) then begin
                self -> loadbackscatter, detectors[icam], ustates[is].prefilter, bgfile = bgf, bpfile = psff
                if(file_test(psff) AND file_test(bgf)) then begin
                  cfg_list[cfg_idx].objects += '        PSF=' + psff + string(10b)
                  cfg_list[cfg_idx].objects += '        BACK_GAIN=' + bgf + string(10b)
                endif else begin
                  print, inam, ' : No backscatter files found for prefilter: ' $
                         , ustates[is].prefilter
                endelse
              endif
              if(n_elements(nfac) gt 1) then $
                 cfg_list[cfg_idx].objects += '        NF=' + red_stri(nfac[1]) + string(10b)
              if keyword_set(redux) && max(nremove) gt 0 then $
                 cfg_list[cfg_idx].objects += '        DISCARD=' $
                                              + strjoin(strtrim(nremove,2),',') + string(10b)
              cfg_list[cfg_idx].objects += '        INCOMPLETE' + string(10b)
              cfg_list[cfg_idx].objects += '    }' + string(10b)
              cfg_list[cfg_idx].objects += '}' + string(10b)
              
              if(keyword_set(wb_states)) then begin
                
                ;; select WB files with the same framenumbers
                self->selectfiles, prefilter=upref[ipref], scan=scannumber, $
                                   files=ref_files, states=ref_states, selected=ref_sel $
                                   , framenumbers = state_list[state_idx].framenumber
                
                if( max(ref_sel) lt 0 ) then continue
                
                self -> get_calib, ref_states[ref_sel[0]] $
                                   , gainname = gainname, darkname = darkname, status = status
                if( status lt 0 ) then continue
                
                fullname = file_readlink(ref_states[ref_sel[0]].filename)
                img_dir = file_dirname(fullname,/mark)
                filename = file_basename(fullname)
                pos = STREGEX(filename, '[0-9]{7}', length=len)
                fn_template = strmid(filename, 0, pos) + '%07d' + strmid(filename, pos+len)
                
                cfg_list[cfg_idx].objects += 'object{' + string(10b)
                cfg_list[cfg_idx].objects += '    WAVELENGTH=' $
                                             + strtrim(ustates[is].pf_wavelength,2) + string(10b)
                if(n_elements(weight) gt 2) then $
                   cfg.objects += '    WEIGHT=' + strtrim(weight[3],2) + string(10b) $
                else cfg_list[cfg_idx].objects += '    WEIGHT=0.00' + string(10b)
                cfg_list[cfg_idx].objects += '    OUTPUT_FILE=results/'+detectors[refcam] + '_' $
                                             + date_obs+'T'+folder_tag $
                                             + '_' + scanstring + '_'+ustates[is].fullstate + string(10b)
                cfg_list[cfg_idx].objects += '    channel{' + string(10b)
                cfg_list[cfg_idx].objects += '        IMAGE_DATA_DIR=' + img_dir + string(10b)
                cfg_list[cfg_idx].objects += '        FILENAME_TEMPLATE=' + fn_template + string(10b)
                cfg_list[cfg_idx].objects += '        GAIN_FILE=' + gainname + string(10b)
                cfg_list[cfg_idx].objects += '        DARK_TEMPLATE=' + darkname + string(10b)
                cfg_list[cfg_idx].objects += '        DARK_NUM=0000001' + string(10b)
                cfg_list[cfg_idx].objects += '        ALIGN_CLIP=' + strjoin(strtrim(ref_clip,2),',') + string(10b)
                if( ustates[is].prefilter EQ '8542' OR $
                    ustates[is].prefilter EQ '7772' ) AND ~keyword_set(no_descatter) then begin
                  self -> loadbackscatter, detectors[refcam], ustates[is].prefilter $
                                           , bgfile = bgf, bpfile = psff
                  if(file_test(psff) AND file_test(bgf)) then begin
                    cfg_list[cfg_idx].objects += '        PSF=' + psff + string(10b)
                    cfg_list[cfg_idx].objects += '        BACK_GAIN=' + bgf + string(10b)
                  endif else begin
                    print, inam, ' : No backscatter files found for prefilter: ' $
                           , ustates[is].prefilter
                  endelse
                endif
                if(n_elements(nfac) gt 2) then cfg_list[cfg_idx].objects $
                   += '        NF=' + red_stri(nfac[2]) + string(10b)
                if keyword_set(redux) && max(nremove) gt 0 then $
                   cfg_list[cfg_idx].objects += '        DISCARD=' $
                                                + strjoin(strtrim(nremove,2),',') + string(10b)
                cfg_list[cfg_idx].objects += '        INCOMPLETE' + string(10b)
                cfg_list[cfg_idx].objects += '    }' + string(10b)
                cfg_list[cfg_idx].objects += '}' + string(10b)
              endif
              
            endfor
            
          endfor                ; prefilter loop
          
        endfor                  ; scan loop
        
      endif                     ; icam ne refcam
      
    endfor                      ; cam loop
    
  endfor                        ; dir loop
  
  for ic=0, n_elements(cfg_list)-1 do begin
    
    if( ~file_test(cfg_list[ic].dir, /directory) ) then begin
      file_mkdir, cfg_list[ic].dir+'/data/'
      file_mkdir, cfg_list[ic].dir+'/results/'
    endif
    
;    number_str = string(cfg_list[ic].first_file, format='(I07)') $
;                 + '-' + string(cfg_list[ic].last_file,format='(I07)')
;        cfg_list[ic].globals += 'IMAGE_NUMS=' + number_str +
;        string(10b) + string(10b)
    frmnums = *(cfg_list[ic].framenumbers)
    frmnums = frmnums[uniq(frmnums, sort(frmnums))]
    cfg_list[ic].globals += 'IMAGE_NUMS=' + red_collapserange(frmnums, ld='', rd='')
    
    print,'Writing: ', cfg_list[ic].file
    openw, lun, cfg_list[ic].file, /get_lun, width=2500
    printf, lun, cfg_list[ic].objects + cfg_list[ic].globals
    free_lun, lun
    
    ptr_free, cfg_list[ic].framenumbers

  endfor                        ; cfg loop
  

return
  
  
  
  
  for idir = 0L, Ndirs - 1 do begin

     data_dir = dirs[idir]
     folder_tag = file_basename(data_dir)
     searchcamdir = data_dir + '/' + searchcam + '/'

;     search = self.out_dir+'/data/'+folder_tag+'/'+self.camt
     search = searchcamdir + '*' + searchdet + '*'
     files = file_search(search, count = Nfiles) 
     
     IF Nfiles EQ 0 THEN BEGIN
         print, inam + ' : ERROR -> no frames found : '+search
         print, inam + '   Did you run link_data?'
         return
     ENDIF 

     files = red_sortfiles(temporary(files))
     
     ;; Get image unique states
     self -> extractstates, files, states
;     stat = red_getstates(files, /LINKS)
     
     ;; skip leading frames?
     IF nremove GT 0 THEN red_flagtuning, stat, nremove

     ;; Get unique prefilters
     upref = states[uniq(states.prefilter, sort(states.prefilter))].prefilter
     Nprefs = n_elements(upref)

     ;; Get scan numbers
     uscan = states[uniq(states.scannumber, sort(states.scannumber))].scannumber
     Nscans = n_elements(uscan)

     ;;states = stat.hscan+'.'+stat.state
     ;;pos = uniq(states, sort(states))
     ;;ustat = stat.state[pos]
     ;;ustatp = stat.pref[pos]
     ;;                           ;ustats = stat.scan[pos]

     ;;ntt = n_elements(ustat)
     ;;hscans = stat.hscan[pos]

     ;; Create a reduc file per prefilter and scan number?
     outdir0 = self.out_dir + '/' + momfbddir + '/' + folder_tag

     ;; Choose offset state
     for iscan = 0L, Nscans-1 do begin

        IF n_elements(escan) NE 0 THEN IF iscan NE escan THEN CONTINUE 

        scannumber = uscan[iscan]
        scanstring = string(scannumber,format='(I05)')

        for ipref = 0L, Nprefs-1 do begin

           if(keyword_set(pref)) then begin
              if(upref[ipref] NE pref) then begin
                 print, inam + ' : Skipping prefilter -> ' + upref[ipref]
                 continue
              endif
           endif

           ;; Load align clips
;            clipfile = self.out_dir + '/calib/align_clips.'+upref[ipref]+'.sav'
;            IF(~file_test(clipfile)) THEN BEGIN
;               print, inam + ' : ERROR -> align_clip file not found'
;               print, inam + ' : -> you must run red::getalignclips first!'
;               continue
;            endif
;            restore, clipfile
;            wclip = acl[0]
;            tclip = acl[1]
;            rclip = acl[2]

           wclip = align[0].clip
           xsz = abs(wclip[0,0]-wclip[1,0])+1
           ysz = abs(wclip[2,0]-wclip[3,0])+1
           this_margin = max([0, min([xsz/3, ysz/3, margin])])  ; prevent silly margin values
           ; generate patch positions with margin
           sim_x = rdx_segment( this_margin, xsz-this_margin, numpoints, /momfbd )
           sim_y = rdx_segment( this_margin, ysz-this_margin, numpoints, /momfbd )
           sim_x_string = strjoin(strtrim(sim_x,2), ',')
           sim_y_string = strjoin(strtrim(sim_y,2), ',')

           lam = strmid(string(float(upref[ipref]) * 1.e-10), 2)

           cfg_file = 'momfbd.reduc.'+upref[ipref]+'.'+scanstring+'.cfg'
           outdir = outdir0 + '/'+upref[ipref]+'/cfg/'
           rdir = outdir + 'results/'
           ddir = outdir + 'data/'
           if( ~file_test(rdir, /directory) ) then begin
               file_mkdir, rdir
               file_mkdir, ddir
           endif
           if(n_elements(lun) gt 0) then free_lun, lun
           openw, lun, outdir + cfg_file, /get_lun, width=2500

           ;; Image numbers
           numpos = where((states.scannumber eq uscan[iscan]) AND (states.skip eq 0B) AND (states.prefilter eq upref[ipref]), ncount)
           if( ncount eq 0 ) then continue
           n0 = min(states[numpos].framenumber)
           n1 = max(states[numpos].framenumber)
           numbers = states[numpos].framenumber
           number_range = [ min(states[numpos].framenumber), max(states[numpos].framenumber)]
           number_str = string(number_range[0],format='(I07)') + '-' + string(number_range[1],format='(I07)')
           ;nall = strjoin(strtrim(states[numpos].framenumber,2),',')
           print, inam+' : Prefilter = '+upref[ipref]+' -> scan = '+scanstring+' -> image range = ['+number_str+']'

           self -> get_calib, states[numpos[0]], gainname = gainname, darkname = darkname, status = status
           if status ne 0 then begin
               print, inam+' : no dark/gain found for camera: ', refcam_name
               continue
           endif
           
           ;; WB anchor channel
           printf, lun, 'object{'
           printf, lun, '  WAVELENGTH=' + lam
           printf, lun, '  OUTPUT_FILE=results/'+align[0].state1.detector+'.'+scanstring+'.'+upref[ipref]
           if(n_elements(weight) eq 3) then printf, lun, '  WEIGHT='+string(weight[0])
           printf, lun, '  channel{'
           printf, lun, '    IMAGE_DATA_DIR='+self.out_dir+'/data/'+folder_tag+ '/' +align[0].state1.camera+'_nostate/'
           printf, lun, '    FILENAME_TEMPLATE='+align[0].state1.detector+'.'+scanstring+'.'+upref[ipref]+'.%07d'
           printf, lun, '    GAIN_FILE=' + gainname
           printf, lun, '    DARK_TEMPLATE=' + darkname
           printf, lun, '    DARK_NUM=0000001'
           printf, lun, '    ALIGN_CLIP=' + strjoin(strtrim(wclip,2),',')
           if (upref[ipref] EQ '8542' OR upref[ipref] EQ '7772' ) AND ~keyword_set(no_descatter) then begin
              self -> loadbackscatter, align[0].state1.detector, upref[ipref], bgfile = bgf, bpfile = psff
;              psff = self.descatter_dir+'/'+align[0].state1.detector+'.psf.f0'
;              bgf = self.descatter_dir+'/'+align[0].state1.detector+'.backgain.f0'
;              if(file_test(psff) AND file_test(bgf)) then begin
              printf, lun, '    PSF='+psff
              printf, lun, '    BACK_GAIN='+bgf
;              endif
           endif 

           if(keyword_set(div)) then begin
              printf, lun, '    DIVERSITY='+string(div[0])+' mm'
           endif
           
           if( align[0].xoffs_file ne '' && file_test(align[0].xoffs_file)) then printf, lun, '    XOFFSET='+align[0].xoffs_file
           if( align[0].yoffs_file ne '' && file_test(align[0].yoffs_file)) then printf, lun, '    YOFFSET='+align[0].yoffs_file

           
           if(n_elements(nfac) gt 0) then printf,lun,'    NF=',red_stri(nfac[0])
           printf, lun, '  }'
           printf, lun, '}'

stop
    ;; Loop all wavelengths
           pos1 = where((ustatp eq upref[ipref]), count)
           if(count eq 0) then continue
           ustat1 = ustat[pos1]

           for ii = 0L, count - 1 do BEGIN
              
              ;; External states?
              if(keyword_set(state)) then begin
                 dum = where(state eq ustat1[ii], cstate)
                 if(cstate eq 0) then continue
                 print, inam+' : found '+state+' -> scan = '+scanstring
              endif

              self -> whichoffset, ustat1[ii], xoff = xoff, yoff = yoff

              ;; Trans. camera
              istate = red_encode_scan(hscans[pos1[ii]], scan)+'.'+ustat1[ii]

              ;; lc4?
              tmp = strsplit(istate,'.', /extract)
              ntmp = n_elements(tmp)

              idx = strsplit(ustat1[ii],'.')
              nidx = n_elements(idx)
              iwavt = strmid(ustat1[ii], idx[0], idx[nidx-1]-1)

              if(keyword_set(skip)) then begin
                 dum = where(iwavt eq skip, ccout)
                 if ccout ne 0 then begin
                    print, inam+' : skipping state -> '+ustat1[ii]
                    continue
                 endif
              endif

              printf, lun, 'object{'
              printf, lun, '  WAVELENGTH=' + lam
              printf, lun, '  OUTPUT_FILE=results/'+self.camttag+'.'+istate 
              if(n_elements(weight) eq 3) then printf, lun, '  WEIGHT='+string(weight[1])
              printf, lun, '  channel{'
              printf, lun, '    IMAGE_DATA_DIR='+self.out_dir+'/data/'+folder_tag+ '/' +self.camt+'/'
              printf, lun, '    FILENAME_TEMPLATE='+self.camttag+'.'+istate+'.%07d'

              if(~keyword_set(unpol)) then begin
                 if(keyword_set(oldgains)) then begin
                    search = self.out_dir+'/gaintables/'+self.camttag + '.' + ustat1[ii] + '*.gain'
                 endif else begin
                    search = self.out_dir+'/gaintables/'+folder_tag+'/'+self.camttag + '.' + istate+'.gain'
                 endelse
              endif Else begin

                 search = self.out_dir+'/gaintables/'+self.camttag + $
                          '.' + strmid(ustat1[ii], idx[0], $
                                       idx[nidx-1])+ '*unpol.gain'
              endelse
              printf, lun, '    GAIN_FILE=' + file_search(search)
              printf, lun, '    DARK_TEMPLATE='+self.out_dir+'/darks/'+self.camttag+'.summed.0000001'
              printf, lun, '    DARK_NUM=0000001'
              printf, lun, '    ' + tclip

              xofile = self.out_dir+'/calib/'+self.camttag+'.'+xoff
              yofile = self.out_dir+'/calib/'+self.camttag+'.'+yoff
              if(file_test(xofile)) then printf, lun, '    XOFFSET='+xofile
              if(file_test(yofile)) then printf, lun, '    YOFFSET='+yofile

              if (upref[ipref] EQ '8542' OR upref[ipref] EQ '7772' ) AND ~keyword_set(no_descatter) then begin
                 self -> loadbackscatter, self.camttag, upref[ipref], bgfile = bgf, bpfile = psff
;                 psff = self.descatter_dir+'/'+self.camttag+'.psf.f0'
;                 bgf = self.descatter_dir+'/'+self.camttag+'.backgain.f0'
;                 if(file_test(psff) AND file_test(bgf)) then begin
                 printf, lun, '    PSF='+psff
                 printf, lun, '    BACK_GAIN='+bgf
;              endif
              endif 

              if(keyword_set(div)) then begin
                 printf, lun, '    DIVERSITY='+string(div[1])+' mm'
              endif
              if(n_elements(nfac) gt 0) then printf,lun,'    NF=',red_stri(nfac[1])
           
              printf, lun, '    INCOMPLETE'
              printf, lun, '  }'
              printf, lun, '}'  

              ;; Reflected camera
              printf, lun, 'object{'
              printf, lun, '  WAVELENGTH=' + lam
              printf, lun, '  OUTPUT_FILE=results/'+self.camrtag+'.'+istate 
              if(n_elements(weight) eq 3) then printf, lun, '  WEIGHT='+string(weight[2])
              printf, lun, '  channel{'
              printf, lun, '    IMAGE_DATA_DIR='+self.out_dir+'/data/'+folder_tag+ '/' +self.camr+'/'
              printf, lun, '    FILENAME_TEMPLATE='+self.camrtag+'.'+istate+'.%07d'
                                ;   printf, lun, '    DIVERSITY=0.0 mm' 
              if(~keyword_set(unpol)) then begin
                 if(keyword_set(oldgains)) then begin
                    search = self.out_dir+'/gaintables/'+self.camrtag + '.' + ustat1[ii] + '*.gain'
                 endif else begin
                    search = self.out_dir+'/gaintables/'+folder_tag+'/'+self.camrtag + '.' + istate+'.gain'
                 endelse
              endif Else begin
                 idx = strsplit(ustat1[ii],'.')
                 nidx = n_elements(idx)
                 search = file_search(self.out_dir+'/gaintables/'+self.camrtag + $
                                      '.' + strmid(ustat1[ii], idx[0], $
                                                   idx[nidx-1])+ '*unpol.gain')
                                ;if tmp[ntmp-1] eq 'lc4' then search = self.out_dir+'/gaintables/'+$
                                ;                                      self.camrtag + '.' + ustat[pos[ii]] + $
                                ;                                      '*.gain'
              endelse
              printf, lun, '    GAIN_FILE=' + file_search(search)
              printf, lun, '    DARK_TEMPLATE='+self.out_dir+'/darks/'+self.camrtag+'.summed.0000001'
              printf, lun, '    DARK_NUM=0000001'
              printf, lun, '    ' + rclip
              xofile = self.out_dir+'/calib/'+self.camrtag+'.'+xoff
              yofile = self.out_dir+'/calib/'+self.camrtag+'.'+yoff
              if(file_test(xofile)) then printf, lun, '    XOFFSET='+xofile
              if(file_test(yofile)) then printf, lun, '    YOFFSET='+yofile
                                ;
              if (upref[ipref] EQ '8542' OR upref[ipref] EQ '7772' ) AND ~keyword_set(no_descatter) then begin
                 self -> loadbackscatter, self.camrtag, upref[ipref], bgfile = bgf, bpfile = psff
;                 psff = self.descatter_dir+'/'+self.camrtag+'.psf.f0'
;                 bgf = self.descatter_dir+'/'+self.camrtag+'.backgain.f0'
;                 if(file_test(psff) AND file_test(bgf)) then begin
                 printf, lun, '    PSF='+psff
                 printf, lun, '    BACK_GAIN='+bgf
;                 endif
              endif 

              if(keyword_set(div)) then begin
                 printf, lun, '    DIVERSITY='+string(div[2])+' mm'
              endif
              if(n_elements(nfac) gt 0) then printf,lun,'    NF=',red_stri(nfac[2])

              printf, lun, '    INCOMPLETE'
              printf, lun, '  }'
              printf, lun, '}'

              ;; WB with states (for de-warping to the anchor, only to
              ;; remove rubbersheet when differential seeing is
              ;; strong)
              if(keyword_set(wb_states)) then begin
                 printf, lun, 'object{'
                 printf, lun, '  WAVELENGTH=' + lam
                 printf, lun, '  WEIGHT=0.00'
                 printf, lun, '  OUTPUT_FILE=results/'+align[0].state1.detector+'.'+istate 
                 printf, lun, '  channel{'
                 printf, lun, '    IMAGE_DATA_DIR='+self.out_dir+'/data/'+folder_tag+ '/' +align[0].state1.camera+'/'
                 printf, lun, '    FILENAME_TEMPLATE='+align[0].state1.detector+'.'+istate+'.%07d'
                                ; printf, lun, '    DIVERSITY=0.0 mm'
                 printf, lun, '    GAIN_FILE=' + file_search(self.out_dir+'/gaintables/'+align[0].state1.detector + $
                                                             '.' + upref[ipref] + '*.gain')
                 printf, lun, '    DARK_TEMPLATE='+self.out_dir+'/darks/'+align[0].state1.detector+'.summed.0000001'
                 printf, lun, '    DARK_NUM=0000001'
                 printf, lun, '    ' + wclip
                 
                 if (upref[ipref] EQ '8542' OR upref[ipref] EQ '7772' ) AND ~keyword_set(no_descatter) then begin
                    self -> loadbackscatter, align[0].state1.detector, upref[ipref], bgfile = bgf, bpfile = psff
;                    psff = self.descatter_dir+'/'+align[0].state1.detector+'.psf.f0'
;                    bgf = self.descatter_dir+'/'+align[0].state1.detector+'.backgain.f0'
;                    if(file_test(psff) AND file_test(bgf)) then begin
                    printf, lun, '    PSF='+psff
                    printf, lun, '    BACK_GAIN='+bgf
;                    endif
                 endif 

                 if(keyword_set(div)) then begin
                    printf, lun, '    DIVERSITY='+string(div[0])+' mm'
                 endif
                 xofile = self.out_dir+'/calib/'+align[0].state1.detector+'.'+xoff
                 yofile = self.out_dir+'/calib/'+align[0].state1.detector+'.'+yoff
                 if(file_test(xofile)) then printf, lun, '    XOFFSET='+xofile
                 if(file_test(yofile)) then printf, lun, '    YOFFSET='+yofile
                 if(n_elements(nfac) gt 0) then printf,lun,'    NF=',red_stri(nfac[0])

                 printf, lun, '    INCOMPLETE'
                 printf, lun, '  }'
                 printf, lun, '}'
              endif          
           endfor               ; ii

           ;; Global keywords
           printf, lun, 'PROG_DATA_DIR=./data/'
           printf, lun, 'DATE_OBS='+date_obs
           printf, lun, 'IMAGE_NUMS='+nall       ;;  n0+'-'+n1
           printf, lun, 'BASIS=Karhunen-Loeve'
           printf, lun, 'MODES='+modes
           printf, lun, 'NUM_POINTS='+strtrim(numpoints,2)
           printf, lun, 'TELESCOPE_D=0.97'
           printf, lun, 'ARCSECPERPIX='+self.image_scale
           printf, lun, 'PIXELSIZE=16.0E-6'
           printf, lun, 'GETSTEP=getstep_conjugate_gradient'
           printf, lun, 'GRADIENT=gradient_diff'
           printf, lun, 'MAX_LOCAL_SHIFT='+string(maxshift,format='(I0)')
           printf, lun, 'NEW_CONSTRAINTS'
           printf, lun, 'FILE_TYPE='+self.filetype
           if self.filetype eq 'ANA' then begin
               printf, lun, 'DATA_TYPE=FLOAT'
           endif 
           printf, lun, 'FAST_QR'
           IF self.filetype EQ 'MOMFBD' THEN BEGIN
               printf, lun, 'GET_PSF'
               printf, lun, 'GET_PSF_AVG'
           ENDIF
           printf, lun, 'FPMETHOD=horint'
           printf, lun, 'SIM_X='+sim_x_string
           printf, lun, 'SIM_Y='+sim_y_string

           ;; External keywords?
           if(keyword_set(global_keywords)) then begin
              nk = n_elements(global_keywords)
              for ki = 0L, nk -1 do printf, lun, global_keywords[ki]
           endif

           free_lun, lun
           
        endfor                  ; ipref
     endfor                     ; iscan
  endfor                        ; idir

  print, inam+' : done!'
  
end
