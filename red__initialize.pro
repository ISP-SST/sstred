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
; :Returns:
; 
; 
; :Params:
; 
;    filename : 
;   
;   
;   
; 
; :Keywords:
; 
;    develop : in, optional, type=boolean
; 
;       Run in developer mode.
; 
; :History:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2013-12-10 : PS  adapt for multiple flat_dir
;
;   2014-01-08 : MGL. Initialize new fields: isodate, log_dir, telog, pinhole_spacing.
; 
;   2014-01-10 : PS  New config variable filtype
;
;   2014-01-22 : MGL. Adapt to string functions moved to the str_
;                namespace.
;
;   2014-03-21 : MGL. Allow for multiple dark_dir.
;
;   2014-10-16 : MGL. Set up for Coyote graphics.
;
;   2016-02-15 : MGL. Make descatter_dir local.
;
;   2016-08-23 : THI. Rename camtag to detector and channel to camera,
;                so the names match those of the corresponding SolarNet
;                keywords.
;
;   2017-01-25 : MGL. Added (nominal) diversity.
;
;   2017-03-09 : MGL. Version info. New keyword "develop".
;
;   2017-03-10 : MGL. New member "developer_mode". 
;
;   2017-04-26 : MGL. Use red_mpfit_version. 
;
;   2017-06-19 : MGL. If multiple coyote directories, give more info
;                about them. 
;
;   2017-06-20 : MGL. Clean up key-value handling.
;
;   2017-11-08 : THI: Bump the rdx-version dependence. 
;
;   2018-05-31 : MGL. Set self.refcam to camera ending in '-W'. 
;
;   2020-03-11 : MGL. Default for self.direction.
;
;-
pro red::initialize, filename, develop = develop

  ;; Are we running in developer mode?
  self.developer_mode = keyword_set(develop)

  ;; Test file
  if(~file_test(filename)) then begin
    print, 'red::initialize : ERROR : file not found: ', filename
    return
  endif
  self.filename = filename
  self.filetype = 'MOMFBD'
  
  ;; Init vars
  self.polcal_dir = ''
  
  self.dodata = 1B
  self.doflat = 1B
  self.dopinh = 1B
  self.dodark = 1B
  self.dopolcal = 1B
  self.dodescatter = 1B

  ;; open file and get fields
  openr, lun, filename, /get_lun
  nl = 0L
  while(~eof(lun)) do begin

    ;; Read a line
    line = ''
    readf, lun, line                           ; read line

    ;; Ignore empty or commented out lines
    if strmid(line, 0, 1) eq '#' or strlen(strtrim(line, 2)) eq 0 then continue
    keyval = red_parseparameter(line)
    key = keyval[0]
    value = red_cleanvalue(keyval[1])

    ;; Get fields
    case key of
      'root_dir': begin
        self.root_dir = value
        IF strmid(self.root_dir, strlen(self.root_dir)-1) NE '/' THEN $
           self.root_dir += '/'
      end
      'camsz'           : self.camsz = value
      'diversity'       : self.diversity = value
      'direction'       : self.direction = value
      'filetype'        : self.filetype = value
      'image_scale'     : self.image_scale = value
      'isodate'         : self.isodate = value
      'out_dir'         : self.out_dir = value
      'pinhole_spacing' : self.pinhole_spacing = float(value)
      'pixel_size'      : self.pixel_size = value
      'polcal_dir'      : self.polcal_dir = value
      'prefilter_dir'   : self.prefilter_dir = value
      'telescope_d'     : self.telescope_d = value
;        'descatter_dir' : self.descatter_dir = value
      'dark_dir': begin
        if ptr_valid(self.dark_dir) then red_append, *self.dark_dir, value $
        else self.dark_dir = ptr_new(value, /no_copy)
      end
      'flat_dir': begin
        ;;red_append, self.flat_dir, value, /point
        ;;help, self.flat_dir
        ;;stop
        if ptr_valid(self.flat_dir) then red_append, *self.flat_dir, value $
        else self.flat_dir = ptr_new(value, /no_copy)
      end
      'pinh_dir': begin
        if ptr_valid(self.pinh_dirs) then red_append, *self.pinh_dirs, value $
        else self.pinh_dirs = ptr_new(value, /no_copy)
      end
      'data_dir': begin
        ;;red_append, self.data_dirs, value, /point
        if ptr_valid(self.data_dirs) then red_append, *self.data_dirs, value $
        else self.data_dirs = ptr_new(value, /no_copy)
      end
      'camera': begin
        if ptr_valid(self.cameras) then red_append, *self.cameras, value $
        else self.cameras = ptr_new(value, /no_copy)
      end
      else: begin
        print, 'red::initialize : Skipping line ',$
               strcompress(string(nl),/remove_all), ': ',line
      end
    endcase
    nl+= 1
  endwhile
  free_lun, lun

  if(strlen(self.root_dir) gt 0) then begin
    print, 'red::initialize : root_dir = '+ self.root_dir
    IF ptr_valid(self.flat_dir) THEN FOR i = 0, n_elements(*self.flat_dir)-1 DO $
       IF strmid((*self.flat_dir)[i], 0, 1) NE '/' AND strlen((*self.flat_dir)[i]) GT 0 $
       THEN (*self.flat_dir)[i] = self.root_dir + (*self.flat_dir)[i]
    IF ptr_valid(self.dark_dir) THEN FOR i = 0, n_elements(*self.dark_dir)-1 DO $
       IF strmid((*self.dark_dir)[i], 0, 1) NE '/' AND strlen((*self.dark_dir)[i]) GT 0 $
       THEN (*self.dark_dir)[i] = self.root_dir + (*self.dark_dir)[i]
    IF ptr_valid(self.data_dirs) THEN FOR i = 0, n_elements(*self.data_dirs)-1 DO $
       IF strmid((*self.data_dirs)[i], 0, 1) NE '/' AND strlen((*self.data_dirs)[i]) GT 0 $
       THEN (*self.data_dirs)[i] = self.root_dir + (*self.data_dirs)[i]
    IF ptr_valid(self.pinh_dirs) THEN FOR i = 0, n_elements(*self.pinh_dirs)-1 DO $
       IF strmid((*self.pinh_dirs)[i], 0, 1) NE '/' AND strlen((*self.pinh_dirs)[i]) GT 0 $
       THEN (*self.pinh_dirs)[i] = self.root_dir + (*self.pinh_dirs)[i]
    ;; if(strmid(self.dark_dir, 0, 1) NE '/' AND strlen(self.dark_dir) gt 0) then $
    ;;   self.dark_dir = self.root_dir + self.dark_dir
    if(strmid(self.prefilter_dir, 0, 1) NE '/' AND strlen(self.prefilter_dir) gt 0) then $
       self.prefilter_dir = self.root_dir + self.prefilter_dir
    ;; if(strmid(self.data_list[self.ndir-1], 0, 1) NE '/' AND strlen(self.data_list[self.ndir-1]) gt 0) then $
    ;;    self.data_list[0:self.ndir-1] = self.root_dir + self.data_list[0:self.ndir-1]
    if(strmid(self.polcal_dir, 0, 1) NE '/' AND strlen(self.polcal_dir) gt 0) then $
       self.polcal_dir = self.root_dir + self.polcal_dir
    ;; if(strmid(self.descatter_dir, 0, 1) NE '/' AND  strlen(self.descatter_dir) gt 0) then $
    ;;     self.descatter_dir = self.root_dir + self.descatter_dir
  endif
  if(strlen(self.telescope_d) eq 0) then self.telescope_d = '0.970'
  print, 'red::initialize : telescope_d = '+self.telescope_d 

  if(strlen(self.image_scale) eq 0) then self.image_scale = '0.0592'
  print, 'red::initialize : image_scale = '+self.image_scale 

  if(strlen(self.pixel_size) eq 0) then self.pixel_size = '16.0E-6'
  print, 'red::initialize : pixel_size = '+self.pixel_size 

  if(strlen(self.camsz) eq 0) then self.camsz = '1024'
  print, 'red::initialize : camsz = '+self.camsz

  if self.pinhole_spacing eq 0.0 then begin
    ;; Default value (5.12 arcseconds between pinholes) was
    ;; measured in 2013 by comparison with SDO/HMI images.
    self.pinhole_spacing = 5.12
  endif 
  print, 'red::initialize : pinhole_spacing = '+strtrim(self.pinhole_spacing, 2)

  ;; check available fields
;  if(self.descatter_dir eq '') then begin
;     print, 'red::initialize : WARNING : descatter_dir is undefined!'
;     self.dodescatter = 0B
;  endif
;  if(self.dark_dir eq '') then begin
  IF ~ptr_valid(self.dark_dir) then begin
    print, 'red::initialize : WARNING : dark_dir is undefined!'
    self.dodark = 0B
  endif
  IF ~ptr_valid(self.flat_dir) then begin
    print, 'red::initialize : WARNING : flat_dir is undefined!'
    self.doflat = 0B
  endif
  IF ~ptr_valid(self.data_dirs) then begin
    print, 'red::initialize : WARNING : data_dir is undefined!'
  endif
  if(self.out_dir eq '') then begin
    print, 'red::initialize : ERROR : out_dir is undefined!'
    return
  endif
  ;;if(self.data_dir eq '') then begin
  ;;   print, 'red::initialize : WARNING : data_dir is undefined!'
  ;;   self.dodata = 0B
  ;;endif
  if( ~ptr_valid(self.pinh_dirs) ) then begin
    print, 'red::initialize : WARNING : pinh_dir is undefined!'
    self.dopinh = 0B
  endif
  if(self.polcal_dir eq '') then begin
    print, 'red::initialize : WARNING : polcal_dir is undefined!'
    self.dopolcal = 0B
  endif
  if(self.isodate eq '') then begin
    print, 'red::initialize : WARNING : isodate is undefined!'
    print, '                  Try to get it from PWD!'
    date = stregex(getenv('PWD'),'[12][0-9][0-9][0-9][-.][0-1][0-9][-.][0-3][0-9]',/extr)
    if date eq '' then begin
      print, 'red::initialize : WARNING : No recognizable date in PWD. Giving up.'
    endif else begin
      ;; Do a red_strreplace in case the found date uses dots rather
      ;; than dashes.
      self.isodate = red_strreplace(date, '.', '-', n = 2)
    endelse
  endif
  
  ;; Fields that depend on fields defined above:
  self.log_dir = self.out_dir+'/downloads/sstlogs/'
  self.telog = self.log_dir+'positionLog_'+red_strreplace(self.isodate, '-', '.', n = 2)+'_final'
  self.descatter_dir = self.out_dir+'/downloads/backscatter/'  
  if ptr_valid(self.cameras) then begin
    self.refcam = where(strmatch(*self.cameras,'*-W'), Nmatch)
    if Nmatch eq 0 then stop
  endif
  
  ;; print fields
  IF ptr_valid(self.dark_dir) THEN BEGIN
    if(self.dopolcal) then print, 'red::initialize : polcal_dir = '+ self.polcal_dir
    nn = n_elements(*self.dark_dir)
    IF(nn EQ 1) THEN BEGIN
      print, 'red::initialize : dark_dir = '+ (*self.dark_dir)[0]
    ENDIF ELSE BEGIN
      print, 'red::initialize : dark_dirs :'
      FOR k = 0, nn-1 DO print, string(k, format='(I5)') + ' -> ' + (*self.dark_dir)[k]
    ENDELSE
  ENDIF
  IF ptr_valid(self.flat_dir) THEN BEGIN
    nn = n_elements(*self.flat_dir)
    IF(nn EQ 1) THEN BEGIN
      print, 'red::initialize : flat_dir = '+ (*self.flat_dir)[0]
    ENDIF ELSE BEGIN
      print, 'red::initialize : flat_dirs :'
      FOR k = 0, nn-1 DO print, string(k, format='(I5)') + ' -> ' + (*self.flat_dir)[k]
    ENDELSE
  ENDIF
  IF ptr_valid(self.cameras) THEN BEGIN
    nn = n_elements(*self.cameras)
    IF(nn EQ 1) THEN BEGIN
      print, 'red::initialize : camera = '+ (*self.cameras)[0]
    ENDIF ELSE BEGIN
      print, 'red::initialize : cameras :'
      FOR k = 0, nn-1 DO print, string(k, format='(I5)') + ' -> ' + (*self.cameras)[k]
    ENDELSE
  ENDIF
  if ptr_valid(self.pinh_dirs) then begin
    nn = n_elements(*self.pinh_dirs)
    if(nn eq 1) then begin
      print, 'red::initialize : pinh_dir = '+ (*self.pinh_dirs)[0]
    endif else begin
      print, 'red::initialize : pinh_dirs :'
      for k = 0, nn-1 do print, string(k, format='(I5)') + ' -> ' + (*self.pinh_dirs)[k]
    endelse
  endif
  ;;if(self.doflat) then print, 'red::initialize : flat_dir = '+ strjoin(*self.flat_dir, '  ')
  if ptr_valid(self.data_dirs) then begin
    nn = n_elements(*self.data_dirs)
    if(nn eq 1) then begin
      print, 'red::initialize : data_dir = '+ (*self.data_dirs)[0]
    endif else begin
      print, 'red::initialize : data_dirs :'
      for k = 0, nn-1 do print, string(k, format='(I5)') + ' -> ' + (*self.data_dirs)[k]
    endelse
  endif
  if(self.dodescatter) then print, 'red::initialize : descatter_dir = '+ self.descatter_dir
  if(self.filetype) then print, 'red::initialize : filetype = '+ self.filetype
  print, 'red::initialize : out_dir = '+ self.out_dir



  ;; Versions and libraries info ----------------------------------------------------

  paths = strsplit(!path,":",/extract) 
  git_describe_command = 'git describe --always --abbrev=12 --long --dirty=\ \(Modified\)'
  git_count_command = 'git rev-list HEAD --count'
  git_diff_command = 'git diff HEAD'
  git_status_command = 'git status'

  
  ;; Pipeline version
  srcdir = file_dirname( routine_filepath("red::initialize"), /mark )
  spawn, 'cd '+srcdir+'; ' + git_describe_command, pipeline_gitoutput
  pipeline_gitoutput = red_strreplace(pipeline_gitoutput, 'release/', '')
  spawn, 'cd '+srcdir+'; ' + git_status_command, pipeline_status 
  if strmatch(pipeline_gitoutput, '*(Modified)') then $
     self.version_problems += 'The pipeline is modified. '
  self.version_pipeline = strjoin((strsplit(pipeline_gitoutput, '-', /extract))[0:1], '-')

  ;; Pipeline branch
  spawn, 'cd '+srcdir+'; ' + git_status_command, pipeline_status
  self.version_pipeline_branch = (strsplit(pipeline_status[0],/extract))[-1]
  if strmatch(self.version_pipeline_branch,'*_dev') then  self.version_problems += 'Running the development branch ' $
     + self.version_pipeline_branch
  
  ;; Redux dlm version. We require that the ANA and MOMFBD dlms are
  ;; part of the rdx dlm and the same version.
  help,/dlm,'rdx', output = rdx_dlm_version
  dlmpos = strpos(rdx_dlm_version[1], 'release/')
  self.version_reduxdlm = strjoin((strsplit(strmid(rdx_dlm_version[1],dlmpos+8), '-', /extract))[0:1], '-')
  help,/dlm,'ana' , output = ana_dlm_version
  dlmpos = strpos(ana_dlm_version[1], 'release/')
  ana_dlm_version = strjoin((strsplit(strmid(ana_dlm_version[1],dlmpos+8), '-', /extract))[0:1], '-')
  help,/dlm,'momfbd', output = momfbd_dlm_version
  dlmpos = strpos(momfbd_dlm_version[1], 'release/')
  momfbd_dlm_version = strjoin((strsplit(strmid(momfbd_dlm_version[1],dlmpos+8), '-', /extract))[0:1], '-')

  if ana_dlm_version ne self.version_reduxdlm then self.version_problems += 'ANA DLM not identical to redux DLM'
  if momfbd_dlm_version ne self.version_reduxdlm then self.version_problems += 'MOMFBD DLM not identical to redux DLM'
  
  rdx_dlm_required_versionstr = '1.0.0-1'
  rdx_dlm_required_version = fix(strsplit(rdx_dlm_required_versionstr, '-.', /EXTRACT))
  rdx, version=rdx_dlm_version
  rdx_dlm_version = fix(strsplit(rdx_dlm_version, '-.', /EXTRACT))
  n_elem = min([n_elements(rdx_dlm_required_version),n_elements(rdx_dlm_version)])
  for i=0,n_elem-1 do begin
    if rdx_dlm_version[i] gt rdx_dlm_required_version[i] then break
    if rdx_dlm_version[i] lt rdx_dlm_required_version[i] then begin
      self.version_problems += 'The RDX DLM should be updated to version >= ' + rdx_dlm_required_versionstr
      break
    endif
  endfor

  ;; Coyote library version
  coyotepaths = paths(where(strmatch(paths,'*coyote'), Ncoyote))
  case Ncoyote of
    0: begin
      print, 'The Coyote library does not seem to be in your !path.'
      stop
    end
    1: begin
      ;coyotedir = file_dirname( filepath(root_dir = coyotepaths[0], "cgcolor"), /mark )
      spawn, 'cd '+coyotepaths[0]+'; ' + git_count_command, coyote_gitoutput
      coyote_gitoutput = red_strreplace(coyote_gitoutput, 'release/', '')
      ;; if strmatch(coyote_gitoutput, '(Modified)') then $
      ;;   self.version_problems += 'The Coyote library is modified. '
      ;; self.version_coyote = strjoin((strsplit(coyote_gitoutput, '-', /extract))[0:1], '-')
      self.version_coyote = coyote_gitoutput
      spawn, 'cd '+coyotepaths[0]+'; ' + git_diff_command, coyote_gitoutput
      if size(coyote_gitoutput, /n_dim) gt 0 then self.version_problems $
         += 'The Coyote library is not the latest version. '
    end
    else: begin
      self.version_coyote = 'Undefined'
      self.version_problems += 'Multiple Coyote directories. '
    end
  endcase


  ;; IDLastro library version
  idlastropaths = paths(where(strmatch(paths, '*IDLAstro/pro'), Nwhere))
  case Nwhere of
    0: begin
      print, 'The IDLAstro library does not seem to be in your !path.'
      stop
    end
    1: begin
      spawn, 'cd '+idlastropaths[0]+'; ' + git_count_command, idlastro_gitoutput
      idlastro_gitoutput = red_strreplace(idlastro_gitoutput, 'release/', '')
      ;; if strmatch(idlastro_gitoutput, '(Modified)') then $
      ;;  self.version_problems += 'The IDLAstro library is modified. '
      ;; self.version_idlastro = strjoin((strsplit(idlastro_gitoutput, '-', /extract))[0:1], '-')
      self.version_idlastro = idlastro_gitoutput
      spawn, 'cd '+idlastropaths[0]+'; ' + git_diff_command, idlastro_gitoutput
      if size(idlastro_gitoutput, /n_dim) gt 0 then self.version_problems $
         += 'The IDLAstro library is not the latest version. '
    end
    else: begin
      self.version_idlastro = 'Undefined.'
      self.version_problems += 'Multiple IDLAstro directories. '
    end
  endcase


  ;; mpfit version  
  mpfitpaths = paths(where(strmatch(paths, '*mpfit'), Nwhere))
  case Nwhere of
    0: begin
      print, 'The mpfit library does not seem to be in your !path.'
      stop
    end
    1: begin
      red_mpfit_version, mpfitpaths, local_version = local_version
      self.version_mpfit = local_version 
      ;; mpfit is not under version control so we will use the time of
      ;; the latest chenge, as defined by the $Id string.
;      spawn, 'grep "\$Id" '+mpfitpaths[0]+'/*.pro', mpfit_spawnoutput
;      timestamps = STREGEX(mpfit_spawnoutput,'[0-9][0-9][0-9][0-9]/[0-9][0-9]/[0-9][0-9] ' $
;                           + '[0-9][0-9]:[0-9][0-9]:[0-9][0-9]',/EXTRACT)
;      for i = 0, n_elements(timestamps)-1 do timestamps[i] = red_strreplace(timestamps[i], '/', '-', n = 2)
;      for i = 0, n_elements(timestamps)-1 do timestamps[i] = red_strreplace(timestamps[i], ' ', 'T')
;      self.version_mpfit = (timestamps(sort(timestamps)))[n_elements(timestamps)-1]
    end
    else: begin
      self.version_mpfit = 'Undefined'
      self.version_problems += 'Multiple mpfit directories. '
    end
  endcase

  ;; Report problems
  if strlen(self.version_problems) gt 0 then begin
    print
    print, 'Problem(s) with your installation:'
    print, self.version_problems
    print

    if Ncoyote gt 1 then begin
      print, 'You seem to have multiple copies of the coyote library in your path.'
      print, 'Here is the output of cgFindCoyoteFiles:'
      print
      cgFindCoyoteFiles
    endif

    print
    if self.developer_mode then begin
      print, 'You are running ' + strlowcase(tag_names(self,/STRUCTURE_NAME)) + 'red in developer mode, your output'
      print, 'will be marked as not conforming to the SOLARNET standard.'
      print
    endif else begin
      print, 'You can proceed with your processing but your output will be marked as'
      print, 'not conforming to the SOLARNET standard.'
      print
      print, 'Just restart with the /develop keyword set:'
      print
      print, '     a = ' + strlowcase(tag_names(self,/STRUCTURE_NAME)) + "red('config.txt', /develop)"
      print
      retall
    endelse
  endif

  cgWindow_SetDefs, PS_Decomposed=1
  
  return
end
