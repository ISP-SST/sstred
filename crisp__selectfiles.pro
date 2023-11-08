; docformat = 'rst'

;+
;  Select (a subset) of files and return their filenames and state information
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;     Tomas Hillberg, ISP
; 
; 
; :Params:
; 
; 
; :Keywords:
; 
;     cam : in, optional, type=string
;   
;         Name of the camera channel. Ex. 'Crisp-W'
;   
;     dirs : in, optional, type=strarr
;
;         Locations where to find files.
; 
;     files : in, optional, type=strarr
;
;         File-list. No file search will be done
; 
;     states : in, optional, type=strarr
;   
;         State information. Extractstates will only be called if
;         files were updated.
;   
;     prefilter : in, optional, type=strarr
;   
;         only return files matching at least one of these prefilters
;   
;     ustat : in, optional, type=strarr
;   
;         only return files matching at least one of these states
;   
;     flat : in, optional, type=boolean
;
;        If no dirs were specified, use self.flat_dir
; 
;     dark : in, optional, type=boolean
;
;        If no dirs were specified, use self.dark_dir
; 
;     nremove : in, optional, type=int
;
;        Skip this many frames after each state-change
; 
;     force : in, optional, type=boolean
;
;        Re-populate files/states.
;
;     scan : in, optional, type=intarr
;
;        Only return files matching these scan numbers.
; 
;     polcal : in, optional, type=boolean
; 
;        Set this to add polcal-specific items in the states, qw and
;        lp. 
; 
;     selected : out, optional, type=intarr
;
;         Return the indices for the selection. If this is not specified,
;         the files/states arrays will be over-written to only contain the
;         selection.
;
;     count : out, optional, type=integer
;
;         The number of selected files.
;
;     complement : out, optional, type=intarr
;
;         The complement of selected.
;
;     ncomplement : out, optional, type=integer
;
;         The number of non-selected files.
; 
; :History:
; 
;   2016-05-19 : First version.
; 
;   2016-08-23 : THI. Rename camtag to detector and channel to camera,
;                so the names match those of the corresponding SolarNet
;                keywords.
; 
;   2017-06-28 : MGL. Adapt to new version of crisp::selectfiles.
;
;   2017-06-29 : MGL. Bugfix: states.pref --> states.prefilter.
;
;   2017-07-06 : MGL. New keyword polcal. Do not sort files.
;
;   2018-04-18 : MGL. New keywords count, complement, and ncomplement.
; 
;   2018-05-25 : MGL. Add selection by scannumber, framenumbers, and
;                fpi_states. 
; 
;-
pro crisp::selectfiles, cam = cam $
                        , complement = complement $
                        , count = count $
                        , dark = dark $
                        , dirs = dirs $
                        , files = files $
                        , flat = flat $
                        , force = force $
                        , fpi_states = fpi_states $
                        , framenumbers = framenumbers $
                        , ncomplement = ncomplement $
                        , nremove = nremove $
                        , polcal = polcal $
                        , prefilter = prefilter $
                        , scan = scan $
                        , selected = selected $
                        , states = states $
                        , strip_settings = strip_settings $
                        , subdir = subdir $
                        , timestamps = timestamps $
                        , ustat = ustat 
  
  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)                                      
  
  ;; Unless we select any
  count = 0L                  
  ncomplement = n_elements(files)

  if( keyword_set(force) || n_elements(files) eq 0 ) then begin
    
    if( n_elements(cam) ne 1 ) then begin
      print,inam+': Only a single cam supported.'
      return
    endif
    detector = self->RED::getdetector(cam)
    
    if( n_elements(subdir) ne 1 ) then subdir = cam

    file_template = subdir + '/' + strtrim(detector, 2)+ '*'
    if n_elements(scan) eq 1 then file_template += string(scan, format='(I05)') + '*'
    
    if( n_elements(dirs) gt 0 ) then dirs = [dirs] $ ; ensure it's an array, even with 1 element
    else begin 
      if( keyword_set(dark) && ptr_valid(self.dark_dir) ) then dirs = *self.dark_dir $
      else if( keyword_set(flat) && ptr_valid(self.flat_dir) ) then dirs = *self.flat_dir
    endelse
    
    path_spec = dirs + '/' + file_template

    files = file_search(path_spec)
    files = files(where( strpos(files, '.lcd.') LT 0, nf) )
;    files = red_sortfiles(files) ; Slow! Move sorting to when we have a single state?
    force = 1                   ; we have new files, force extractstates
  endif
  
  if( n_elements(files) eq 0 || files[0] eq '' ) then begin
    return
  endif
  
  if( n_elements(files) eq 1 ) then files = [files]

  if( n_elements(force) gt 0 || n_elements(states) eq 0 ) then begin
;    self->extractstates, files, states, /basename, /cam, /prefilter, /fullstate
    self->extractstates, files, states, strip_settings = strip_settings, polcal = polcal, cam = cam
                                ; TODO: this is a really ugly way to drop the WB states, think of something better
;    wb_cams = (strmatch( states.detector, self->getdetector('Crisp-W')) $
;               or strmatch( states.detector, self->getdetector('Crisp-D')))
;    pos = where(wb_cams gt 0)
;                                ; replace NB state-info with prefilter for the WB cameras
;    if( min(pos) ge 0 ) then states[pos].fullstate = states[pos].prefilter
  endif

  if( n_elements(states) eq 1 ) then states = [states]
  
  states.skip *= 0              ; always clear selection, so repeated calls with the same files/states are possible
  if( keyword_set(nremove) ) then self->skip, states, nremove

  Ncam = n_elements(cam)
  if( Ncam gt 0 ) then begin
    selected = [states.skip] * 0
    camt = [cam]                ; make sure it's an array
    for ip = 0, Ncam-1 do begin
      pos = where(states.camera eq camt[ip],count)
      if( count ne 0 ) then begin
        selected[pos] = 1
      endif
    endfor
    pos = where(selected lt 1,count)
    if( count ne 0 ) then states[pos].skip = 1
  endif
  
  Npref = n_elements(prefilter)
  if( Npref gt 0 ) then begin
    selected = states.skip * 0
    prefilter = [prefilter]     ; make sure it's an array
    for ip = 0, Npref-1 do begin
      pos = where(states.prefilter eq prefilter[ip])
      if( min(pos) ge 0 ) then begin
        selected[pos] = 1
      endif
                                ; stop
    endfor                      ; ip
    states[where(selected lt 1)].skip = 1
  endif
  
  Nframes = n_elements(framenumbers)
  if( Nframes gt 0 ) then begin
    selected = [states.skip] * 0
    tframes = [framenumbers]    ; make sure it's an array
    for is = 0, Nframes-1 do begin
      pos = where(states.framenumber eq tframes[is],count)
      if( count ne 0 ) then selected[pos] = 1
    endfor
    pos = where(selected lt 1,count)
    if( count ne 0 ) then states[pos].skip = 1
  endif
  
  Nscan = n_elements(scan)
  if( Nscan gt 0 ) then begin
    selected = [states.skip] * 0
    tscans = [scan]             ; make sure it's an array
    for is = 0, Nscan-1 do begin
      pos = where(states.scannumber eq tscans[is],count)
      if( count ne 0 ) then selected[pos] = 1
    endfor
    pos = where(selected lt 1,count)
    if( count ne 0 ) then states[pos].skip = 1
  endif
  
  Nstates = n_elements(ustat)
  if( Nstates gt 0 ) then begin
    selected = states.skip * 0
    tstates = [ustat]           ; make sure it's an array
    for ip = 0, Nstates-1 do begin
      pos = where(states.fullstate eq tstates[ip],count)
      if( count ne 0 ) then selected[pos] = 1
    endfor
    pos = where(selected lt 1,count)
    if( count ne 0 ) then states[pos].skip = 1
  endif

  Nfpi = n_elements(fpi_states)
  if( Nfpi gt 0 ) then begin
    selected = states.skip * 0
    tfpi = [fpi_states]         ; make sure it's an array
    for ip = 0, Nfpi-1 do begin
      pos = where(states.fpi_state eq tfpi[ip])
      if( min(pos) ge 0 ) then selected[pos] = 1
    endfor
    states[where(selected lt 1)].skip = 1
  endif

  Ntimestamps = n_elements(timestamps)
  if( Ntimestamps gt 0 ) then begin
    selected = states.skip * 0
    ts = [timestamps]           ; make sure it's an array
    for ip = 0, Ntimestamps-1 do begin
      pos = where(strmatch(states.filename,'*'+ts[ip]+'*'),count)
      if( count ne 0 ) then selected[pos] = 1
    endfor
    pos = where(selected lt 1,count)
    if( count ne 0 ) then states[pos].skip = 1
  endif
  
  selected = where(states.skip lt 1, count $
                   , complement = complement, Ncomplement = Ncomplement)

  if arg_present(selected) then begin
    if count eq 0 then undefine,selected ; don't return -1
    return
  endif
  
  ;; if keyword selected is not present, return selected subsets as new files/states
  if count ne 0 then begin
    states = states[selected]
    files = states.filename
  endif else begin              ; return empty files/states
    dummy = temporary(states)
    dummy = temporary(files)
  endelse

end
