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
; :Returns:
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
;-
pro crisp::selectfiles, cam = cam $
                        , dirs = dirs $
                        , files = files $
                        , states = states $
                        , prefilter = prefilter $
                        , ustat = ustat $
                        , flat = flat $
                        , dark = dark $
                        , nremove = nremove $
                        , force = force $
                        , selected = selected $
                        , strip_settings = strip_settings $
                        , polcal = polcal

  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])
 
  if( n_elements(force) gt 0 || n_elements(files) eq 0 ) then begin

    if( n_elements(cam) ne 1 ) then begin
      print,inam+': Only a single cam supported.'
      return
    endif
    detector = self->RED::getdetector(cam)
    
    file_template = cam + '/' + strtrim(detector, 2)+ '*'
    
    if( n_elements(dirs) gt 0 ) then dirs = [dirs] $ ; ensure it's an array, even with 1 element
    else begin 
      if( keyword_set(dark) && ptr_valid(self.dark_dir) ) then dirs = *self.dark_dir $
      else if( keyword_set(flat) && ptr_valid(self.flat_dir) ) then dirs = *self.flat_dir
    endelse
    
    path_spec = dirs + '/' + file_template

    print, 212
    files = file_search(path_spec)
    files = files(where( strpos(files, '.lcd.') LT 0, nf) )
;    files = red_sortfiles(files) ; Slow! Move sorting to when we have a single state?
    force = 1                   ; we have new files, force extractstates
  endif
  
  if( n_elements(files) eq 0 || files[0] eq '' ) then begin
    print, 111
    return
  endif

  if( n_elements(force) gt 0 || n_elements(states) eq 0 ) then begin
;    self->extractstates, files, states, /basename, /cam, /prefilter, /fullstate
    self->extractstates, files, states, strip_wb=strip_wb, strip_settings = strip_settings, polcal = polcal
                                ; TODO: this is a really ugly way to drop the WB states, think of something better
;    wb_cams = (strmatch( states.detector, self->getdetector('Crisp-W')) $
;               or strmatch( states.detector, self->getdetector('Crisp-D')))
;    pos = where(wb_cams gt 0)
;                                ; replace NB state-info with prefilter for the WB cameras
;    if( min(pos) ge 0 ) then states[pos].fullstate = states[pos].prefilter
  endif

  states.skip *= 0              ; always clear selection, so repeated calls with the same files/states are possible
  if( keyword_set(nremove) ) then self->skip, states, nremove

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
  
  Nstates = n_elements(ustat)
  if( Nstates gt 0 ) then begin
    selected = states.skip * 0
    tstates = [ustat]           ; make sure it's an array
    for ip = 0, Nstates-1 do begin
      pos = where(states.fullstate eq tstates[ip])
      if( min(pos) ge 0 ) then selected[pos] = 1
    endfor
    states[where(selected lt 1)].skip = 1
  endif

  selected = where( states.skip lt 1 )
  
  if arg_present(selected) then return
  
                                ; if keyword selected is not present, return selected subsets as new files/states
  if( min(selected) ge 0 ) then begin
    states = states[selected]
    files = states.filename
  endif else begin              ; return empty files/states
    dummy = temporary(states)
    dummy = temporary(files)
  endelse

end
