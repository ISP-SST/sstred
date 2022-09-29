; docformat = 'rst'

;+
;  Select (a subset) of files and return their filenames and state information
; 
; :Categories:
;
;    CHROMIS pipeline
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
;   2016-06-01 : THI. Allow filtering on camera tag.
;                Pass keyword strip_settings to extractstates.
;
;   2016-06-02 : MGL. Remove some keywords to extractstates.
;
;   2016-06-09 : MGL. Bugfix: states.pref --> states.prefilter. 
;
;   2016-08-23 : THI. Rename camtag to detector and channel to camera,
;                so the names match those of the corresponding SolarNet
;                keywords.
;
;   2016-08-31 : MGL. Update file_template for data files not starting
;                with the detector name.
; 
;   2016-09-23 : THI. Add selection by scannumber
; 
;   2016-10-13 : THI. Add selection by framenumbers and option strip_wb to
;                be forwarded to extractstates
;
;   2016-10-27 : MGL. New keywords count, complement, and ncomplement.
;
;   2017-12-20 : MGL. Always do extractstates, if the files are the
;                same this is a fast operation now because of caching.
;
;   2017-12-21 : MGL. Do extractstates only if needed.
; 
;   2018-05-30 : MGL. Add selection by fpi_states.
; 
;-
pro chromis::selectfiles, cam = cam $ 
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
                          , prefilter = prefilter $
                          , scan = scan $
                          , selected = selected $
                          , states = states $
                          , strip_settings = strip_settings $
                          , subdir = subdir $
                          , timestamps = timestamps $
                          , ustat = ustat

  compile_opt idl2
  
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

    file_template = subdir + '/*' + strtrim(detector, 2)+ '*'

    if( n_elements(dirs) gt 0 ) then dirs = [dirs] $ ; ensure it's an array, even with 1 element
    else begin 
      if( keyword_set(dark) && ptr_valid(self.dark_dir) ) then dirs = *self.dark_dir $
      else if( keyword_set(flat) && ptr_valid(self.flat_dir) ) then dirs = *self.flat_dir
    endelse
    
    path_spec = dirs + '/' + file_template
    
    files = file_search(path_spec)
    files = files[where( strpos(files, '.lcd.') LT 0, nf)]
    files = strtrim(files,2)
    idx = where( files ne '' )
    if min(idx) ge 0 then files = files[ idx ] $
    else return
;;    files = red_sortfiles(files) ; slow
    force = 1                   ; we have new files, force extractstates
  endif
  
  if( n_elements(files) eq 0 || files[0] eq '' ) then begin
    return
  endif
  
  if( n_elements(files) eq 1 ) then files = [files]

  ;; Get new states if necessary
  if keyword_set(force) $                             ; If we choose to
     || n_elements(files) ne n_elements(states) $     ; If numbers don't agree
     || ~min(file_same(files,states.filename)) then $ ; If the actual files are not the same
        self->extractstates, files, states, strip_settings = strip_settings, cam = cam

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
                                ; Prefilter is not allowed to select empty strings.
    idx = where(prefilter ne '',count)
    if count ne 0 then begin
      prefilter = prefilter[idx] 
      
      selected = [states.skip] * 0
      tpref = [prefilter]       ; make sure it's an array
      for ip = 0, Npref-1 do begin
        undefine, pos
        for is = 0, n_elements(states)-1 do begin
          if( self->match_prefilters(states[is].prefilter, tpref) ) then red_append, pos, is
        endfor
        if( n_elements(pos) gt 0 ) then begin
          selected[pos] = 1
        endif
      endfor
      pos = where(selected lt 1,count)
      if( count ne 0 ) then states[pos].skip = 1
    endif
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
    selected = [states.skip] * 0
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
  
                                ; if keyword selected is not present, return selected subsets as new files/states
  if count ne 0 then begin
    states = states[selected]
    files = states.filename
  endif else begin              ; return empty files/states
    dummy = temporary(states)
    dummy = temporary(files)
  endelse

end

a = chromisred(/dev)

files = file_search('/scratch/mats/2016.09.19/CHROMIS-jaime_recipe/data/09:28:36/Chromis-N/*', count = Nfiles1)
print, Nfiles1
a -> selectfiles, files = files, states = states, scan = [3, 4, 5], sel = sel

print, states[sel].tun_wavelength

files = file_search('/scratch/mats/2016.09.19/CHROMIS-jaime_recipe/data/09:28:36/Chromis-W/*', count = Nfiles2)
print, nfiles2
a -> selectfiles, files = files, states = states, scan = [3, 4, 5], sel = sel
print, states[sel].tun_wavelength


end
