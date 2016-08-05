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
;   2016-06-01 : THI. Allow filtering on camera tag.
;                Pass keyword strip_settings to extractstates.
;
;   2016-06-02 : MGL. Remove some keywords to extractstates.
;
;   2016-06-09 : MGL. Bugfix: states.pref --> states.prefilter. 
; 
;-
pro chromis::selectfiles, cam = cam $
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
                      , strip_settings = strip_settings

    inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])
    
    if( n_elements(force) gt 0 || n_elements(files) eq 0 ) then begin

        if( n_elements(cam) ne 1 ) then begin
            print,inam+': Only a single cam supported.'
            return
        endif
        camtag = self->RED::getcamtag(cam)
        
        file_template = cam + '/' + camtag + '*'
        
        if( n_elements(dirs) gt 0 ) then dirs = [dirs] $    ; ensure it's an array, even with 1 element
        else begin 
            if( keyword_set(dark) && ptr_valid(self.dark_dir) ) then dirs = *self.dark_dir $
            else if( keyword_set(flat) && ptr_valid(self.flat_dir) ) then dirs = *self.flat_dir
        endelse
        
        path_spec = dirs + '/' + file_template
    
        files = file_search(path_spec)
        files = files(where( strpos(files, '.lcd.') LT 0, nf) )
        files = red_sortfiles(files)
        force = 1   ; we have new files, force extractstates
    endif
    
    if( n_elements(files) eq 0 || files[0] eq '' ) then begin
        return
    endif
    
    if( n_elements(files) eq 1 ) then files = [files]
    
    if( n_elements(force) gt 0 || n_elements(states) eq 0 ) then begin
        self->extractstates, files, states, /strip_wb, strip_settings = strip_settings
    endif

    if( n_elements(states) eq 1 ) then states = [states]
    
    states.skip *= 0    ; always clear selection, so repeated calls with the same files/states are possible
    if( keyword_set(nremove) ) then self->skip, states, nremove

    Ncam = n_elements(cam)
    if( Ncam gt 0 ) then begin
        selected = [states.skip] * 0
        camt = [cam]    ; make sure it's an array
        for ip = 0, Ncam-1 do begin
            pos = where(states.channel eq camt[ip])
            if( min(pos) ge 0 ) then begin
                selected[pos] = 1
            endif
        endfor
        pos = where(selected lt 1)
        if( min(pos) ge 0 ) then states[pos].skip = 1
    endif
    
    Npref = n_elements(prefilter)
    if( Npref gt 0 ) then begin
        selected = [states.skip] * 0
        tpref = [prefilter]    ; make sure it's an array
        for ip = 0, Npref-1 do begin
            pos = where(states.prefilter eq tpref[ip])
            if( min(pos) ge 0 ) then begin
                selected[pos] = 1
            endif
        endfor
        pos = where(selected lt 1)
        if( min(pos) ge 0 ) then states[pos].skip = 1
    endif
    
    Nstates = n_elements(ustat)
    if( Nstates gt 0 ) then begin
        selected = [states.skip] * 0
        tstates = [ustat]    ; make sure it's an array
        for ip = 0, Nstates-1 do begin
            pos = where(states.fullstate eq tstates[ip])
            if( min(pos) ge 0 ) then selected[pos] = 1
        endfor
        pos = where(selected lt 1)
        if( min(pos) ge 0 ) then states[pos].skip = 1
    endif

    selected = where( states.skip lt 1 )
    
    if arg_present(selected) then return
    
    ; if keyword selected is not present, return selected subsets as new files/states
    if( min(selected) ge 0 ) then begin
        states = states[selected]
        files = states.filename
    endif else begin    ; return empty files/states
        dummy = temporary(states)
        dummy = temporary(files)
    endelse

end