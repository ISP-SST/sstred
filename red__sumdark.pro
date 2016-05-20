; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :author:
; 
; 
; 
; 
; :returns:
; 
; 
; :Params:
; 
; 
; :Keywords:
; 
;    overwrite  : 
;   
;   
;   
;    check  : 
;   
;   
;    cams : in, optional, type=strarr
; 
;      A list of cameras (or rather camera subdirs).
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2013-08-23 : MGL. Added support for logging.
; 
;   2013-08-27 : MGL. Let the subprogram find out its own name. Added
;                self.done to the selfinfo.
; 
;   2013-09-19 : MGL. Make the sum long integer.
; 
;   2014-03-21 : MGL. Allow for multiple dark_dir. Make it work for
;                blue cameras. New keyword cams.
;
;   2014-04-07 : THI. Bugfix: look for darks in dark_dir.
;
;   2016-05-17 : THI. Use cam_channels by default, added keyword dirs
;                to use specific dark-folder. Re-write to sum files from
;                multiple folders at once.
;
;-
pro red::sumdark, overwrite = overwrite, $
                  check = check, $
                  cams = cams, $
                  dirs = dirs

    ;; Defaults
    if n_elements(overwrite) eq 0 then overwrite = 0

    if n_elements(check) eq 0 then check = 0

    if n_elements(dirs) gt 0 then dirs = [dirs] $
    else if ptr_valid(self.dark_dir) then dirs = *self.dark_dir

    if n_elements(cams) gt 0 then cams = [cams] $
    else if ptr_valid(self.cam_channels) then cams = *self.cam_channels

    ;; Name of this method
    inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

    ;; Logging
    help, /obj, self, output = selfinfo1
    help, /struct, self.done, output = selfinfo2 
    red_writelog, selfinfo = [selfinfo1, selfinfo2]

    Ncams = n_elements(cams)
    if( Ncams eq 0) then begin
        print, inam+' : ERROR : undefined cams (and cam_channels)'
        return
    endif

    Ndirs = n_elements(dirs)
    if( Ndirs eq 0) then begin
        print, inam+' : ERROR : no dark directories defined'
        return
    endif else begin
        if Ndirs gt 1 then dirstr = '['+ strjoin(dirs,';') + ']' $
        else dirstr = dirs[0]
    endelse

    outdir = self.out_dir+'/darks/'
    darkname = self->getdark('', summed_name=sdarkname)
    file_mkdir, file_dirname(darkname)
    file_mkdir, file_dirname(sdarkname)

    for ic = 0L, Ncams-1 do begin

        cam = cams[ic]
        camtag = self->getcamtag( cam )

        darkname = self->getdark(cam, summed_name=sdarkname)
        if( ~keyword_set(overwrite) && file_test(darkname) && file_test(sflatname)) then begin
            print, inam+' : file exists: ' + darkname + ' , skipping! (run sumdark, /overwrite to recreate)'
            continue
        endif

        self->selectfiles, cam=cam, dirs=dirs, $
                         files=files, states=states, /force

        nf = n_elements(files)
        if( nf eq 0 || files[0] eq '') then begin
            print, inam+' : '+cam+': no files found in: '+dirstr
            print, inam+' : '+cam+': skipping camera!'
            continue
        endif else begin
            print, inam+' : Found '+red_stri(nf)+' files in: '+ dirstr + '/' + cam + '/'
        endelse

        ; TODO: get/modify/put header in a generic way
        head = fzhead(files[0])
        if rdx_hasopencv() then begin
            dark = rdx_sumfiles(files, check = check, summed = darksum, nsum=nsum, verbose=2)
        endif else begin
            dark = red_sumfiles(files, check = check, summed = darksum, nsum=nsum)
        endelse
        ;; Normalize dark
        ;;dark = float(tmp / count)

        outheader = red_dark_h2h(files[nf-1], head, nsum)

        ; TODO: output should be written through class-specific methods
        print, inam+' : saving ', darkname
        fzwrite, dark, darkname, 'mean dark'

        print, inam+' : saving ', sdarkname
        fzwrite, long(darksum), sdarkname, outheader

    endfor                        ; (ic loop)

    self.done.sumdark = 1B

    return

end
