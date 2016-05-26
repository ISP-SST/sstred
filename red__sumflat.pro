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
; 
; :Keywords:
; 
;    overwrite  : 
;   
;   
;   
;    ustat  : 
;   
;   
;   
;    old  : 
;   
;   
;   
;    remove  : 
;   
;   
;   
;    check : in, optional, type=boolean
;
;
;    sum_in_idl : in, optional, type=boolean
;
;      Bypass rdx_sumfiles.
;   
;   
;   
; 
; 
; :History:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2013-08-27 : MGL. Added support for logging. Let the subprogram
;                find out its own name. Added self.done to the
;                selfinfo.
; 
;   2013-08-29 : MGL. Remove *lcd* files from file list.
; 
;   2013-10-30 : MGL. Get correct prefilter state for WB flats in the
;                presence of the focus file name field. 
; 
;   2013-12-10 : PS  Adapt for multiple flat directories; add lim
;                keyword (passthrough to red_sumfiles)
;
;   2013-12-13 : PS  Only store raw sums if requested.  Save
;                unnormalized sum and add correct number of summed
;                files in header
;
;   2014-01-23 : MGL. Use red_extractstates instead of red_getstates
;                and local extraction of info from file names.
;
;   2014-10-06 : JdlCR. Prefilter keyword, please DO NOT REMOVE!
;
;   2016-05-19 : THI: Re-write to use overloaded class methods.
; 
;   2016-05-25 : MGL. New keyword sum_in_idl. Get darks from files
;                with state info in the name.
;
;-
pro red::sumflat, overwrite = overwrite, ustat = ustat, old = old, $
                  remove = remove, cams = cams, check = check, lim = lim, $
                  store_rawsum = store_rawsum, prefilter = prefilter, $
                  dirs = dirs, $
                  sum_in_idl = sum_in_idl

    ;; Defaults
    if( n_elements(overwrite) eq 0 ) then overwrite = 0

    if( n_elements(check) eq 0 ) then check = 0

    if( n_elements(dirs) gt 0 ) then dirs = [dirs] $
    else if ptr_valid(self.flat_dir) then dirs = *self.flat_dir

    if( n_elements(cams) gt 0 ) then cams = [cams] $
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
        print, inam+' : ERROR : no flat directories defined'
        return
    endif else begin
        if Ndirs gt 1 then dirstr = '['+ strjoin(dirs,';') + ']' $
        else dirstr = dirs[0]
    endelse


    flatname = self->getflat( '', state='', summed_name=sflatname )
    file_mkdir, file_dirname(flatname)
    if keyword_set(store_rawsum) then file_mkdir, file_dirname(sflatname)

    ;; cameras
    for ic = 0L, Ncams-1 do begin

        cam = cams[ic]
        camtag = self->getcamtag( cam )

        self->selectfiles, cam=cam, dirs=dirs, prefilter=prefilter, ustat=ustat, $
                         files=files, states=states, nremove=remove, /force

        nf = n_elements(files)
        if( nf eq 0 || files[0] eq '') then begin
            print, inam+' : '+cam+': no files found in: '+dirstr
            print, inam+' : '+cam+': skipping camera!'
            continue
        endif else begin
            print, inam+' : Found '+red_stri(nf)+' files in: '+ dirstr + '/' + cam + '/'
        endelse

        gain = states.gain[uniq(states.gain, sort(states.gain))]
        exposure = states.exposure[uniq(states.exposure, sort(states.exposure))]
        state = string(exposure*1000, format = '(f4.2)')+'ms_G'+string(gain, format = '(f05.2)')

        dname = self->getdark(cam, state = state, data=dd)
        if( n_elements(dd) eq 0 ) then begin
            print, inam+' : no dark found for camera ', cam
            continue
        endif
stop
        state_list = [states[uniq(states.fullstate, sort(states.fullstate))].fullstate]

        ns = n_elements(state_list)

        ;; Loop over states and sum
        for ss = 0L, ns - 1 do begin

            flatname = self->getflat( cam, state=state_list[ss], summed_name=sflatname )

            ;; If file does not exist, do sum!
            if( ~keyword_set(overwrite) && file_test(flatname) ) then begin
                if (~keyword_set(store_rawsum) || file_test(sflatname)) then begin ; only skip if rawsum also exists
                   print, inam+' : file exists: ' + flatname + ' , skipping! (run sumflat, /overwrite to recreate)'
                   continue
                endif
            endif

            self->selectfiles, prefilter=prefilter, ustat=state_list[ss], $
                           files=files, states=states, selected=sel

            if( min(sel) lt 0 ) then begin
                print, inam+' : '+cam+': no files found for state: '+state_list[ss]
                continue
            endif

            print, inam+' : adding flats for state -> ' + state_list[ss]
            if(keyword_set(check)) then begin
                openw, lun, self.out_dir + '/flats/'+camtag+'.'+state_list[ss]+'.discarded.txt', width = 500, /get_lun
            endif

            ;; Sum files
            if keyword_set(check) then begin
                ;;; if summing from two directories, same frame
                ;;; numbers from two directories are consecutive ans
                ;;; produce false drop info.  Re-sort them
                tmplist = files[sel]
                tmplist = tmplist(sort(tmplist))
                if( ~keyword_set(sum_in_idl) and rdx_hasopencv() ) then begin
                    flat = rdx_sumfiles(tmplist, time_ave = time_ave, check = check, $
                                      lun = lun, lim = lim, summed = summed, nsum = nsum, verbose=2)
                endif else begin
                    flat = red_sumfiles(tmplist, time_ave = time_ave, check = check, $
                                      lun = lun, lim = lim, summed = summed, nsum = nsum)
                endelse
            endif else begin 
                if( ~keyword_set(sum_in_idl) and rdx_hasopencv() ) then begin
                    flat = rdx_sumfiles(files[sel], time_ave = time_ave, check = check, $
                                      lun = lun, lim = lim, summed = summed, nsum = nsum, verbose=2)
                endif else begin
                    flat = red_sumfiles(files[sel], time_ave = time_ave, check = check, $
                                      lun = lun, lim = lim, summed = summed, nsum = nsum)
                endelse
            endelse

            ;; Remove dark
            flat -= dd
            flat1 = long(temporary(summed))

            ;; Output the raw (if requested) and averaged flats
            ; TODO: output should be written through class-specific methods
            if keyword_set(store_rawsum) then begin
                headerout = 't='+time_ave+' n_sum='+red_stri(nsum)
                print, inam+' : saving ' + sflatname
                fzwrite, flat1, sflatname, headerout
            endif

            headerout = 't='+time_ave+' n_aver='+red_stri(nsum)+' darkcorrected'
            print, inam+' : saving ' + flatname
            fzwrite, float(flat), flatname, headerout

            if keyword_set(check) then free_lun, lun


        endfor  ; states

    endfor  ;  cameras

end  
