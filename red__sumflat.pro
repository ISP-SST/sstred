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
; :Keywords:
; 
;    overwrite  : 
;   
;   
;   
;    ustat  : 
;
;    nthreads  : in, optional, type=integer
;   
;       The number of threads to use for summing.
;   
;    remove  : 
;   
;   
;   
;    check : in, optional, type=boolean
;
;
;    sum_in_rdx : in, optional, type=boolean
;
;      Use rdx_sumfiles.
;   
;    filter  : in, optional, type=int, default=3
;
;       Size of the medianfilter to use when checking data.
;   
;    outdir : in, optional, type=string
;   
;       Write the summed files in this directory instead of in the
;       directory given by get_calib.
;   
;    softlink : in, optional, type=boolean
;   
;       Generate softlinks to the summed files from the file names
;       given by get_calib. 
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
;   2016-05-26 : MGL. Use get_calib method. Base log file name on
;                flatname.
; 
;   2016-05-27 : MGL. Various fixes.
; 
;   2016-05-29 : MGL. Make a FITS header for the summed flat. Save
;                both ANA and FITS format files and do it with
;                red_writedata.
; 
;   2016-05-30 : MGL. Improve the headers.
; 
;   2016-06-09 : MGL. Make the output directory.
; 
;   2016-06-29 : MGL. New keyword outdir.
;
;   2016-08-23 : THI. Rename camtag to detector and channel to camera,
;                so the names match those of the corresponding SolarNet
;                keywords.
; 
;   2016-09-19 : MGL. Changed format strings so exposure times header
;                keywords do not overflow.
;
;   2016-09-21 : THI. Make the size of the medianfilter a parameter.
; 
;   2016-09-22 : MGL. Base output header on relevant original file.
; 
;   2017-03-13 : MGL. Use red_sumheaders.
;
;   2017-07-06 : THI. Get framenumbers and timestamps from rdx_sumfiles
;                and pass them on to red_sumheaders.
;
;   2017-08-07 : MGL. New keyword nthreads.
; 
;   2016-08-10 : MGL. New keyword softlink.
; 
;   2016-08-11 : MGL. When writing the rawsum, write also a FITS
;                version of it.
; 
;   2016-08-23 : MGL. When writing the rawsum, write only the FITS
;                version.
; 
;   2018-04-18 : MGL. Write only in FITS format.
;
;   2025-03-27 : MGL. New keyword dark_timestamp.
;
;-
pro red::sumflat, overwrite = overwrite, $
                  dark_timestamp = dark_timestamp, $
                  ustat = ustat, $
                  nthreads = nthreads, $
                  remove = remove, $
                  cams = cams, $
                  check = check, $
                  lim = lim, $
                  store_rawsum = store_rawsum, $
                  prefilter = prefilter, $
                  dirs = dirs, $
                  outdir = outdir, $
                  softlink = softlink, $
                  sum_in_rdx = sum_in_rdx, $
                  filter = filter

  ;; Defaults
;  if( n_elements(overwrite) eq 0 ) then overwrite = 0
;  if( n_elements(check) eq 0 ) then check = 0
  if( n_elements(dirs) gt 0 ) then dirs = [dirs] $
  else if ptr_valid(self.flat_dir) then dirs = *self.flat_dir
  if( n_elements(cams) gt 0 ) then cams = [cams] $
  else if ptr_valid(self.cameras) then cams = *self.cameras

  ;; Prepare for logging (after setting of defaults).
  ;; Set up a dictionary with all parameters that are in use
  red_make_prpara, prpara, check
  red_make_prpara, prpara, overwrite
  red_make_prpara, prpara, store_rawsum
  red_make_prpara, prpara, sum_in_rdx
  red_make_prpara, prpara, cams
  red_make_prpara, prpara, dirs
  red_make_prpara, prpara, filter
  red_make_prpara, prpara, lim
  red_make_prpara, prpara, outdir
  red_make_prpara, prpara, prefilter
  red_make_prpara, prpara, ustat
 

  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)

  Ncams = n_elements(cams)
  if( Ncams eq 0) then begin
    print, inam+' : ERROR : undefined cams (and cameras)'
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
  
  ;; cameras
  for icam = 0L, Ncams-1 do begin

    cam = cams[icam]
    
    ;; Search files before doing selectfiles so we can limit to the
    ;; wanted prefilter (if used).
    if n_elements(prefilter) ne 0 and strmatch(cam,'Crisp-?') then begin
      files = file_search(dirs + '/' + cam + '/' + '*.'+prefilter+'.*', count = Nfiles)
    endif else undefine, files

    self->selectfiles, cam=cam, dirs=dirs, prefilter=prefilter, ustat=ustat, $
                       files=files, states=states, nremove=remove, /force

    ;; From the above call we expect states to be an array of structures.

    nf = n_elements(files)
    if( nf eq 0 || files[0] eq '') then begin
      print, inam+' : '+cam+': no files found in: '+dirstr
      print, inam+' : '+cam+': skipping camera!'
      continue
    endif else begin
      print, inam+' : Found '+red_stri(nf)+' files in: '+ dirstr + '/' + cam + '/'
    endelse

    state_list = [states[uniq(states.fullstate, sort(states.fullstate))].fullstate]
    state_list = strtrim(state_list,2)
    idx = where( state_list ne '' )
    if min(idx) ge 0 then state_list = state_list[ idx ] $
    else continue

    Nstates = n_elements(state_list)

    ;; Loop over states and sum
    for istate = 0L, Nstates - 1 do begin

      self->selectfiles, prefilter=prefilter, ustat=state_list[istate], $
                         files=files, states=states, selected=sel

      ;; Get the flat file name for the selected state
      self -> get_calib, states[sel[0]] $
                         , flatname = flatname, sflatname = sflatname, status = status

      if n_elements(outdir) ne 0 then begin
        origname = flatname     ; The original flatname, might softlink to this
        sorigname = sflatname   ; The original sflatname, might softlink to this
        flatname = outdir + '/' + file_basename(flatname)
        sflatname = outdir + '/' + file_basename(sflatname)
        sourcename = red_strreplace(flatname,red_strreplace(file_dirname(origname),self.out_dir+'/','')+'/','')
        ssourcename = red_strreplace(sflatname,red_strreplace(file_dirname(sorigname),self.out_dir+'/','')+'/','')
      endif 

      file_mkdir, file_dirname(flatname)

      ;; If file does not exist, do sum!
      if( ~keyword_set(overwrite) && file_test(flatname) ) then begin
        if (~keyword_set(store_rawsum) $
            || file_test(sflatname)) then begin ; only skip if rawsum also exists
          print, inam+' : file exists: ' + flatname $
                 + ' , skipping! (run sumflat, /overwrite to recreate)'
          continue
        endif
      endif
      
      ;; Read the dark frame 
      self -> get_calib, states[sel[0]], darkdata = dd, status = status, timestamp = dark_timestamp
      if status ne 0 then begin
        print, inam+' : no dark found for camera ', cam
        retall
      endif

      if( min(sel) lt 0 ) then begin
        print, inam+' : '+cam+': no files found for state: '+state_list[istate]
        continue
      endif

      print, inam+' : summing flats for state -> ' + state_list[istate]
      print, inam+' : to be saved in ' + flatname
      if(keyword_set(check)) then openw, lun, flatname + '_discarded.txt' $
                                         , width = 500, /get_lun
      
      ;; Sum files
      if keyword_set(check) then begin
        ;; If summing from two directories, same frame numbers from
        ;; two directories are consecutive and produce false drop
        ;; info. Re-sort them.
        tmplist = files[sel]
        tmplist = tmplist(sort(tmplist))

        if( keyword_set(sum_in_rdx) and rdx_hasopencv() ) then begin
          flat = rdx_sumfiles(tmplist, check=check, lun=lun $
                              , nthreads = nthreads $
                              , lim=lim, summed=summed, nsum=nsum, filter=filter $
                              , discarded = discarded, framenumbers = framenumbers $
                              , time_beg=time_beg, time_end=time_end, time_avg=time_avg $
                              , verbose=2)
        endif else begin
          flat = red_sumfiles(tmplist, check = check, lun = lun $
                              , nthreads = nthreads $
;                                 , time_avg = time_avg, time_beg = time_beg, time_end = time_end $
                              , lim = lim, summed = summed, nsum = nsum, filter = filter)
        endelse
      endif else begin 
        if( keyword_set(sum_in_rdx) and rdx_hasopencv() ) then begin
          flat = rdx_sumfiles(files[sel], time_avg = time_avg, check = check $
                              , nthreads = nthreads $
                              , lim = lim, summed = summed, nsum = nsum, filter = filter $
                              , verbose=2)
        endif else begin
          flat = red_sumfiles(files[sel], check = check $
                              , nthreads = nthreads $
;                                  time_avg = time_avg, time_beg = time_beg, time_end = time_end, $
                              , lim = lim, summed = summed, nsum = nsum, filter = filter)
        endelse
      endelse

      ;; Subtract dark and make floating point
      flat = float(flat-dd)

      ;; Make FITS headers 
      head  = red_sumheaders(files[sel], flat, nsum=nsum, framenumbers=framenumbers, $
                             time_beg=time_beg, time_end=time_end, time_avg=time_avg )
      if keyword_set(store_rawsum) then $
         shead = red_sumheaders(files[sel], summed, nsum=nsum, framenumbers=framenumbers, $
                             time_beg=time_beg, time_end=time_end, time_avg=time_avg)
      
      ;; Add some more info here, see SOLARNET deliverable D20.4 or
      ;; later versions of that document. 

      self -> headerinfo_addstep, head, prstep = 'SUMMING' $
                                  , prproc = inam, prpara = prpara
      if keyword_set(store_rawsum) then $
         self -> headerinfo_addstep, shead, prstep = 'SUMMING' $
                                     , prproc = inam, prpara = prpara
      
      ;; Write FITS format averaged flat
      print, inam+' : saving ', flatname
      fxaddpar, head, 'FILENAME', file_basename(flatname), after = 'DATE'
      red_writedata, flatname, flat, header=head, filetype='FITS', overwrite = overwrite

      
      ;; Output the raw (if requested) flats
      if keyword_set(store_rawsum) then begin
        ;; FITS
        print, inam+' : saving ', sflatname
        fxaddpar, shead, 'FILENAME', file_basename(sflatname), after = 'DATE'
        red_writedata, sflatname, long(temporary(summed)), header=shead, filetype='FITS', overwrite = overwrite
      endif

      if keyword_set(check) then free_lun, lun

      if keyword_set(softlink) and n_elements(outdir) ne 0 then begin
        file_delete, origname, /allow_nonexistent
        file_link, sourcename, origname
        if keyword_set(store_rawsum) then begin
          file_delete, sorigname, /allow_nonexistent
          file_link, ssourcename, sorigname
        endif
      endif

    endfor                      ; istate

  endfor                        ; icam

end  
