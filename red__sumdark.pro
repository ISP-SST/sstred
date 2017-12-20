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
;    check  : 
;   
;   
;    cams : in, optional, type=strarr
; 
;      A list of cameras (or rather camera subdirs).
;
;    nthreads  : in, optional, type=integer
;   
;       The number of threads to use for summing.
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
; :History:
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
;   2016-05-17 : THI. Use cameras by default, added keyword dirs to
;                use specific dark-folder. Re-write to sum files from
;                multiple folders at once.
; 
;   2016-05-18 : MGL. New keyword sum_in_idl. Started working on
;                making a header for the dark frame.
; 
;   2016-05-23 : MGL. Make the fz outheader here instead of calling
;                red_dark_h2h. Remove some unused code.
;
;   2016-05-25 : MGL. Don't store the non-normalized darks, momfbd
;                doesn't need them. Do store darks in FITS format.
;                Don't use red_sxaddpar. Write darks out as float, not
;                double.
;
;   2016-05-26 : MGL. Use get_calib method to get the file name. Don't
;                call get_calib one extra time just to get the
;                directory.
;
;   2016-06-01 : THI. Loop over states (gain & exposure)
;
;   2016-06-09 : MGL. Tell red_writedata when you want to store as
;                FITS.
;
;   2016-08-23 : THI. Rename camtag to detector and channel to camera,
;                so the names match those of the corresponding
;                SOLARNET keywords.
;
;   2016-09-21 : MGL. Tell user what darkname is being summed for
;                instead of what state.
;
;   2016-09-21 : THI. Make the size of the medianfilter a parameter.
;
;   2016-09-21 : MGL. Put DATE-BEG, DATE-END, DATE-AVE in output
;                header.
;
;   2017-03-13 : MGL. Use red_sumheaders.
;
;   2017-03-16 : MGL. Add FITS headers about this processing step.
;
;   2017-07-06 : THI. Get framenumbers and timestamps from rdx_sumfiles
;                and pass them on to red_sumheaders.
;
;   2017-08-07 : MGL. New keyword nthreads.
; 
;   2016-08-10 : MGL. New keywords outdir and softlink.
;
;-
pro red::sumdark, overwrite = overwrite, $
                  check = check, $
                  cams = cams, $
                  dirs = dirs, $
                  outdir = outdir, $
                  softlink = softlink, $
                  nthreads = nthreads, $
                  sum_in_rdx = sum_in_rdx, $
                  filter = filter

  ;; Defaults
;  if n_elements(overwrite) eq 0 then overwrite = 0
;  if n_elements(check) eq 0 then check = 0
  if n_elements(dirs) gt 0 then dirs = [dirs] $
  else if ptr_valid(self.dark_dir) then dirs = *self.dark_dir
  if n_elements(cams) gt 0 then cams = [cams] $
  else if ptr_valid(self.cameras) then cams = *self.cameras

  ;; Prepare for logging (after setting of defaults).
  ;; Set up a dictionary with all parameters that are in use
  red_make_prpara, prpara, overwrite
  red_make_prpara, prpara, check
  red_make_prpara, prpara, sum_in_rdx
  red_make_prpara, prpara, cams
  red_make_prpara, prpara, dirs
  red_make_prpara, prpara, filter
  red_make_prpara, prpara, outdir
 
  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

;  ;; Logging
;  help, /obj, self, output = selfinfo1
;  help, /struct, self.done, output = selfinfo2 
;  red_writelog, selfinfo = [selfinfo1, selfinfo2]

  Ncams = n_elements(cams)
  if( Ncams eq 0) then begin
    print, inam+' : ERROR : undefined cams (and cameras)'
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

  for icam = 0L, Ncams-1 do begin

    cam = cams[icam]

    self->selectfiles, cam=cam, dirs=dirs, ustat=ustat, $
                       files=files, states=states, /force

    nf = n_elements(files)
    if( nf eq 0 || files[0] eq '') then begin
      print, inam+' : '+cam+': no files found in: '+dirstr
      print, inam+' : '+cam+': skipping camera!'
      continue
    endif else begin
      print, inam+' : Found '+red_stri(nf)+' files in: '+ dirstr + '/' + cam + '/'
    endelse

    state_list = [states[uniq(states.fullstate, sort(states.fullstate))].fullstate]

    Nstates = n_elements(state_list)

    ;; Loop over states and sum
    for istate = 0L, Nstates - 1 do begin

      self->selectfiles, ustat=state_list[istate], $
                         files=files, states=states, selected=sel

      ;; Get the name of the darkfile
      self -> get_calib, states[sel[0]], darkname = darkname, status = status
      if status ne 0 then stop
      
      if n_elements(outdir) ne 0 then begin
        origname = darkname     ; The original darkname, might softlink to this
        darkname = outdir + '/' + file_basename(darkname)
        sourcename = red_strreplace(darkname,red_strreplace(file_dirname(origname),self.out_dir+'/','')+'/','')
      endif
      
      file_mkdir, file_dirname(darkname)

      if( ~keyword_set(overwrite) && file_test(darkname) && file_test(darkname+'.fits')) then begin
        print, inam+' : file exists: ' + darkname + ' , skipping! (run sumdark, /overwrite to recreate)'
        continue
      endif

      if( min(sel) lt 0 ) then begin
        print, inam+' : '+cam+': no files found for state: '+state_list[istate]
        continue
      endif

      print, inam+' : summing darks -> ' + file_basename(darkname)
      if(keyword_set(check)) then begin
        openw, lun, darkname + '_discarded.txt', width = 500, /get_lun
      endif
      
      ;; Do the summing
      if rdx_hasopencv() and keyword_set(sum_in_rdx) then begin
        dark = rdx_sumfiles(files[sel], check=check, lun=lun, summed=darksum, nsum=nsum, filter=filter $
                            , nthreads = nthreads $
                            , framenumbers=framenumbers $
                            , time_beg=time_beg, time_end=time_end, time_avg=time_avg, verbose=2)
      endif else begin
        dark = red_sumfiles(files[sel], check = check, lun = lun, summed = darksum $
                            , nthreads = nthreads $
                            , nsum=nsum, filter=filter)
        ;;$
        ;;                    , time_avg = time_avg, time_beg = time_beg, time_end = time_end)
      endelse

      ;; The momfbd code can't read doubles.
      dark = float(dark)      
      
      ;; Make FITS headers 
      head = red_sumheaders(files[sel], dark, nsum=nsum, framenumbers=framenumbers, $
                            time_beg=time_beg, time_end=time_end, time_avg=time_avg )
;      headsum  = red_sumheaders(files[sel], darksum, nsum = nsum)
      
      ;; Add some more info here, see SOLARNET deliverable D20.4 or
      ;; later versions of that document.

      self -> headerinfo_addstep, head, prstep = 'Dark summing' $
                                  , prproc = inam, prpara = prpara

      ;; Write ANA format dark
      print, inam+' : saving ', darkname
      fxaddpar, head, 'FILENAME', file_basename(darkname), after = 'DATE'
      red_writedata, darkname, dark, header=head, filetype='ana', overwrite = overwrite

      ;; Write FITS format dark
      print, inam+' : saving ', darkname+'.fits'
      fxaddpar, head, 'FILENAME', file_basename(darkname)+'.fits'
      red_writedata, darkname+'.fits', dark, header=head, filetype='fits', overwrite = overwrite
      
      if keyword_set(check) then begin
        free_lun, lun
      endif

      if keyword_set(softlink) and n_elements(outdir) ne 0 then begin
        file_delete, origname, /allow_nonexistent
        file_link, sourcename, origname
        file_delete, origname+'.fits', /allow_nonexistent
        file_link, sourcename+'.fits', origname+'.fits'
      endif

    endfor                      ; states

  endfor                        ; icam 

  return

end
