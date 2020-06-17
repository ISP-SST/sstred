; docformat = 'rst'

;+
; Sum pinhole images.
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
; :Keywords:
; 
;    nthreads : in, optional, type=integer, default=2
;   
;      The number of threads to use for summing and bad-pixel filling.
;   
;    no_descatter : in, optional, type=boolean 
;   
;      Don't do back-scatter compensation.
;   
;    ustat : in, optional, type=string
;   
;      Do only for this state.
;   
;    prefilter : in, optional, type=string
;   
;       Set this keyword to the prefilter you want to sum pinholes
;       from. Otherwise, sum pinholes for all prefilters.
;   
;    pinhole_align : in, optional, type=boolean
; 
;       If true, then perform subpixel alignment of pinholes before
;       summing. 
; 
;    brightest_only : in, optional, type=boolean
;
;       Set this to only sum (one of) the brightest tunings for each
;       prefilter. 
;
;    lc_ignore : in, optional, type=boolean
;
;       Set this to treat all lc states as if they were the same.
; 
; :History:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2013-08-15 : MGL. Added keyword pinhole_align. Process cameras in
;                loops instead of separately.
; 
;   2013-08-27 : MGL. Added support for logging. Let the subprogram
;                find out its own name.
; 
;   2013-08-29 : MGL. Made use of nthread consistent. Loop over
;                cameras rather than treating them one by one.
; 
;   2013-08-30 : MGL. When keyword pinhole_align is set, request that
;                red_sumfiles do: dark, flat, fillpix, alignment
;                before summing.
; 
;   2013-09-02 : MGL. Two new keywords, brightest_only and lc_ignore
;                (but they don't actually do anything yet).
; 
;   2013-09-03 : MGL. Fixed descattering bug.
; 
;   2013-09-04 : MGL. Store also in floating point. This is the
;                version used in red::pinholecalib.pro.
; 
;   2016-02-15 : MGL. Use loadbackscatter. Remove keyword descatter,
;                new keyword no_descatter.
; 
;   2016-05-31 : THI. Re-write to use the class-methods.
; 
;   2016-08-23 : THI. Rename camtag to detector and channel to camera,
;                so the names match those of the corresponding SolarNet
;                keywords.
; 
;   2017-04-13 : MGL. Construct SOLARNET metadata fits header.
; 
;   2017-07-06 : THI. Get framenumbers and timestamps from rdx_sumfiles
;                and pass them on to red_sumheaders.
; 
;   2016-08-10 : MGL. New keyword outdir.
;
;
;-
pro red::sumpinh, nthreads = nthreads $
                  , no_descatter = no_descatter $
                  , ustat = ustat $
                  , prefilter = prefilter $
                  , cams = cams $
                  , dirs = dirs $
                  , outdir = outdir $
                  , pinhole_align = pinhole_align $
                  , brightest_only = brightest_only $
                  , lc_ignore = lc_ignore $
                  , overwrite = overwrite $
                  , sum_in_rdx = sum_in_rdx

  if( n_elements(dirs) gt 0 ) then dirs = [dirs] $
  else if ptr_valid(self.pinh_dirs) then dirs = *self.pinh_dirs
  
  if(n_elements(cams) eq 0 and ptr_valid(self.cameras)) then cams = *self.cameras
  if(n_elements(cams) eq 1) then cams = [cams]

  ;; Prepare for logging (after setting of defaults).
  ;; Set up a dictionary with all parameters that are in use
  red_make_prpara, prpara, no_descatter
  red_make_prpara, prpara, overwrite
  red_make_prpara, prpara, sum_in_rdx
  red_make_prpara, prpara, pinhole_align
  red_make_prpara, prpara, brightest_only 
  red_make_prpara, prpara, cams
  red_make_prpara, prpara, dirs
  red_make_prpara, prpara, prefilter
  red_make_prpara, prpara, ustat
  red_make_prpara, prpara, outdir


  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])
  
;  ;; Logging
;  help, /obj, self, output = selfinfo 
;  red_writelog, selfinfo = selfinfo

  if n_elements(nthreads) eq 0 then nthread = 2 else nthread = nthreads

  Ncams = n_elements(cams)
  if( Ncams eq 0) then begin
    print, inam+' : ERROR : undefined cams (and cameras)'
    return
  endif

  Ndirs = n_elements(dirs)
  if( Ndirs eq 0) then begin
    print, inam+' : ERROR : no pinhole directories defined'
    return
  endif else begin
    if Ndirs gt 1 then begin
      dirstr = '['+ strjoin(dirs,';') + ']'
      print,'WARNING: sumpinh was called with multiple directories.'
      print,'Only the first directory will be used.'
      print,"To use a particular directory, use: a->sumpinh,dir='path/to/phdata'"
      dirs = [dirs[0]]
    endif else dirstr = dirs[0]
    dirstr = dirs[0]            ; remove when we can properly deal with multiple directories.
  endelse

  if ~file_test(dirs,/directory) then begin
    print, inam + ' : ERROR : "'+dirs+'" is not a directory'
    return
  endif
  
  ;; Loop over cameras
  for icam = 0, Ncams-1 do begin

    cam = cams[icam]
    detector = self->getdetector( cam )

    if 0 then begin
      ;; Needed for old data? E.g., 2012?
      if n_elements(prefilter) gt 0 then prefin = prefilter else undefine, prefin
      files = red_raw_search(dirs[0], count = nf) ;, prefilters = prefin)
      self -> extractstates, files, states
      self->selectfiles, cam=cam, dirs=dirs, prefilter=prefin, ustat=ustat, $
                         files=files, states=states ;, /force
    endif else begin
      self->selectfiles, cam=cam, dirs=dirs, prefilter=prefin, ustat=ustat, $
                         files=files, states=states, /force
    endelse
    
    nf = n_elements(files)
    if( nf eq 0 || files[0] eq '') then begin
      print, inam+' : '+cam+': no files found in: '+dirstr
      print, inam+' : '+cam+': skipping camera!'
      continue
    endif else begin
      print, inam+' : Found '+red_stri(nf)+' files in: '+ dirstr + '/' + cam + '/'
    endelse

    ;; We sum separately based on fpi_state, this should match the
    ;; definition of pinhname as generated by get_calib.
    state_list = states[uniq(states.fpi_state, sort(states.fpi_state))]

    Nstates = n_elements(state_list)
    ;; Loop over states and sum
    for istate = 0L, Nstates - 1 do begin
      
      this_state = state_list[istate]
      sel = where(states.fpi_state eq this_state.fpi_state)
      
      ;; Get the pinholes file name for the selected state, as well as
      ;; dark and flat data.
      self -> get_calib, states[sel[0]], darkdata = dd, flatdata = ff, $
                         darkname = darkname,  flatname = flatname, $
                         pinhname = pinhname, status = status
      if( status ne 0 ) then begin
        print, inam+' : failed to load calibration data for:', states[sel[0]].filename
        print, 'One of these seems to be missing:'
        print, 'darkname=', darkname
        print, 'flatname=', flatname
        print, 'pinhname=', pinhname
        continue
      endif

      if n_elements(outdir) ne 0 then pinhname = outdir + '/' + file_basename(pinhname)

      ;; If file does not exist, do sum!
      if( ~keyword_set(overwrite) && file_test(pinhname) ) then begin
        print, inam+' : file exists: ' + pinhname + ' , skipping! (run sumpinh, /overwrite to recreate)'
        continue
      endif

      if( n_elements(sel) lt 1 || min(sel) lt 0 ) then begin
        print, inam+' : '+cam+': no files found for state: '+state_list[istate].fullstate
        continue
      endif
      
      pref = states[sel[0]].prefilter
      print, inam+' : adding pinholes for state -> ' + state_list[istate].fpi_state

      DoBackscatter = 0
      if (~keyword_set(no_descatter) AND self.dodescatter AND (pref eq '8542' OR pref eq '7772')) then begin
        self -> loadbackscatter, detector, pref, bgt, Psft
        DoBackscatter = 1
      endif
      if DoBackscatter gt 0 then begin
        ff = rdx_descatter(ff, bgt, Psft, /verbose, nthreads = nthread)
      endif
      gain = self->flat2gain(ff)

      ;; Sum files

      ;; Dark and flat correction and bad-pixel filling done by
      ;; red_sumfiles on each frame before alignment.
      if rdx_hasopencv() and keyword_set(sum_in_rdx) then begin
        psum = rdx_sumfiles(files[sel], pinhole_align=pinhole_align, dark=dd, $
                            nthreads = nthreads, $
                            gain=gain, backscatter_gain=bgt, backscatter_psf=Psft, $
                            nsum=nsum, framenumbers=framenumbers, time_beg=time_beg, $
                            time_end=time_end, time_avg=time_avg, verbose=2 )
      endif else begin
        psum = red_sumfiles(files[sel], pinhole_align=pinhole_align, dark=dd, gain=gain, $
                            nthreads = nthreads, $
                            backscatter_gain=bgt, backscatter_psf=Psft, nsum=nsum)
      endelse

      ;; Make FITS headers 
      head  = red_sumheaders(files[sel], psum, nsum=nsum, framenumbers=framenumbers, $
                             time_beg=time_beg, time_end=time_end, time_avg=time_avg)
      fxaddpar, head, 'FILENAME', file_basename(pinhname), after = 'DATE'

      ;; Add some more info here, see SOLARNET deliverable D20.4 or
      ;; later versions of that document.

      self -> headerinfo_addstep, head, prstep = 'SUMMING' $
                                  , prproc = inam, prpara = prpara

      file_mkdir, file_dirname(pinhname)

      print, inam+' : saving ' + pinhname
      if keyword_set(pinhole_align) then begin
        red_writedata, pinhname, psum, header=head, overwrite = overwrite
      endif else begin
        red_writedata, pinhname, fix(round(10. * psum)), header=head, overwrite = overwrite
      endelse

    endfor                      ; istate

  endfor                        ; icam

end
