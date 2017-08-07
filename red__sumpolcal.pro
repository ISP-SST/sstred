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
;    remove  : 
;   
;   
;   
;    ucam  : 
;   
;   
;   
;    check : in, optional, type=boolean
;
;    nthreads  : in, optional, type=integer
;   
;       The number of threads to use for summing and backscatter correction.
;
;
;    sum_in_rdx : in, optional, type=boolean
;
;      Use rdx_sumfiles.
;
;    outdir : in, optional, type=string
;   
;   
; 
; 
; :History:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2013-07-09 : MGL. Worked around outdir1, outdir2, and lun when not
;                specifying /old.
; 
; 
;   2013-07-10 : MGL. Worked around fzhead bug by using fzread
;                instead. 
;
;   2013-08-27 : MGL. Added support for logging. Let the subprogram
;                find out its own name.
;
;   2017-07-06 : MGL. Adapt to new pipeline mechanisms. Add keywords
;                dirs, outdir, sum_in_rdx. Add metadata handling.
;                Remove keyword old.
;
;   2017-08-07 : MGL. New keyword nthreads.
; 
;
;-
pro red::sumpolcal, check=check $
                    , dirs = dirs $
                    , nthreads = nthreads $
                    , outdir = outdir $
                    , overwrite = overwrite $
                    , remove = remove $
                    , sum_in_rdx = sum_in_rdx $
                    , ucam = ucam 
           
  ;; Defaults
  if( n_elements(dirs) gt 0 ) then dirs = [dirs] $
  else if ptr_valid(self.polcal_dir) then dirs = *self.polcal_dir

  ;; Prepare for logging (after setting of defaults).
  ;; Set up a dictionary with all parameters that are in use
  prpara = dictionary()
  ;; Boolean keywords
  if keyword_set(check) then prpara['check'] = check
  if keyword_set(sum_in_rdx) then prpara['sum_in_rdx'] = sum_in_rdx
  ;; Non-boolean keywords
  if n_elements(dirs) ne 0 then prpara['dirs'] = dirs
  if n_elements(outdir) ne 0 then prpara['outdir'] = outdir
  if n_elements(ucam) ne 0 then prpara['ucam'] = ucam

  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  self -> getdetectors 
  cams = *self.cameras
  detectors = *self.detectors
  Ncams = n_elements(cams)
  if Ncams eq 0 then begin
    print, inam+' : ERROR : undefined cams (and cameras)'
    return
  endif

  ;; Loop cameras
  for icam = 0, n_elements(cams)-1 do begin

    cam = cams[icam]
    
    if strmatch(cam,'*-[DW]') then continue ; Wideband cameras don't have polcal data
    
    if keyword_set(ucam) then begin
      if cam ne ucam then begin
        print, inam + ' : skipping '+cam
        continue
      endif
    endif

    
    self->selectfiles, cam=cam, dirs=dirs, prefilter=prefilter, ustat=ustat, $
                       files=files, states=states, nremove=remove, /force, /polcal

    
    state_list = [states[uniq(states.fullstate, sort(states.fullstate))].fullstate]

    Nstates = n_elements(state_list)

;    ;; Find files
;    spawn, 'find ' + self.polcal_dir + '/' + cam + "/ | grep im.ex | grep -v '.lcd.'", files
;    
;    nt = n_elements(files)
;    if(files[0] eq '') then begin
;      print, inam+' : no files found in '+self.polcal_dir+'/'+cam+', skipping camera!'
;      continue
;    endif
;
;    ;; Sort files based on image number
;    files = red_sortfiles(files)


    
;    ;; Get image states
;    pstat = red_getstates_polcal(files)

    ; if(keyword_set(remove)) then pstat.star[*] = red_flagchange(pstat.qw)

    ;; Unique states
;    ustate = pstat.state[uniq(pstat.state, sort(pstat.state))]
;    Nstates = n_elements(ustate)

    print, inam+' : found '+red_stri(Nstates)+' states for '+cam

    camtag = (strsplit(file_basename(files[0]), '.',/extract))[0]

    for istate = 0L, Nstates - 1 do begin 
 
      self->selectfiles, prefilter=prefilter, ustat=state_list[istate], $
                         files=files, states=states, selected=sel, /polcal

;      pos = where((pstat.state eq ustate[istate]) AND (pstat.star eq 0B), count)
;      if(count eq 0) then continue


      ;; Get the polcal file name for the selected state
      self -> get_calib, states[sel[0]] $
                         , polsname = polsname, status = status
;      polsname = polsname[0]

      if n_elements(outdir) ne 0 then begin
        polsname = outdir + '/' + cam + '/' + file_basename(polsname)
      endif

      file_mkdir, file_dirname(polsname)

      print, inam+' : summing polcal for state -> ' + state_list[istate]
      print, inam+' : to be saved in ' + polsname
      if(keyword_set(check)) then openw, lun, red_strreplace(polsname,'.fits','_discarded.txt') $
                                         , width = 500, /get_lun
      
      filelist = files[sel]
      filelist = filelist(sort(filelist))

      print, inam+' : summing frames for '+cam+' -> '+state_list[istate]
      if keyword_set(sum_in_rdx) and rdx_hasopencv() then begin
        pcal = rdx_sumfiles(filelist, lun = lun, lim = lim $
                            , nthreads = nthreads $
                            , nsum = nsum, filter = filter $
                            , check = check, discarded = discarded, framenumbers = framenumbers $
                            , time_beg = time_beg, time_end = time_end, time_avg = time_avg $
                            , verbose=2)
      endif else begin
        pcal = red_sumfiles(filelist, check = check, lun = lun, lim = lim $
                            , nthreads = nthreads $
;                                  time_avg = time_avg, time_beg = time_beg, time_end = time_end, $
                            , nsum = nsum, filter = filter)
      endelse

      ;; Make FITS headers 
      head  = red_sumheaders(filelist, pcal, nsum=nsum, framenumbers=framenumbers, $
                             time_beg=time_beg, time_end=time_end, time_avg=time_avg )

      ;; Add some more info here, see SOLARNET deliverable D20.4 or
      ;; later versions of that document. 
      red_fitsaddpar, head, 'STATE', state_list[istate], 'Polcal state'
      self -> headerinfo_addstep, head, prstep = 'Polcal summing' $
                                  , prproc = inam, prpara = prpara

;      ;; Write ANA format pcal
;      print, inam+' : saving ', polcname
;      fxaddpar, head, 'FILENAME', file_basename(polcname), after = 'DATE'
;      red_writedata, polcname, polc, header=head, filetype='ANA', overwrite = overwrite

      ;; Write FITS format pcal
      print, inam+' : saving ', polsname
      red_fitsaddpar, head, 'FILENAME', file_basename(polsname), after = 'DATE'
      red_writedata, polsname, pcal, header=head, filetype='FITS', overwrite = overwrite

      if keyword_set(check) then free_lun, lun

    endfor                      ; istate

  endfor                        ; icam

end
