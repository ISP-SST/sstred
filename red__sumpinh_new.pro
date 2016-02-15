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
;       The number of threads to use for backscattering and bad-pixel filling.
;   
;    no_descatter : in, optional, type=boolean 
;   
;      Don't do back-scatter compensation.
;   
;    ustat : in, optional, type=strarr
;   
;      Process these states. Overrides keywords all_tunings and
;      all_lcstates. Do not use together with keyword pref.
;   
;    pref : in, optional, type=string
;   
;       Set this keyword to the prefilter you want to sum pinholes
;       from. Otherwise, sum pinholes for all prefilters. Do not use
;       together with keyword ustat.
;   
;    pinhole_align : in, optional, type=boolean
; 
;       If true, then perform subpixel alignment of pinholes before
;       summing. 
; 
;    all_tunings : in, optional, type=boolean
;
;       Set this to process all lambda tunings. Otherwise process only
;       the brightest tuning for each prefilter.
;
;    all_lcstates : in, optional, type=boolean
;
;       Set this to process all lc states one by one, otherwise just
;       pick one.
;
;    overwrite : in, optional, type=boolean
;
;       Overwrite existing output files.
; 
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
;   2014-01-13 : MGL. New boolean for whether or not to do descatter.
;                Added documentation and changed some variable names.
;                Now stop if descatter is wanted but the backscatter
;                psf or bgain does not exist.
; 
;   2014-01-14 : MGL. Removed the (never used) keywords brightest_only
;                and lc_ignore. Instead new keywords expressing the
;                opposite choice (and therefore opposite default):
;                all_tunings and all_lcstates. Removed alternate
;                outdir. Implemented selection of brightest tuning and
;                picking a single LC state.
; 
;   2014-01-15 : MGL. Make an outer prefilter loop. The prefilter must
;                be known at the point where we select states.
;                Simplify use of nthreads keyword.
; 
;   2014-01-16 : MGL. New keyword: overwrite. Save gain in the pinh/
;                subdirectory to be used for display during the
;                getalignclips step.
; 
;   2016-02-15 : MGL. Use loadbackscatter. Remove keyword descatter,
;                new keyword no_descatter.
; 
; 
; 
; 
;-
pro red::sumpinh_new, nthreads = nthreads $
                      , no_descatter = no_descatter $
                      , ustat = ustat $
                      , pref = epref $
                      , pinhole_align = pinhole_align $
                      , all_tunings = all_tunings $
                      , all_lcstates = all_lcstates $
                      , overwrite = overwrite

  ;; TODO:
  ;;
  ;; * Test with /descatter, particularly together with
  ;;   /pinhole_align.
  ;;
  ;; * For wide lines like Halpha, implement mechanism for not using a
  ;;   single tuning by default. Perhaps have some member of the
  ;;   object be a list of "wide" lines that require either
  ;;   /all_tunings or a selection of tunings and interpolation.


  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  ;; Logging
  help, /obj, self, output = selfinfo 
  red_writelog, selfinfo = selfinfo
  
  if ~self.dopinh then begin
     print, inam + ' : ERROR : undefined pinh_dir'
     return
  endif
  if n_elements(nthreads) eq 0 then nthreads = 2

  ;; Output dir
  outdir = self.out_dir+ '/' + 'pinh/'
  file_mkdir, outdir

  ;; Cameras
  cams = [self.camt, self.camr, self.camwb]
  Ncams = n_elements(cams)

  ;; Selected both prefilter and state?
  if n_elements(epref) ne 0 and n_elements(ustat) ne 0 then begin
      print, inam + ' : ERROR : Use only one of keywords pref and ustat.'
      retall
  endif

  ;; Loop over cameras
  for icam = 0, Ncams-1 do begin

     ;; Create file list for this camera
     spawn, 'find ' + self.pinh_dir + '/' + cams[icam] + "/ | grep im.ex |  grep -v '.lcd.'", files
     Nfiles = n_elements(files)
     ;;nr = n_elements(rfiles)
     ;;nw = n_elements(wfiles)  

     if (files[0] eq '') then begin
        print, inam + ' : files not found for this camera:'
        print, '   '+cams[icam]+' -> '+ Nfiles
        return
     endif
     
     ;; Sort files
     files = red_sortfiles(files)

     ;; Get states for the cameras
     print, inam + ' : extracting states for camera '+cams[icam] 
     stat = red_getstates(files)

     ;; Flagging first image after tuning (will use transmitted camera for the states)
     state = stat.state
;     rstate = rstat.state
;     wstate = wstat.state

     ;; Deal with prefilters
     
     prefilters = stat.pref[uniq(stat.pref, sort(stat.pref))]
     if n_elements(epref) gt 0 then begin
        ;; Do any files match the pref keyword (variable epref)?
        indx = where(prefilters eq epref)
        if max(indx) eq -1 then begin
           ;; No matching prefilters
           print, inam + ' ERROR : No files match the pref keyword: ', epref
           retall
        endif else begin
           ;; Process only the wanted prefilter
           prefilters = [epref]
        endelse
     endif
     if n_elements(ustat) gt 0 then begin
        ;; Do any files match the prefilters in ustat?
        mask = bytarr(n_elements(ustat))
        for ii = 0, n_elements(ustat)-1 do begin
           upref = (strsplit(ustat[ii], '.', /extract))[0]
           mask[ii] = max(where(prefilters eq upref)) ne -1 ; 1 if match
        endfor
        if min(mask) eq 0 then begin
           ;; Mismatch
           print, inam + ' ERROR : No files match at least one the prefilters in the ustat keyword: ' $
                  , ustat[where(mask)]
           retall
        endif
     endif
     Nprefilters = n_elements(prefilters)
     
     ;; At this point we should know that we will only be looping over
     ;; prefilters that correspond to at least some files.
     
     for ipref = 0, Nprefilters-1 do begin

        pref = prefilters[ipref]
  
        ;; Decide which states to process. Also make name tags for the
        ;; file names.
        if n_elements(ustat) ne 0 then begin
           ;; If ustat is given as a keyword, then use it.
           selected_states = ustat
           name_tags = '.' + ustat
        endif else begin
           ;; Base selction on status of keywords all_tunings and all_lcstates.
           possible_indx = uniq(state, sort(state))
           possible_states = state[possible_indx]
           ;; Filter valid states according to current prefilter (pref)
           filtered_states = where(strmatch(possible_states, '*'+pref+'*'))
           if n_elements(filtered_states) eq 0 or filtered_states[0] eq -1 then continue
           possible_states = possible_states[filtered_states]
           possible_indx = possible_indx[filtered_states]
           if keyword_set(all_tunings) and keyword_set(all_lcstates) then begin
              ;; Process all possible states
              selected_states = possible_states
              name_tags = '.' + selected_states ; All state info in the file name tag
           endif else if keyword_set(all_tunings) then begin
              ;; Pick one LC state, process all tunings
              picked_lcstate = (strsplit(possible_states[0],'.',/extract))[2]
              lcindx = possible_indx(where(strmatch(possible_states, '*'+picked_lcstate)))
              selected_states = state[lcindx]
              name_tags = '.' + stat.pref[lcindx] + '.' + stat.wav[lcindx] ; Only prefilter and tuning
           endif else begin
              ;; Pick brightest tuning based on distance from line core.
              ;; Note that there can be multiple line cores (6301 & 6302
              ;; lines).
              tmp = max(abs(stat.dwav[possible_indx]-float(stat.wav[possible_indx])), ml)           
              if keyword_set(all_lcstates) then begin
                 ;; Process all LC states with this tuning
                 picked_tuning = strjoin((strsplit(state[possible_indx[ml]],'.',/extract))[0:1],'.')
                 tunindx = possible_indx(where(strmatch(possible_states, picked_tuning+'*')))
                 selected_states = state[tunindx]
                 name_tags = '.' + stat.pref[tunindx] + '.' + stat.lc[tunindx] ; Only prefilter and LC states
              endif else begin
                 ;; Just go with the LC state we happened to pick
                 selected_states = state[possible_indx[ml]] 
                 name_tags = '.' + stat.pref[possible_indx[ml]] ; Only the prefilter
              endelse
           endelse
        endelse
        
        ;; Get camera tag
        cam = red_camtag(files[0])
        
        ;; Load darks
        dark = f0(self.out_dir + '/darks/'+cam+'.dark')
        
        ;; Loop over selected_states
        firsttime = 1B
        Nstates = n_elements(selected_states)
        for ii = 0L, Nstates-1 do begin
           
;        pref = (strsplit(selected_states[ii], '.',/extract))[0]
;
;        if n_elements(epref) gt 0 then begin
;           if pref ne epref then continue
;        endif
           
           ;; File names for saving
           pname = cam + name_tags[ii] + '.pinh'
           pfname = cam + name_tags[ii] + '.fpinh'
           
           AlreadyDone =  file_test(outdir+pname) and (file_test(outdir+pfname) or ~keyword_set(pinhole_align))
           if AlreadyDone and ~keyword_set(overwrite) then begin
              
              print, inam + ' : Already done:'
              print, pname
              if keyword_set(pinhole_align) then begin
                 print, pfname
              endif
              
              continue
              
           endif                ; AlreadyDone
           
           if cams[icam] eq self.camwb then begin
              flatf = self.out_dir + '/flats/' + strjoin([cam, pref, 'flat'],'.')
           endif else begin
              flatf = self.out_dir + '/flats/' + strjoin([cam, selected_states[ii], 'flat'],'.')
           endelse
           
           if(file_test(flatf)) then begin
              flat = f0(flatf)
           endif else begin
              print, inam + ' : ERROR -> flat not found for '+cam+'.'+selected_states[ii]
              stop
           endelse
           
           ;; Boolean for whether we want to do descatter for this
           ;; prefilter
           DoDescatter = ~keyword_set(no_descatter) AND self.dodescatter AND (pref eq '8542' OR pref eq '7772')
           
           ;; Descatter data?
           if DoDescatter then begin
              
              if(firsttime) then begin
                 
                 self -> loadbackscatter, cam, pref, descatter_bgain, descatter_psf
                                  
;                 descatter_psf_name =  self.descatter_dir+ '/' + cam + '.psf.f0'
;                 descatter_bgain_name = self.descatter_dir+ '/' + cam + '.backgain.f0'
;                 
;                 if file_test(descatter_psf_name) and file_test(descatter_bgain_name) then begin
;                    descatter_psf   = f0(descatter_psf_name)
;                    descatter_bgain = f0(descatter_bgain_name)
;                 endif else begin
;                    print, inam + ' : ERROR -> backscatter psf and/or gain not found:'
;                    print, descatter_psf_name
;                    print, descatter_bgain_name
;                    stop    
;                 endelse
                 
                 firsttime = 0B
                 
              endif             ; firsttime
              
              flat = red_cdescatter(flat, descatter_bgain, descatter_psf, /verbose, nthreads = nthreads)
              
           endif                ; DoDescatter
           
           gain = red_flat2gain(flat)
           
           print, inam+' : summing pinh for state -> ' + selected_states[ii]
           pos = where((state eq selected_states[ii]), count)
           if count eq 0 then begin
              
              print, inam + ' : WARNING-> No files found for camera '+cams[icam]+' -> '+ selected_states[ii]
              
           endif else begin
              
              if keyword_set(pinhole_align) then begin
                 
                 ;; Dark and flat correction and bad-pixel filling done
                 ;; by red_sumfiles on each frame before alignment.
                 c = red_sumfiles(stat.files[pos], /pinhole_align, dark = dark, gain = gain)
                 
                 ;; Descatter not taken care of yet
                 
              endif else begin
                 
                 ;; Dark and flat correction, possibly descattering, and
                 ;; then bad-pixel filling done here after summing.
                 c = red_sumfiles(stat.files[pos]) - dark
                 
                 if DoDescatter then begin
                    c = red_cdescatter(c, descatter_bgain, descatter_psf, /verbose, nthreads = nthreads)
                 endif
                 
                 print, inam + ' : Filling pixels'
                 c = red_fillpix(temporary(c) * gain, mask=gain ne 0, nthreads = nthreads)
                 
              endelse
              
              ;; Save
              head = 'n_aver=' + red_stri(count)
              
              ;;namout = cam+'.' +selected_states[ii]+'.pinh'
              print, inam+' : saving ' + outdir + pname
              fzwrite, fix(round(10. * c)), outdir+pname, head

              if keyword_set(pinhole_align) then begin
                 ;;namout = cam+'.' +selected_states[ii]+'.fpinh'
                 print, inam+' : saving ' + outdir + pfname
                 fzwrite, c, outdir+pfname, head
              endif
              
              ;; Save gain to be used for displaying during the
              ;; getalignclips step.
              fzwrite, gain, outdir+cam + name_tags[ii] + '.gain', ' '

           endelse              ; count
        endfor                  ; ii
     endfor                     ; ipref
  endfor                        ; icam

end
