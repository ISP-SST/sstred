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
;      The number of threads to use for bad-pixel filling.
;   
;    descatter : in, optional, type=boolean 
;   
;      Do back-scatter compensation.
;   
;    ustat : in, optional, type=string
;   
;      Do only for this state.
;   
;    pref : in, optional, type=string
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
;-
pro red::sumpinh, nthreads = nthreads $
                  , descatter = descatter $
                  , ustat = ustat $
                  , pref = epref $
                  , pinhole_align = pinhole_align $
                  , brightest_only = brightest_only $
                  , lc_ignore = lc_ignore

  if ~keyword_set(pinhole_align) then pinhole_align = 0
  
  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  ;; Logging
  help, /obj, self, output = selfinfo 
  red_writelog, selfinfo = selfinfo
  
  if ~self.dopinh then begin
     print, inam + ' : ERROR : undefined pinh_dir'
     return
  endif
  if ~keyword_set(nthreads) then nthread = 2 else nthread = nthreads

  ;; Output dir
  if keyword_set(pinhole_align) then begin
     outdir = self.out_dir+ '/' + 'pinh_align/'
  endif else begin
     outdir = self.out_dir+ '/' + 'pinh/'
  endelse
  file_mkdir, outdir

  ;; Cameras
  cams = [self.camt, self.camr, self.camwb]
  Ncams = n_elements(cams)

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
     if n_elements(ustat) eq 0 then ustat = state[uniq(state, sort(state))]

     ;; Get camera tag
     cam = red_camtag(files[0])
     
     ;; Load darks
     dark = f0(self.out_dir + '/darks/'+cam+'.dark')

     ;; Loop
     firsttime = 1B
     ns = n_elements(ustat)
     for ii = 0L, ns-1 do BEGIN
        IF(n_elements(epref) GT 0) THEN BEGIN
           IF ((strsplit(ustat[ii],'.',/extract))[0] NE epref) THEN CONTINUE
        endif

        pref = (strsplit(ustat[ii], '.',/extract))[0]
        if cams[icam] eq self.camwb then begin
           flatf = self.out_dir + '/flats/' + strjoin([cam, pref, 'flat'],'.')
        endif else begin
           flatf = self.out_dir + '/flats/' + strjoin([cam, ustat[ii], 'flat'],'.')
        endelse

        if(file_test(flatf)) then begin
           flat = f0(flatf)
        endif else begin
           print, inam + ' : ERROR -> flat not found for '+cam+'.'+ustat[ii]
           stop
        endelse

        ;; Descatter data?
        if(keyword_set(descatter) AND self.dodescatter AND (pref eq '8542' OR pref eq '7772')) then begin
           if(firsttime) then begin
              
              ;; This code does not seem to know what prefilter we are
              ;; using! (Unless that knowledge is in descatter_dir.)

              descatter_psf_name =  self.descatter_dir+ '/' + cam + '.psf.f0'
              descatter_bgain_name = self.descatter_dir+ '/' + cam + '.backgain.f0'

;           ptf = self.descatter_dir+ '/' + tcam + '.psf.f0'
;           prf = self.descatter_dir+ '/' + rcam + '.psf.f0'
;           pwf = self.descatter_dir+ '/' + wcam + '.psf.f0'
;           btf = self.descatter_dir+ '/' + tcam + '.backgain.f0'
;           brf = self.descatter_dir+ '/' + rcam + '.backgain.f0'
;           bwf = self.descatter_dir+ '/' + wcam + '.backgain.f0'
              
              if file_test(descatter_psf_name) and file_test(descatter_bgain_name) then begin
                 Psft = f0(descatter_psf_name)
                 bgt = f0(descatter_bgain_name)
                 et = 1
              endif else et = 0
;           if(file_test(prf) and file_test(brf)) then begin
;              Psfr = f0(prf)
;              bgr = f0(brf)
;              er = 1
;           endif else er = 0
;           if(file_test(pwf) and file_test(bwf)) then begin
;              Psfw = f0(pwf)
;              bgw = f0(bwf)
;              ew = 1
;           endif else ew = 0
              
              firsttime = 0B
           endif

           if et then begin
              flat = red_cdescatter(flat, bgt, Psft, /verbose, nthreads = nthread)
           endif
        endif
        
        gain = red_flat2gain(flat)

        print, inam+' : summing pinh for state -> ' + ustat[ii]
        pos = where((state eq ustat[ii]), count)
        if count eq 0 then begin

           print, inam + ' : WARNING-> No files found for camera '+cams[icam]+' -> '+ ustat[ii]

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
              
              if(keyword_set(descatter) AND self.dodescatter AND (pref eq '8542' OR pref eq '7772')) then begin
                 if et then begin
                    c = red_cdescatter(c, bgt, Psft, /verbose, nthreads = nthread)
                 endif
              endif

              print, inam + ' : Filling pixels'
              c = red_fillpix(temporary(c) * gain, mask=gain ne 0, nthreads = nthread)

           endelse

           ;; Save
           head = 'n_aver=' + red_stri(count)

           namout = cam+'.' +ustat[ii]+'.pinh'
           print, inam+' : saving ' + outdir + namout
           fzwrite, fix(round(10. * c)), outdir+namout, head

           namout = cam+'.' +ustat[ii]+'.fpinh'
           print, inam+' : saving ' + outdir + namout
           fzwrite, c, outdir+namout, head

        endelse                 ; count
     endfor                     ; ii
  endfor                        ; icam

end
