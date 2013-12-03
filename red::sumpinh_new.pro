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
;    nthreads  : in, optional, type=integer
;   
;       The number of threads to use for backscattering.
;   
;    descatter :  in, optional, type=boolean
;   
;       Do descatter?
;   
;    ustat : in, optional, type=strarr
;   
;       Only for these states.
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
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2013-08-15 : MGL. Added keyword pinhole_align. Process cameras in
;                loops instead of separately. Added documentation.
; 
; 
; 
;-
pro red::sumpinh_new, nthreads = nthreads, descatter = descatter, ustat = ustat, pref=epref $
                  , pinhole_align = pinhole_align
  
  inam = 'red::sumpinh : '

  if ~self.dopinh then begin
     print, inam + 'ERROR : undefined pinh_dir'
     return
  endif

  if n_elements(nthreads) eq 0 then nthreads = 2 else nthreads = nthreads

  ;; Create file list for each camera
  spawn, 'find ' + self.pinh_dir + '/' + self.camt + "/ | grep im.ex |  grep -v '.lcd.'", tfiles
  spawn, 'find ' + self.pinh_dir + '/' + self.camr + "/ | grep im.ex |  grep -v '.lcd.'", rfiles
  spawn, 'find ' + self.pinh_dir + '/' + self.camwb + "/ | grep im.ex | grep -v '.lcd.'", wfiles
  nt = n_elements(tfiles)
  nr = n_elements(rfiles)
  nw = n_elements(wfiles)
  

  if (tfiles[0] eq '') OR (rfiles[0] eq '') OR (wfiles[0] eq '') then begin
     print, inam + 'files not found for all cameras:'
     print, '   '+self.camt+' -> '+ nt
     print, '   '+self.camr+' -> '+ nr
     print, '   '+self.camw+' -> '+ nw
     return
  endif
  ;;
  ;; sort files
  tfiles = red_sortfiles(tfiles)
  rfiles = red_sortfiles(rfiles)
  wfiles = red_sortfiles(wfiles)
  ;;
  ;; Get states for the cameras
  print, inam + 'extracting states for all cameras' 
  tstat = red_getstates(tfiles)
  rstat = red_getstates(rfiles)
  wstat = red_getstates(wfiles)
  stats = [wstat, tstat, rstat]

  ;;
  ;; Flagging first image after tuning (will use tstat for the states)
  ;;
  state = tstat.state
  if ~keyword_set(ustat) then ustat = state[uniq(state, sort(state))]
  
  ;;
  ;; Output dir
  ;;
  outdir = self.out_dir+ '/pinh/'
  file_mkdir, outdir

  ;;
  ;; Get camera tags
  ;;
  tcam = red_camtag(tfiles[0])
  rcam = red_camtag(rfiles[0])
  wcam = red_camtag(wfiles[0])
  cams = [wcam, tcam, rcam]
  Ncam = (size(cams, /dim))[0]

  xsz = self.camsz              ; CCD size in pixels
  ysz = self.camsz              ; CCD size in pixels

  darks = fltarr(xsz, ysz, Ncam)
  flats = fltarr(xsz, ysz, Ncam)
  gains = fltarr(xsz, ysz, Ncam)
  dnames = self.out_dir + '/darks/'+cams+'.dark' ; Dark file names
  fnames = strarr(Ncam)                          ; Flat file names

  ;;
  ;; Read darks
  ;;
  for icam = 0, Ncam-1 do begin
     if file_test(dnames[icam]) then begin
        darks[*, *, icam] = f0(dnames[icam])
     endif else begin
        print, inam + 'ERROR -> dark not found:'
        print, dnames[icam]
        stop
    endelse
  endfor

  ;;
  ;; Loop over states
  ;;
  ns = n_elements(ustat)
  for ii = 0L, ns-1 do begin
     
     pref = (strsplit(ustat[ii], '.',/extract))[0]
     fnames[0] = self.out_dir + 'flats/' + strjoin([cams[0], pref, 'flat'],'.')
     for icam = 1, Ncam-1 do begin
        fnames[icam] = self.out_dir + 'flats/' + strjoin([cams[icam], ustat[ii], 'flat'],'.')
     endfor

     ;; Boolean for if we want to do descatter for this prefilter
     DoDescatter = keyword_set(descatter) AND self.dodescatter AND (pref eq '8542' OR pref eq '7772')

     ;;
     ;; Read descatter data?
     ;;
     if DoDescatter then begin
        backscatter_psfs = fltarr(xsz, ysz, Ncam) ; Right size?
        backscatter_gain = fltarr(xsz, ysz, Ncam) ; Right size?
        for icam = 0, Ncam-1 do begin
           bpname = self.descatter_dir+ '/' + cams[icam] + '.psf.f0'
           bgname = self.descatter_dir+ '/' + cams[icam] + '.backgain.f0'
           if file_test(bpname) and file_test(bgname) then begin
              backscatter_psfs[*, *, icam] = f0(bpname)
              backscatter_gain[*, *, icam] = f0(bgname)
           endif else begin
              print, inam + 'ERROR -> backscatter psf and/or gain not found:'
              print, bpname
              print, bgname
              stop    
           endelse
        endfor
     endif

     ;;
     ;; Read flats 
     ;;
     for icam = 0, Ncam-1 do begin
        fname = self.out_dir + '/darks/'+cams[icam]+'.dark'
        if file_test(fnames[icam]) then begin
           flats[*, *, icam] = f0(fnames[icam])
        endif else begin
           print, inam + 'ERROR -> flat not found:'
           print, fnames[icam]
           stop
        endelse
        if DoDescatter then begin
           flats[*, *, icam] = red_cdescatter(flats[*, *, icam] $
                                              , backscatter_gain[*, *, icam] $
                                              , backscatter_psfs[*, *, icam] $
                                              , /verbose, nthreads = nthreads)
        endif
        gains[*, *, icam] = red_flat2gain(flats[*, *, icam]) 
     endfor

     ;; Summing pinholes
     pnames = cams+'.' +ustat[ii]+'.pinh'
     pfnames = cams+'.' +ustat[ii]+'.fpinh'
     for icam = 0, Ncam-1 do begin

        if file_test(outdir+pnames[icam]) and file_test(outdir+pfnames[icam]) then begin

           print, inam + 'summed pinhole files already exists:'
           print, pnames[icam], ' and ', pfnames[icam]

        endif else begin

           print, inam + 'make summed pinhole file:'
           print, pnames[icam], ' and ', pfnames[icam]
           
           pos = where((stats[icam].state eq ustat[ii]), count)
           if count eq 0 then begin

              print, inam + 'WARNING-> No files found for the '+cams[icam]+' -> ' + ustat[ii]

           endif else begin
              
              if DoDescatter then begin
                 ;; Let red_sumfiles do the descattering
                 psum = red_sumfiles(stats[icam].files[pos], pinhole_align = keyword_set(pinhole_align) $
                                     , gain = gains[*, *, icam], dark = darks[*, *, icam] $
                                     , backscatter_gain = backscatter_gain[*, *, icam] $
                                     , backscatter_psfs = backscatter_psfs[*, *, icam] $
                                     , nthreads = nthreads)  
;                 psum = red_cdescatter(psum $
;                                       , backscatter_gain[*, *, icam] $
;                                       , backscatter_psfs[*, *, icam] $
;                                       , /verbose, nthreads = nthreads)
              endif else begin
                 psum = red_sumfiles(stats[icam].files[pos], pinhole_align = keyword_set(pinhole_align) $
                                     , gain = gains[*, *, icam], dark = darks[*, *, icam])  
              endelse
              
              ;;
              ;; Save
              ;;
              head = 'n_aver=' + red_stri(count)
              print, inam + 'saving ' + outdir + pnames[icam]
              fzwrite, fix(round(10. * psum)), outdir+pnames[icam], head
              print, inam + 'saving ' + outdir + pfnames[icam]
              fzwrite, psum, outdir+pfnames[icam], head
           endelse

        endelse

     endfor

  endfor
end
