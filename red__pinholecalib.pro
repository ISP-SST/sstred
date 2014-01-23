; docformat = 'rst'

;+
; Do momfbd pinhole calibration from within IDL.
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, Institute for Solar Physics, 2012
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
;    state : in, optional
;
;       Tuning and LC state (If not given, use brightest)
;
;    prefilter : in, optional 
;
;       Run only on this prefilter.
;
;    diversity : in, out
;
;       Diversities (in and out)
;
;    finddiversity : in, optional
;
;       Set this to calibrate diversity
;
;    maxdiversity : in, optional 
;
;       Diversity cannot become larger than this [mm]
;
;    deltadiversity : in, optional 
;
;       Search steps
;
;    test_diversities : in, optional  
;
; 
;
;    test_metrics : in, optional  
; 
;    nslaves : in, optional  
;
;        For momfbd
;
;    nthreads : in, optional  
;
;        For momfbd
;
;    xconvergence : out, optional  
;
;        X tilt convergence
;
;    yconvergence : out, optional  
;
;        Y tilt convergence
;
;    dconvergence : out, optional   
;
;        Diversity convergence
;
;    dvalues : out, optional    
;
;        Diversity convergence
;
;    mconvergence : out, optional     
;
;        Metric convergence (out)
;
;    maxit : in, optional   
;
;        Max # of iterations
;
;    mtol : in, optional   
;
;        Metric tolerance
;
;    dtol : in, optional   
;
;        Diversity tolerance
;
;    ttol : in, optional   
;
;        Tilt tolerance
;
;    strehl : out, optional      
;
;        Strehl number
; 
; :History:
; 
;   2013-09-04 : MGL. Use red_momfbd_setup, not momfbd_setup. Improve
;                the procedure for measuring the residual dark level.
;                Now works with median and stdev instead of
;                biweight_mean and robust_sigma. Implement avoiding
;                certain patch sizes because they crash the momfbd
;                slaves.
; 
;   2013-09-06 : MGL. Warning and stop when keyword finddiversity is
;                set. This part is not ready for common use. Removed
;                keyword cams. Do initialization of the offsets, so we
;                don't have to do the red::getoffsets step anymore.
;
;   2013-09-10 : MGL. Moved writing of subfield images to disk (and
;                the associated dark level refinement) to
;                red_pinh_run_momfbd. 
;
;   2013-09-11 : MGL. Use red_findmax2qi, not findmax2qi. Also
;                red_stats rather than stats.
;
;   2014-01-23 : MGL. Use red_extractstates instead of red_getstates
;                and local extraction of info from file names.
;
;-
pro red::pinholecalib, STATE = state $                   
                       , PREFILTER = prefilter $         
                       , DIVERSITY = diversity $         
                       , FINDDIVERSITY = finddiversity $ 
                       , MAXDIVERSITY = maxdiversity $   
                       , DELTADIVERSITY = deltadiversity $
                       , TEST_DIVERSITIES = test_diversities $ 
                       , TEST_METRICS = test_metrics $ 
                       , NSLAVES = nslaves $           
                       , NTHREADS = nthreads $         
                       , XCONVERGENCE = xconvergence $ 
                       , YCONVERGENCE = yconvergence $ 
                       , DCONVERGENCE = dconvergence $ 
                       , DVALUES = dvalues $           
                       , MCONVERGENCE = mconvergence $ 
                       , MAXIT = maxit $               
                       , MTOL = mtol $                 
                       , DTOL = dtol $                 
                       , TTOL = ttol $                 
                       , STREHL = strehl
  
  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])
  
  ;; Logging
  help, /obj, self, output = selfinfo 
  red_writelog, selfinfo = selfinfo, logfile = logfile

  if n_elements(mtol) eq 0 then mtol = 1e-6 ; Metric tolerance
  if n_elements(dtol) eq 0 then dtol = 1e-6 ; Diversity tolerance
  if n_elements(ttol) eq 0 then ttol = 1e-3 ; Tilt tolerance
  
  if n_elements(nslaves) eq 0 then nslaves=6
  if n_elements(nthreads) eq 0 then nthreads=1
  if n_elements(maxit) eq 0 then maxit=400


  ;; Where the summed pinhole images are:
  pinhdir = self.out_dir + '/pinh/'

  ;; Calib directory
  calibdir = self.out_dir + '/calibnew/'
  if ~file_test(calibdir) then begin
     file_mkdir, calibdir
     ;; Use initial offset files from ordinary calib dir.
;     spawn, 'cp -f '+ self.out_dir + '/calib/*offs ' + calibdir
  endif


  ;; Make temporary work directory, as well as momfbd "data" dir
  workdir = calibdir+'work/'
  file_mkdir, workdir
  file_mkdir, workdir+'data'

  ;; Get camera labels:
  self -> getcamtags, dir = self.pinh_dir 
  camt = self.camttag
  camr = self.camrtag
  camw = self.camwbtag
  cams = [camw, camt, camr] 
  
  xsz = self.camsz              ; CCD size in pixels
  ysz = self.camsz              ; CCD size in pixels

  ;; Did we select prefilter?
  if n_elements(prefilter) eq 0 then begin
     print, 'No prefilter given.' ;    <----------------------------
     retall
  endif
  ;; If not, should we find all prefilters and loop over them? YES

  ;; Get list of states (and transmitted file names)  
  cam = camt                    ; Transmitted camera
  files = file_search(pinhdir+cam+'.'+prefilter+'*.fpinh', count = nf)

  
  ;; These file names have less "fields" than the raw files, so
  ;; we add one field when using red_getstates.
;  states = red_getstates(files+'.dum')
;  ustat = states.state[uniq(states.state, sort(states.state))]
;  red_extractstates, files, /basename, fullstate = fullstate
;  ustat = fullstate[uniq(fullstate, sort(fullstate))]
;  ns = n_elements(ustat)

  ;; We need the total state information of the files we found.
  ;; The number of fields depends on the way the pinholes were summed.
  ;; Everything between the camera name and the .fpinh extension
  ;; comprise this information.
  ustat = red_strreplace(red_strreplace(file_basename(files), cam+'.', ''), '.fpinh', '')
  ns = n_elements(ustat)

;  fullstate = strarr(nf)
;  for ii = 0L, nf -1 do begin
;     tmp = strsplit(file_basename(files[ii]), '._', /extract)
;     n = n_elements(tmp)
;     wav = tmp[2]+'_'+tmp[3]
;     pref = tmp[1]
;     lc = tmp[4]
;     fullstate[ii] = pref+'.'+wav + '.' + lc
;  endfor
;  ustat = fullstate[uniq(fullstate, sort(fullstate))]
;  ns = n_elements(ustat)

  ;; We want to process all existing pinhole images. Unless the
  ;; keyword state is provided, then we want to process only the ones
  ;; matching this keyword.

  if n_elements(state) ne 0 then begin

     ;; Check that the state given in keywords exists. Quit with error
     ;; if it does not.

     if total(strmatch(ustat,state)) eq 0 then begin
        print, inam+' : state not found: '+state
        stop
     endif else begin
        ustat = state
     endelse

  endif 

  for istat = 0, n_elements(ustat)-1 do begin

     stat = ustat[istat]

     ;; File names for the full-FOV pinhole images and offsets.
     pnames = cams + '.' + stat + '.fpinh' 
     xnames = cams + '.' + stat + '.xoffs'
     ynames = cams + '.' + stat + '.yoffs'
     
     ;; File name templates for the subfield-size images and offsets.
     ptemplates = cams + '.' + stat + '.pinh.%03d' 
     xtemplates = cams + '.' + stat + '.xoffs.%03d'
     ytemplates = cams + '.' + stat + '.yoffs.%03d'
     
     Nch = (size(pnames, /dim))[0]
     if n_elements(diversity) eq 0 then diversity = fltarr(Nch)

     ;; Get align_clips and apply them to images. 
     clipfile = self.out_dir + '/calib/align_clips.'+prefilter+'.sav'
     IF(~file_test(clipfile)) THEN BEGIN
        print, inam + ' : ERROR -> align_clip file not found'
        print, inam + ' : -> you must run red::getalignclips first!'
        return
     endif
     restore, clipfile             

     ;; The clipfile provides:      
     ;; ACL             STRING    = Array[3] - CLIP strings for config file
     ;; CL              INT       = Array[4, 3] - CLIPs as numbers
     ;; REFROT          INT       = Rot orientation code
     ;; SSH             INT       = Array[2, 3] ????
     ;; SX              LONG      = clipped X size
     ;; SY              LONG      = clipped Y size
     ;; 
     ;; Put this info in the log file:
     clipinfo = strarr(7)
     clipinfo[0] = 'Clip info:'
     clipinfo[1:3] = ACL
     clipinfo[4] = 'refrot: '+red_stri(refrot)
     clipinfo[5] = 'sx: '+red_stri(sx)
     clipinfo[6] = 'sy: '+red_stri(sy)
     red_writelog, /add, logfile = logfile, top_info_strings = clipinfo



     ;; Read summed pinhole array images for the selected state. These
     ;; images are already corrected for dark and flat as well as
     ;; fillpixed by the sumfiles routine. First image in array is the
     ;; "anchor". Also read the existing offsets. Apply clips.
     images = fltarr(sx, sy, Nch)

     for ich = 0, Nch-1 do begin
        print, 'Reading and clipping ' + pinhdir + pnames[ich]
        images[*, *, ich] = red_clipim(f0(pinhdir + pnames[ich]), cl[*, ich])
     endfor

     ;; Find pinhole grid.
     red_findpinholegrid_new, images[*, *, 0], simx, simy, dx = dx, dy = dy, Npinh = Npinh
     ;; Set subimage size to mean distance between spots.
     sz = round((dx+dy)/2.)
     ;; Make sure sz is an even number:
     sz += (sz mod 2)

     ;; Avoid sizes that make the momfbd slaves crash.
     badsizes = [90]
     while total(sz eq badsizes) gt 0 do sz += 2

     
     ;; Measure positions used for calculating initial offsets.
     d = sz/2
     strehl = fltarr(Nch, Npinh)
     pos = fltarr(2, Nch, Npinh) ; Peak positions within subfield
     for ihole=0,Npinh-1 do begin

        ;; Read out subfield (ix,iy) from images
        subfields = images[simx[ihole]-d:simx[ihole]+d-1,simy[ihole]-d:simy[ihole]+d-1, *]
        
        for ich = 0, Nch-1 do begin

           strehl[ich, ihole] = max(subfields[*, *, ich])/mean(subfields[*, *, ich])
           
           ;; Peak position
           maxloc = (where(subfields[*,*,ich] eq max(subfields[*,*,ich])) )[0]
           pos[*, ich, ihole] = [maxloc mod sz, maxloc/sz]
           
           print, 'pos', [maxloc mod sz, maxloc/sz]

        endfor                  ; ich
     endfor                     ; ihole

     
                                ; oldxoffs = fltarr(sx, sy, Nch)
                                ; oldyoffs = fltarr(sx, sy, Nch)
                                ; oldcalibdir = self.out_dir + '/calib/'
 ;;; Read the existing offset files if present
                                ; for ich = 1, Nch-1 do begin
                                ;    if file_test(oldcalibdir + xnames[ich]) then begin
                                ;       print, 'Reading  and clipping' + oldcalibdir + xnames[ich]
                                ;       oldxoffs[*, *, ich] = f0(oldcalibdir + xnames[ich])
                                ;    endif
                                ;    if file_test(oldcalibdir + ynames[ich]) then begin
                                ;       print, 'Reading  and clipping' + oldcalibdir + ynames[ich]
                                ;       oldyoffs[*, *, ich] = f0(oldcalibdir + ynames[ich])
                                ;    endif
                                ; endfor

     ;; Make the tilts in units of 1/100 pixel.
     xtilts = fltarr(Nch, Npinh)
     xtilts[1, *] = (pos[0, 1, *] - pos[0, 0, *])*100
     xtilts[2, *] = (pos[0, 2, *] - pos[0, 0, *])*100

     ytilts = fltarr(Nch, Npinh)
     ytilts[1, *] = (pos[1, 1, *] - pos[1, 0, *])*100
     ytilts[2, *] = (pos[1, 2, *] - pos[1, 0, *])*100

;; Test!
;  xtilts += +300
;  ytilts += -200

     ;; Make fits to the "tilts" in units of 1/100 pixel.
     red_pinh_make_fits, simx, simy, sx, sy $
                         , xtilts = xtilts $ ; Input
                         , ytilts = ytilts $ ; Input
                         , dxoffs = xoffs $  ; Output
                         , dyoffs = yoffs    ; Output
stop     
;   help, oldxoffs, xoffs 
;   cgplot,oldxoffs[*,*,2],xoffs[*,*,2],psym=3, /aspect
;   cgplot,/over,color='red',[-1000,1000],[-1000,1000]                    
; stop  
     ;; Write the initial offsets to files
     for ich = 1, Nch-1 do begin
        print, 'Writing ' + calibdir + xnames[ich]+'.init'
        fzwrite, fix(round(xoffs[*, *, ich])), calibdir + xnames[ich]+'.init', ''
        print, 'Writing ' + calibdir + ynames[ich]+'.init'
        fzwrite, fix(round(yoffs[*, *, ich])), calibdir + ynames[ich]+'.init', ''
     endfor

     
     ;; Start iterations
     ;; Set up manager and slaves.
     red_momfbd_setup, PORT = port, NSLAVES = nslaves, NMANAGERS=1, /FREE, NTHREADS=nthreads

     xconvergence = fltarr(Nch, MaxIt)
     yconvergence = fltarr(Nch, MaxIt)
     mconvergence = fltarr(Npinhx, Npinhy, MaxIt)
     if keyword_set(finddiversity) then begin
        dconvergence = fltarr(Nch, MaxIt)
        dvalues = fltarr(Nch, MaxIt)
     endif

     ;; Iterate tilts only
     for it = 0, MaxIt-1 do begin

        ;; Run momfbd on the pinholes and read back the results
        red_pinh_run_momfbd, images, xoffs, yoffs, simx, simy, sz $
                             , xtemplates, ytemplates, ptemplates $
                             , stat $
                             , self.telescope_d, self.image_scale, self.pixel_size $
                             , DIVERSITY = diversity $
                             , FINDDIVERSITY = 0 $
                             , WORKDIR = workdir $
                             , PORT = port $
                             , NSLAVES = nslaves $
                             , NTHREADS = nthreads $
                             , XTILTS = dxtilts $ ; Output
                             , YTILTS = dytilts $ ; Output
                             , METRICS = metrics  ; Output

        ;; Updates
        xtilts += dxtilts
        ytilts += dytilts

        oldxoffs = xoffs
        oldyoffs = yoffs

        ;; Make fits to the tilts 
        red_pinh_make_fits, simx, simy, sx, sy $
                            , XTILTS = xtilts $ ; Input
                            , YTILTS = ytilts $ ; Input
                            , dxoffs = xoffs $  ; Output
                            , dyoffs = yoffs    ; Output
        
;     ;; Updates
;     xoffs += dxoffs
;     yoffs += dyoffs

        dxoffs = xoffs-oldxoffs
        dyoffs = yoffs-oldyoffs

        red_stats, xoffs, name = 'xoffs'
        red_stats, dxoffs, name = 'dxoffs'

        red_stats, yoffs, name = 'yoffs'
        red_stats, dyoffs, name = 'dyoffs'


        ;; Convergence
        for ich = 1, Nch-1 do begin
           xconvergence[ich, it] = max(dxoffs[*, *, ich])-min(dxoffs[*, *, ich])
           yconvergence[ich, it] = max(dyoffs[*, *, ich])-min(dyoffs[*, *, ich])
        endfor
        mconvergence[*, *, it] = metrics

        red_pinh_plot_convergence, IT = it $ ; # of iterations so far
                                   , DTOL = dtol $
                                   , TTOL = ttol $
                                   , XCONVERGENCE = xconvergence $ ; X tilt convergence (out)
                                   , YCONVERGENCE = yconvergence $ ; Y tilt convergence (out)
                                   , MCONVERGENCE = mconvergence   ; Metric convergence (out)     
        
        ;; If converged, then get out of the it loop.
        if max(abs(xconvergence[*, it])) lt ttol and $
           max(abs(yconvergence[*, it])) lt ttol then break
        if it gt 0 then begin
           dm = max(mconvergence[*, *, it]-mconvergence[*, *, it-1])
           sm = max(mconvergence[*, *, it]+mconvergence[*, *, it-1])
           if dm/sm lt mtol then break
        endif

     endfor                     ; for it
     

     if keyword_set(finddiversity) then begin

        print, inam+' : Do you want to refine the diversity?'
        print, inam+' : '
        print, inam+' : This part is not well tested because we currently do not have'
        print, inam+' : a CRISP diversity camera. It was written for the purpose of  '
        print, inam+' : finding small diversity differences between nominally focused'
        print, inam+' : cameras. We never were able to make this work in a way we    '
        print, inam+' : trusted.                                                     '
        print, inam+' : '
        print, inam+' : Be prepared to do some debugging work or possibly rewriting  '
        print, inam+' : this functionality completely. /MGL '

        stop

        ;; Try some combinations of diversities, find min(metric).
        makegrid = n_elements(test_diversities) eq 0


        if makegrid then begin
           if n_elements(deltadiversity) eq 0 then deltadiversity = .5 ; focus step length [mm]
           if n_elements(maxdiversity) eq 0 then Nf = 10 else Nf = maxdiversity/deltadiversity 
           test_diversities = fltarr(Nch, Nf+1, 2*Nf+1)
           for ii = 0, Nf do begin
              for jj = -Nf, Nf do begin
                 test_diversities[*, ii, jj+Nf] = diversity + [0., ii*deltadiversity, jj*deltadiversity]
              endfor
           endfor
           domirror = 1
        endif else domirror = 0

        Nf0 = (size(test_diversities, /dim))[1]
        Nf1 = (size(test_diversities, /dim))[2]
        test_metrics = fltarr(Nf0, Nf1)

        ;; This loop assumes there are only two diversities to vary!
        for ii = 0, Nf0-1 do begin
           for jj = 0, Nf1-1 do begin
              
              print, 'Test diversity ', test_diversities[*, ii, jj], ' mm.'
              

              ;; Run momfbd on the pinholes and read back the results
              red_pinh_run_momfbd, images, xoffs, yoffs, simx, simy, sz $
                                   , xtemplates, ytemplates, ptemplates $
                                   , stat $
                                   , self.telescope_d, self.image_scale, self.pixel_size $
                                   , DIVERSITY = test_diversities[*, ii, jj] $
                                   , FINDDIVERSITY = 0 $
                                   , WORKDIR = workdir $
                                   , PORT = port $
                                   , NSLAVES = nslaves $
                                   , NTHREADS = nthreads $
                                   , FOC = foc $    ; Output
                                   , METRICS = metrics ; Output

              test_metrics[ii, jj] = median(metrics)

              print, test_metrics

           endfor
        endfor

        if domirror then begin
           ;; Exploit mirror symmetry
           test_metrics2 = [reverse(test_metrics[1:*,*],1),test_metrics]
           ss = size(test_metrics2, /dim)
           test_diversities2 = fltarr([Nch, ss])
           test_diversities2[*, 0:Nf0-2, *] = reverse(test_diversities[*, 1:*,*], 2)
           test_diversities2[*, Nf0-1:*, *] = test_diversities
        endif else begin
           test_metrics2 = test_metrics
           test_diversities2 = test_diversities
           ss = size(test_metrics2, /dim)
        endelse
        
        
        tmp = min(test_metrics2, minloc)
        ii = minloc MOD ss[0]
        jj = minloc / ss[0]
        print, 'Min metric: ', test_metrics2[ii, jj], ' at ', [ii, jj]
        best_diversity = fltarr(Nch)
        for ich = 1, Nch-1 do best_diversity[ich] = test_diversities2[ich, ii, jj]
        print, 'Best diversity: ', best_diversity
        
        if 1 then begin
           ;; Improve best_diversity by interpolation before considering
           ;; them "found".
           c = test_metrics2[ii-1:ii+1, jj-1:jj+1]
           aaa = red_findmax2qi(-c, /verbose) ; "subpixel" minimum
           cmin = Interpolate(c, aaa[0], aaa[1], cubic = -0.5)
           print, 'Subpixel min metric: ', cmin, ' at ', aaa
           for ich = 1, Nch-1 do best_diversity[ich] $
              = Interpolate(reform(test_diversities2[ich, ii-1:ii+1, jj-1:jj+1]) $
                            , aaa[0], aaa[1], cubic = -0.5)
           
           print, 'Best subpixel diversity: ', best_diversity
           ;;Best subpixel diversity:       0.00000     0.784720     -3.08438
        endif 

        ;; Use the found diversities
        diversity = best_diversity

        
        ;; Now iterate to finalize the offsets
        for it = 0, MaxIt-1 do begin

           ;; Run momfbd on the pinholes and read back the results
           red_pinh_run_momfbd, images, xoffs, yoffs, simx, simy, sz $
                                , xtemplates, ytemplates, ptemplates $
                                , stat $
                                , self.telescope_d, self.image_scale, self.pixel_size $
                                , DIVERSITY = diversity $
;                             , FINDDIVERSITY = finddiversity $
                                , WORKDIR = workdir $
                                , PORT = port $
                                , NSLAVES = nslaves $
                                , NTHREADS = nthreads $
                                , XTILTS = xtilts $ ; Output
                                , YTILTS = ytilts $ ; Output
                                , FOC = foc $       ; Output
                                , METRICS = metrics ; Output

           ;; Make fits to the tilts and average the diversity
           red_pinh_make_fits, simx, simy, sx, sy $
                               , XTILTS = xtilts $ ; Input
                               , YTILTS = ytilts $ ; Input
;                            , FOC = foc $       ; Input
;                            , dfoc = dfoc $     ; Output
;                            , ddiv = ddiv $       ; Output
                               , dxoffs = dxoffs $ ; Output
                               , dyoffs = dyoffs   ; Output
           

           ;; Updates
           xoffs += dxoffs
           yoffs += dyoffs
;        diversity += -ddiv

           ;; Convergence
           for ich = 1, Nch-1 do begin
              xconvergence[ich, it] = max(dxoffs[*, *, ich])-min(dxoffs[*, *, ich])
              yconvergence[ich, it] = max(dyoffs[*, *, ich])-min(dyoffs[*, *, ich])
           endfor
;        dconvergence[*, it] = ddiv
;        dvalues[*, it] = diversity
           mconvergence[*, *, it] = metrics


;        if keyword_set(finddiversity) then begin
;           red_pinh_plot_convergence, IT = it $ ; # of iterations so far
;                                      , DTOL = dtol $
;                                      , TTOL = ttol $
;                                      , XCONVERGENCE = xconvergence $ ; X tilt convergence (out)
;                                      , YCONVERGENCE = yconvergence $ ; Y tilt convergence (out)
;                                      , DCONVERGENCE = dconvergence $ ; Diversity convergence (out)
;                                      , DVALUES = dvalues $      ; Diversity convergence (out)
;                                      , MCONVERGENCE = mconvergence ; Metric convergence (out)     
;        endif else begin
           red_pinh_plot_convergence, IT = it $ ; # of iterations so far
                                      , DTOL = dtol $
                                      , TTOL = ttol $
                                      , XCONVERGENCE = xconvergence $ ; X tilt convergence (out)
                                      , YCONVERGENCE = yconvergence $ ; Y tilt convergence (out)
                                      , MCONVERGENCE = mconvergence   ; Metric convergence (out)     
;        endelse

           ;; If converged, then get out of the it loop.
;        if keyword_set(finddiversity) then begin
;           if max(abs(xconvergence[*, it])) lt ttol and $
;              max(abs(yconvergence[*, it])) lt ttol and $
;              max(abs(dconvergence[*, it])) lt dtol then break
;        endif else begin
;           if max(abs(xconvergence[*, it])) lt ttol and $
;              max(abs(yconvergence[*, it])) lt ttol then break
;        endelse

           if it gt 0 then begin
              dm = max(mconvergence[*, *, it]-mconvergence[*, *, it-1])
              sm = max(mconvergence[*, *, it]+mconvergence[*, *, it-1])
              if dm/sm lt mtol then break
           endif

        endfor                  ; for it

     endif                      ; finddiversity
     

     ;; Kill manager and slaves
     red_momfbd_setup, PORT = port, NMANAGERS=0
     ;;  spawn,'rm -f '+workdir

     ;; Crop convergence arrays.
     xconvergence = xconvergence[*, 0:it-1]
     yconvergence = yconvergence[*, 0:it-1]
     if keyword_set(finddiversity) then dconvergence = dconvergence[*, 0:it-1]
     

     ;; Write the updated offsets to files
     for ich = 1, Nch-1 do begin
        print, 'Writing ' + calibdir + xnames[ich]
        fzwrite, fix(round(xoffs[*, *, ich])), calibdir + xnames[ich], ''
        print, 'Writing ' + calibdir + ynames[ich]
        fzwrite, fix(round(yoffs[*, *, ich])), calibdir + ynames[ich], ''
     endfor

     fzwrite, diversity, calibdir+'diversity.fz', ''

  endfor                        ; istat

end
