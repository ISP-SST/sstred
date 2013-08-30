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
;    cams : in, optional    
;
;        Cam tags in correct order
; 
; 
; :History:
; 
; 
; 
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
                       , STREHL = strehl $             
                       , CAMS = cams                   
  
  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])
  
  ;; Logging (move to end of file so we can log the clipfile?)
  help, /obj, self, output = selfinfo 
  red_writelog, selfinfo = selfinfo

  if n_elements(mtol) eq 0 then mtol = 1e-6 ; Metric tolerance
  if n_elements(dtol) eq 0 then dtol = 1e-6 ; Diversity tolerance
  if n_elements(ttol) eq 0 then ttol = 1e-3 ; Tilt tolerance
  
  if n_elements(nslaves) eq 0 then nslaves=6
  if n_elements(nthreads) eq 0 then nthreads=1


  ;; Where the summed pinhole images are:
  pinhdir = self.out_dir + '/pinh_test/'

  ;; Calib directory
  calibdir = self.out_dir + '/calib_test/'
  if ~file_test(calibdir) then begin
     file_mkdir, calibdir
     ;; Use initial offset files from ordinary calib dir.
     spawn, 'cp -f '+self.out_dir + '/calib/*offs ' + calibdir
  endif


  ;; Make temporary work directory, as well as momfbd "data" dir
  workdir = calibdir+'work/'
  file_mkdir, workdir
  file_mkdir, workdir+'data'

  if n_elements(cams) eq 0 then begin
     ;; Something like this to get camera labels:
     self -> getcamtags, dir = self.pinh_dir 
     cams = [self.camwbtag, self.camttag, self.camrtag] 
  endif

  xsz = self.camsz              ; CCD size in pixels
  ysz = self.camsz              ; CCD size in pixels

  ;; Did we select prefilter?
  if n_elements(prefilter) eq 0 then begin
     print, 'No prefilter given.' ;    <----------------------------
     retall
  endif
  ;; If not, should we find all prefilters and loop over them?

  ;; Get list of states (and transmitted file names)  
  cam = cams[1]                   ; Transmitted camera
  spawn, 'cd ' + pinhdir + '; find  | cut -d/ -f2 | grep "'+cam+'." | grep ".pinh" | grep '+prefilter, files
  nf = n_elements(files)


  ;; These file names have less "fields" than the raw files, so
  ;; we can't use red_getstates. But we'll do something similar.
  fullstate = strarr(nf)
  for ii = 0L, nf -1 do begin
     tmp = strsplit(file_basename(files[ii]), '._', /extract)
     n = n_elements(tmp)
     wav = tmp[2]+'_'+tmp[3]
     pref = tmp[1]
     lc = tmp[4]
     fullstate[ii] = pref+'.'+wav + '.' + lc
  endfor
  ustat = fullstate[uniq(fullstate, sort(fullstate))]
  ns = n_elements(ustat)

  if n_elements(state) ne 0 then begin

     ;; Check that the state given in keywords exists. Quit with error
     ;; if it does not.

     if total(strmatch(ustat,state)) eq 0 then begin
        print, inam+' : state not found: '+state
        stop
     endif else begin
        stat = state
     endelse

  endif else begin

     ;; Pick one of the brightest states to work with. (Assuming we
     ;; don't want to perform the calibration for all states
     ;; independently.)

     maxvals = fltarr(ns)
     for ii = 0L, ns-1 do begin
        fname = pinhdir + cam + '.' + ustat[ii] + '.fpinh'
        maxvals[ii] = max(f0(fname))
     endfor
     tmp = max(maxvals, mloc)
     stat = ustat[mloc]

  endelse


  ;; File names for the full-FOV pinhole images and offsets.
  pnames = cams + '.' + stat + '.fpinh' 
  xnames = cams + '.' + stat + '.xoffs'
  ynames = cams + '.' + stat + '.yoffs'
  
  ;, File name templates for the subfield-size images and offsets.
  ptemplates = cams + '.fpinh.%03d' 
  xtemplates = cams + '.xoffs.%03d'
  ytemplates = cams + '.yoffs.%03d'

  Nch = (size(pnames, /dim))[0]
  if n_elements(diversity) eq 0 then diversity = fltarr(Nch)

  ;; Get align_clips and apply them to images. 
  clipfile = self.out_dir + '/calib/align_clips.'+pref+'.sav'
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
  findpinholegrid, images[*, *, 0], simx, simy
  simx = simx[1:(size(simx, /dim))[0]-2] ; Skip outermost rows and columns?
  simy = simy[1:(size(simx, /dim))[0]-2]
  ;; Set subimage size to mean distance between spots:
  sz = round(median([deriv(simx),deriv(simy)])) ;>128
  ;; Make sure sz is an even number:
  sz += (sz mod 2)



  ;; Remove rows and columns that are too close to the border
  simx = simx[where(simx gt sz/2 and simy lt sx-sz/2-1)]
  simy = simy[where(simy gt sz/2 and simy lt sx-sz/2-1)]
  Npinhx = (size(simx, /dim))[0]
  Npinhy = (size(simy, /dim))[0]


  ;; Preprocess pinholes
  imno = 0
  d=sz/2
  strehl = fltarr(Npinhx, Npinhy, Nch)
  for ix=0,Npinhx-1 do begin
     for iy=0,Npinhy-1 do begin

        ;; Read out subfield (ix,iy) from images
        subfields = images[simx[ix]-d:simx[ix]+d-1,simy[iy]-d:simy[iy]+d-1, *]

        ;; Fine tune dark level and normalization for each subfield.
        ;; Get the bias by fitting a Gaussian to the histogram peak.
        ;; Then write to files with imno in name
        Nbins = 10000
        for ich = 0, Nch-1 do begin
           
           fname = workdir+string(format = '(%"'+ptemplates[ich]+'")', imno)

           if ~file_test(fname) or arg_present(strehl) then begin

              subfields[*, *, ich] = subfields[*, *, ich]/max(subfields[*, *, ich])

              mn = biweight_mean(subfields[*,*,ich])
              st = robust_sigma(subfields[*,*,ich])

              hmax = mn + 400.*st
              hmin = mn - 400.*st

              hh = histogram(subfields[*, *, ich], min = hmin, max = hmax $
                             , Nbins = Nbins, locations = locations)
              binsize = (hmax - hmin) / (Nbins - 1)
              intensities = locations + binsize/2.
              
              plot,intensities,hh,/xstyle,xrange=[-.001,.001], /ystyle
              
              print,max(smooth(hh,15),ml)  
              Nfit = 51
              indx = ml + indgen(Nfit) - (Nfit-1)/2
              
              plot,intensities[indx],hh[indx], /XSTYLE, /YSTYLE
              
              yfit=MPFITPEAK(float(intensities[indx]), float(hh[indx]), a, Nterms = 5)
              oplot,intensities[indx],yfit,color=fsc_color('yellow')
              oplot,[0,0]+a[1],[0,1]*max(hh),color=fsc_color('yellow')
              oplot,[0,0]+mn,[0,1]*max(hh),color=fsc_color('cyan')
              print,intensities[ml]  
              
              ;; Subtract the fitted dark level and renormalize to unit
              ;; max. The momfbd program will then renormalize all
              ;; channels to the same total energy.
              subfields[*, *, ich] = subfields[*, *, ich] - a[1]
              subfields[*, *, ich] = subfields[*, *, ich]/max(subfields[*, *, ich])

              fzwrite, fix(round(subfields[*, *, ich]*15000.)), fname, ''

              strehl[ix, iy, ich] = max(subfields[*, *, ich])/mean(subfields[*, *, ich])

           endif                ; if file_test
        endfor                  ; for ich

        imno += 1

     endfor                     ; for iy
  endfor                        ; for ix

  xoffs = fltarr(sx, sy, Nch)
  yoffs = fltarr(sx, sy, Nch)
  ;; Read the existing offset files if present
  for ich = 1, Nch-1 do begin
     if file_test(calibdir + xnames[ich]) then begin
        print, 'Reading  and clipping' + calibdir + xnames[ich]
        xoffs[*, *, ich] = f0(calibdir + xnames[ich])
     endif
     if file_test(calibdir + ynames[ich]) then begin
        print, 'Reading  and clipping' + calibdir + ynames[ich]
        yoffs[*, *, ich] = f0(calibdir + ynames[ich])
     endif
  endfor

  ;; Start iterations
  ;; Set up manager and slaves.
  momfbd_setup, PORT = port, NSLAVES = nslaves, NMANAGERS=1, /FREE, NTHREADS=nthreads

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
     red_pinh_run_momfbd, xoffs, yoffs, simx, simy, sz $
                          , xtemplates, ytemplates, ptemplates $
                          , stat $
                          , self.telescope_d, self.image_scale, self.pixel_size $
                          , DIVERSITY = diversity $
                          , FINDDIVERSITY = 0 $
                          , WORKDIR = workdir $
                          , PORT = port $
                          , NSLAVES = nslaves $
                          , NTHREADS = nthreads $
                          , XTILTS = xtilts $           ; Output
                          , YTILTS = ytilts $           ; Output
                          , METRICS = metrics           ; Output

     ;; Make fits to the tilts 
     red_pinh_make_fits, simx, simy, sx, sy $
                         , XTILTS = xtilts $ ; Input
                         , YTILTS = ytilts $ ; Input
                         , dxoffs = dxoffs $ ; Output
                         , dyoffs = dyoffs   ; Output
     
     ;; Updates
     xoffs += dxoffs
     yoffs += dyoffs

     ;; Convergence
     for ich = 1, Nch-1 do begin
        xconvergence[ich, it] = max(dxoffs[*, *, ich])-min(dxoffs[*, *, ich])
        yconvergence[ich, it] = max(dyoffs[*, *, ich])-min(dyoffs[*, *, ich])
     endfor
     mconvergence[*, *, it] = metrics

     red_pinh_plot_convergence, IT = it $ ; # of iterations so far
                                , DTOL = dtol $
                                , TTOL = ttol $
                                , XCONVERGENCE = xconvergence $    ; X tilt convergence (out)
                                , YCONVERGENCE = yconvergence $    ; Y tilt convergence (out)
                                , MCONVERGENCE = mconvergence      ; Metric convergence (out)     
     
     ;; If converged, then get out of the it loop.
     if max(abs(xconvergence[*, it])) lt ttol and $
        max(abs(yconvergence[*, it])) lt ttol then break
     if it gt 0 then begin
        dm = max(mconvergence[*, *, it]-mconvergence[*, *, it-1])
        sm = max(mconvergence[*, *, it]+mconvergence[*, *, it-1])
        if dm/sm lt mtol then break
     endif

  endfor                        ; for it
  

  ;; Do we want to refine the diversity?

  if keyword_set(finddiversity) then begin

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
           red_pinh_run_momfbd, xoffs, yoffs, simx, simy, sz $
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
        aaa = findmax2qi(-c, /verbose) ; "subpixel" minimum
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
        red_pinh_run_momfbd, xoffs, yoffs, simx, simy, sz $
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
                            , dyoffs = dyoffs  ; Output
        

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
                                      , MCONVERGENCE = mconvergence ; Metric convergence (out)     
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

     endfor                     ; for it

  endif                         ; finddiversity

  ;; Kill manager and slaves
  momfbd_setup, PORT = port, NMANAGERS=0
  ;;  spawn,'rm -f '+workdir

  ;; Crop convergence arrays.
  xconvergence = xconvergence[*, 0:it-1]
  yconvergence = yconvergence[*, 0:it-1]
  if keyword_set(finddiversity) then dconvergence = dconvergence[*, 0:it-1]
     

  ;; Write the updated offsets to files
  for ich = 1, Nch-1 do begin
     print, 'Writing ' + pinhdir + xnames[ich]
     fzwrite, fix(round(xoffs[*, *, ich])), calibdir + xnames[ich], ''
     print, 'Writing ' + pinhdir + ynames[ich]
     fzwrite, fix(round(yoffs[*, *, ich])), calibdir + ynames[ich], ''
  endfor

  fzwrite, diversity, calibdir+'diversity.fz', ''

end
