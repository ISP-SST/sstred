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
;    edgemargin : in, optional   
;
;        Number of rows/columns to exclude when searching for pinholes.
;
;    subfieldpadding : in, optional   
;
;        Pad this many pixels around the subfields when submitting to momfbd.
;        (which doesn't seem to handle exact image-sizes too well)
;
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
;   2014-01-24 : MGL. Removed diversity part of the code. 
;
;   2014-04-08 : THI. Loop over all prefilters/states unless specified.
;                New keyword show_plots, change default to not display
;                plots. 
;
;   2015-05-07 : MGL. Add 86 to list of bad subfield sizes for momfbd.
;                Also call red_pinh_run_momfbd with margin=2.
;
;   2015-05-07 : THI. Added keywords edgemargin and subfieldpadding.
;
;   2016-02-06 : THI. Bugfix: subimage cutouts could reach outside image border.
;
;   2016-08-23 : THI. Rename camtag to detector and channel to camera,
;                so the names match those of the corresponding SolarNet
;                keywords. 
;
;   2017-03-07 : MGL. Renamed this version to pinholecalib_old.
; 
;-
pro red::pinholecalib_old, STATE = state $                   
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
                       , edgemargin = edgemargin $
                       , subfieldpadding = subfieldpadding $
                       , show_plots = show_plots

    ;; Name of this method
    inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

    ;; Logging
    help, /obj, self, output = selfinfo 
    red_writelog, selfinfo = selfinfo, logfile = logfile

    if n_elements(mtol) eq 0 then mtol = 1e-6 ; Metric tolerance
    if n_elements(dtol) eq 0 then dtol = 1e-6 ; Diversity tolerance
    ;  if n_elements(ttol) eq 0 then ttol = 1e-3 ; Tilt tolerance
    if n_elements(ttol) eq 0 then ttol = 1 ; Tilt tolerance (1/100 of a pixel)

    if n_elements(nslaves) eq 0 then nslaves=6
    if n_elements(nthreads) eq 0 then nthreads=1
    if n_elements(maxit) eq 0 then maxit=400
    if n_elements(edgemargin) eq 0 then edgemargin=0
    if n_elements(subfieldpadding) eq 0 then subfieldpadding=2
    
    edgemargin += subfieldpadding      ; prevent subimage cut-outs from going outside image boundary

    ;; Where the summed pinhole images are:
    pinhdir = self.out_dir + '/pinh_align/'

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
    self -> getdetectors
    camt = self.camttag
    camr = self.camrtag
    camw = self.camwbtag
    cams = [camw, camt, camr] 

    xsz = self.camsz              ; CCD size in pixels
    ysz = self.camsz              ; CCD size in pixels

    ;; Get list of states (and transmitted file names)  
    cam = camt                    ; Transmitted camera
    files = file_search(pinhdir+cam+'.*.fpinh', count = nf)
    red_extractstates, files, /basename, pref = prefilters
    prefilters = prefilters[uniq(prefilters, sort(prefilters))]

    if n_elements(prefilter) ne 0 then begin
        indx = where(prefilters eq prefilter)
        if max(indx) eq -1 then begin
            print, inam+' : WARNING : Keyword prefilter does not match any pinhole file names: ', prefilter
            return
        endif
        prefilters = prefilters[indx]
    endif


    Npref = n_elements(prefilters)
    for ipref = 0, Npref-1 do begin

        files = file_search(pinhdir+cam+'.'+prefilters[ipref]+'*.fpinh', count = nf)

        ;; We need the total state information of the files we found.
        red_extractstates, files, /basename, fullstate = fullstates

        ;; We want to process all existing pinhole images. Unless the
        ;; keyword state is provided, then we want to process only the ones
        ;; matching this keyword.
        if n_elements(state) ne 0 then begin
            ;; Check that the state given in keywords exists. Quit with error
            ;; if it does not.
            indx = where(strmatch(fullstates,state))
            if total(indx) eq 0 then begin
                print, inam+' : state not found: '+state
                stop
            endif else begin
                fullstates = fullstates[indx]
            endelse
        endif 

        nStates = n_elements(fullstates)
        if (nStates eq 0) then continue
        
        for istat = 0, nStates-1 do begin

            ;; Set up manager and slaves.
            red_momfbd_setup, PORT = port, NSLAVES = nslaves, NMANAGERS=1, /FREE, NTHREADS=nthreads

            stat = fullstates[istat]
            red_extractstates, stat, lam = lambda

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
            clipfile = self.out_dir + '/calib/align_clips.'+prefilters[ipref]+'.sav'
            if(~file_test(clipfile)) then begin
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
            red_findpinholegrid_new, images[*, *, 0], simx, simy, dx = dx, dy = dy, Npinh = Npinh, edgemargin = edgemargin
            ;; Set subimage size to mean distance between spots.
            sz = round((dx+dy)/2.)
            ;; Make sure sz is an even number:
            sz += (sz mod 2)

            ;; Avoid sizes that make the momfbd slaves crash.
            badsizes = [86, 90]
            while total(sz eq badsizes) gt 0 do sz -= 2

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
                    ;print, 'pos', [maxloc mod sz, maxloc/sz]
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

            ;; Make fits to the "tilts" in units of 1/100 pixel.
            red_pinh_make_fits, simx, simy, sx, sy $
                             , xtilts = xtilts $ ; Input
                             , ytilts = ytilts $ ; Input
                             , xoffs = xoffs $  ; Output
                             , yoffs = yoffs    ; Output

            ;; Write the initial offsets to files
            for ich = 1, Nch-1 do begin
                print, 'Writing ' + calibdir + xnames[ich]+'.init'
                fzwrite, fix(round(xoffs[*, *, ich])), calibdir + xnames[ich]+'.init', ' '
                print, 'Writing ' + calibdir + ynames[ich]+'.init'
                fzwrite, fix(round(yoffs[*, *, ich])), calibdir + ynames[ich]+'.init', ' '
            endfor

            ;; Start iterations
            xconvergence = fltarr(Nch, MaxIt)
            yconvergence = fltarr(Nch, MaxIt)
            mconvergence = fltarr(Npinh, MaxIt)
            if keyword_set(finddiversity) then begin
                dconvergence = fltarr(Nch, MaxIt)
                dvalues = fltarr(Nch, MaxIt)
            endif

            ;; Iterate tilts only
            for it = 0, MaxIt-1 do begin

                ;; Run momfbd on the pinholes and read back the results.
                ;; Momfbd is run with the CALIBRATE keyword, so the output
                ;; tilts are in pixels. No need to convert from radians.
                red_pinh_run_momfbd, images, xoffs, yoffs, simx, simy, sz $
                                     , xtemplates, ytemplates, ptemplates $
                                     , stat $
                                     , self.telescope_d $
                                     , self.image_scale $
                                     , self.pixel_size $
                                     , DIVERSITY = diversity $
                                     , FINDDIVERSITY = 0 $
                                     , WORKDIR = workdir $
                                     , PORT = port $
                                     , NSLAVES = nslaves $
                                     , NTHREADS = nthreads $
                                     , XTILTS = dxtilts $  ; Output
                                     , YTILTS = dytilts $  ; Output
                                     , METRICS = metrics $ ; Output
                                     , subfieldpadding = subfieldpadding $
                                     , show_plots=keyword_set(show_plots)

                ;; Updates, note 100ths of a pixel
                xtilts += dxtilts*100
                ytilts += dytilts*100

                oldxoffs = xoffs
                oldyoffs = yoffs

                ;; Make fits to the tilts 
                red_pinh_make_fits, simx, simy, sx, sy $
                                    , XTILTS = xtilts $ ; Input
                                    , YTILTS = ytilts $ ; Input
                                    , xoffs = xoffs $  ; Output
                                    , yoffs = yoffs    ; Output

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

                xconvergence[*, it] = max(max(dxoffs, dim=1), dim=1) - min(min(dxoffs, dim=1), dim=1)
                yconvergence[*, it] = max(max(dyoffs, dim=1), dim=1) - min(min(dyoffs, dim=1), dim=1)
                mconvergence[*, it] = metrics

                if keyword_set(show_plots) then begin 
                    red_pinh_plot_convergence, IT = it $ ; # of iterations so far
                                           , DTOL = dtol $
                                           , TTOL = ttol $
                                           , XCONVERGENCE = xconvergence $ ; X tilt convergence (out)
                                           , YCONVERGENCE = yconvergence $ ; Y tilt convergence (out)
                                           , MCONVERGENCE = mconvergence   ; Metric convergence (out)     
                endif

                ;; If converged, then get out of the it loop.
                if max(abs(xconvergence[*, it])) lt ttol and $
                   max(abs(yconvergence[*, it])) lt ttol then begin
                   print, inam + ' : Converged due to small shifts.'
                   break
                endif
                ;        if it gt 0 then begin
                ;           dm = max(mconvergence[*, it]-mconvergence[*, it-1])
                ;           sm = max(mconvergence[*, it]+mconvergence[*, it-1])
                ;           if dm/sm lt mtol then begin
                ;              print, inam + ' : Converged due to small change in metric.'
                ;              break
                ;           endif
                ;        endif

            endfor                     ; it


            ;; Crop convergence arrays.
            xconvergence = xconvergence[*, 0:it-1]
            yconvergence = yconvergence[*, 0:it-1]
            mconvergence = mconvergence[*, 0:it-1]
            if keyword_set(finddiversity) then dconvergence = dconvergence[*, 0:it-1]

            ;; Write the updated offsets to files
            for ich = 1, Nch-1 do begin
                print, 'Writing ' + calibdir + xnames[ich]
                fzwrite, fix(round(xoffs[*, *, ich])), calibdir + xnames[ich], ' '
                print, 'Writing ' + calibdir + ynames[ich]
                fzwrite, fix(round(yoffs[*, *, ich])), calibdir + ynames[ich], ' '
            endfor         

           ;; Kill manager and slaves
           red_momfbd_setup, PORT = port, NMANAGERS=0

        endfor                       ; istat
        
    endfor                       ; ipref

    ;;spawn,'rm -f '+workdir

    if keyword_set(finddiversity) then begin
        print, inam+' : Do you want to refine the diversity?'
        print, inam+' : '
        print, inam+' : There is some code in red__pinholecalib_diversity.pro, that  '
        print, inam+' : could be used as inspiration for what should go here. That   '
        print, inam+' : code is not well tested because we currently do not have a   '
        print, inam+' : CRISP diversity camera. It was written for the purpose of    '
        print, inam+' : finding small diversity differences between nominally focused'
        print, inam+' : cameras. We never were able to make this work in a way we    '
        print, inam+' : trusted.                                                     '
        print, inam+' : '
        print, inam+' : Be prepared to do some debugging work or possibly rewriting  '
        print, inam+' : of that code. /MGL '
    endif                      ; finddiversity


end
