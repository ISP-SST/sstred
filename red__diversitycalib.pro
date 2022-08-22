; docformat = 'rst'

;+
;   Find the amount of defocusing of the diversity camera. 
;
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, ISP, 2017
;
; 
; 
; :Keywords:
;    
;    threshold : in, optional, type=float, default=0.25
;
;       Threshold for identifying a strong enough pinhole.
;
;    nref : in, optional, type=integer, default=5
;
;      How many of the strongest pinholes to use for finding the
;      approximate transform. Afterwards a refinement is made using
;      >80% of the detected pinholes.
;
;    pref : in, optional, type=string
;
;      Indicate the prefilter you want to calculate the clips for,
;      Default is to do it for all prefilters there is data for.
;
; :History:
;
;      2017-01-25 : MGL. First version.
;
;-
pro red::diversitycalib, threshold = threshold $
                         , nref = nref $
                         , pref = pref $
                         , dir = dir $
                         , nslaves = nslaves $           
                         , nthreads = nthreads $         
                         , verbose = verbose

  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)                                

  ;; Logging
  help, /obj, self, output = selfinfo 
  red_writelog, selfinfo = selfinfo, logfile = logfile

  show_plots = 1                ; Should be keyword!

  if n_elements(nslaves) eq 0 then nslaves=6
  if n_elements(nthreads) eq 0 then nthreads=1

  if(n_elements(threshold) eq 0) then threshold = 0.25
  if(n_elements(nref) eq 0) then nref = 5
  if( n_elements(dir) gt 0 ) then dir = [dir] $
  else if ptr_valid(self.pinh_dirs) then dir = *self.pinh_dirs
  self -> getdetectors, dir = dir
  
  if(n_elements(verbose) eq 0) THEN verbose = 0

  if n_elements(mtol) eq 0 then mtol = 1e-6 ; Metric tolerance
  if n_elements(dtol) eq 0 then dtol = 1e-6 ; Diversity tolerance
  ;;  if n_elements(ttol) eq 0 then ttol = 1e-3 ; Tilt tolerance
  if n_elements(ttol) eq 0 then ttol = 5 ; Tilt tolerance (1/100 of a pixel)


  ;; Find the diversity camera
  pdindx = where(strmatch(*self.cameras,'*-D'), Npd)
  case Npd of
    0: begin
      print, 'No pd camera: '
      print, *self.cameras
      stop
    end
    1: begin
      pd_cam = (*self.cameras)[pdindx[0]]
      pd_detector = (*self.detectors)[pdindx[0]]
    end
    else: begin
      print, 'More than one -D camera:'
      print, *self.cameras
      stop
    end
  endcase


  ;; Find the reference camera
  refindx = where(strmatch(*self.cameras,'*-W'), Nref)
  case Nref of
    0: begin
      print, 'No wideband camera: '
      print, *self.cameras
      stop
    end
    1: begin
      ref_cam = (*self.cameras)[refindx[0]]
      ref_detector = (*self.detectors)[refindx[0]]
    end
    else: begin
      print, 'More than one -W camera:'
      print, *self.cameras
      stop
    end
  endcase

  ref_caminfo = red_camerainfo(ref_detector)

  ph_dir = self.out_dir+'/pinhs/'
  ph_dir = red_strreplace(ph_dir,'//','/')

  output_dir = self.out_dir+'/diversity/'
  file_mkdir, output_dir

  workdir = output_dir+'work/'
  file_mkdir, workdir
  file_mkdir, workdir+'data'

  alignmentsfile = self.out_dir+'/calib/alignments.sav'
  if file_test(alignmentsfile) eq 0 then begin
    print, inam + ' : ' + alignmentsfile + ' not found'
    print, '  Did you run a -> pinholecalib_thi?'
  endif
  restore, alignmentsfile

  pfiles = file_search( ph_dir + '*.pinh' )
  nf = n_elements(pfiles)
  
  self -> selectfiles, files=pfiles, states=states, cam = ref_cam, selected=refselection, count = Nref
  self -> selectfiles, files=pfiles, states=states, cam = pd_cam,  selected=pdselection,  count = Npd

  if Nref ne Npd then begin
    print, inam+' : Number of ref files and pd files do not match.'
    stop
  endif
  
  refstates = states[refselection]
  pdstates = states[pdselection]
  
  reffiles = pfiles[refselection]
  pdfiles = pfiles[pdselection]
  
  diversity_measured_all = fltarr(Nref)

  
 
;  for iref = 0, Nref-1 do begin

  for iref = 0, 0 do begin

    ipd = where(refstates[iref].fpi_state eq pdstates.fpi_state, Nmatch)
    if Nmatch eq 0 then continue
    
    print, reffiles[iref]
    print, pdfiles[ipd]

    im_ref = red_readdata(reffiles[iref])
    im_pd  = red_readdata(pdfiles[ipd])

    ;; Dark level refinement?

    if 0 then begin
      ;; Remap the pd image
      ialign = where(file_basename(reffiles[iref]) eq file_basename(alignments.state1.filename) $
                     and file_basename(pdfiles[iref]) eq file_basename(alignments.state2.filename) $
                     , Nmatch)

      if Nmatch ne 1 then stop

      map = alignments[ialign[0]].map
      im_pd_remapped = rdx_img_project(map, im_pd) ; Note: not perfect alignment for ref/pd!
    endif else begin
      ;; Clip both images
      clip_ref = [62,1859,6,1198]
      clip_pd  = [1849,52,4,1196]
      im_ref         = red_clipim(im_ref, clip_ref)
      im_pd_remapped = red_clipim(im_pd, clip_pd)
    endelse
    dims = size(im_ref, /dim)

    ;; Find pinholes
    red_findpinholegrid_new, im_ref, px, py , dx = dx, dy = dy $
                             , Npinh = Npinh, thres = threshold
    
;    ;; Construct a mask for dark level refinement fitting
;    xx = findgen(dims) mod dims[0]
;    yy = transpose(findgen(reverse(dims)) mod dims[1])
;
;    pxm = px+dx/2
;    pym = py+dy/2
;    pindx = where(pxm lt max(px) and pym lt max(py))
;    pxm = pxm[pindx]
;    pym = pym[pindx]
;    mask = fltarr(dims)
;    mrd = 15                     ; mask radius
;    for ipinh = 0, n_elements(pindx)-1 do mask += (sqrt((xx-pxm[ipinh])^2+(yy-pym[ipinh])^2) le mrd)
;    mindx = where(mask)
;
;    ;; Dark correction for the reference image
;    ccc = planefit(xx[mindx], yy[mindx], im_ref[mindx], 0.)
;    darkplane = ccc(0)+ccc(1)*xx+ccc(2)*yy
;    im_ref -= darkplane
;
;    ;; Dark correction for the diversity image
;    ccc_pd = planefit(xx[mindx], yy[mindx], im_pd_remapped[mindx], 0.)
;    darkplane_pd = ccc_pd(0)+ccc_pd(1)*xx+ccc_pd(2)*yy
;    im_pd_remapped -= darkplane_pd

    images = [[[im_ref]],[[im_pd_remapped]]]

    ;; Set subimage size to mean distance between spots.
    sz = round((dx+dy)/2.)
    ;; Make sure sz is an even number:
    sz -= (sz mod 2)
    ;; Avoid sizes that make the momfbd slaves crash.
    badsizes = [86, 90]
    while total(sz eq badsizes) gt 0 do sz -= 2
    

    ;; Measure positions used for calculating initial offsets.
    d = sz/2
    Nch = 2                     ; Two channels, ref and pd
    pos = fltarr(2, Nch, Npinh) ; Peak positions within subfield
    for ihole=0,Npinh-1 do begin
      ;; Read out subfield (ix,iy) from images
      subfields = images[px[ihole]-d:px[ihole]+d-1,py[ihole]-d:py[ihole]+d-1, *]
      for ich = 0, Nch-1 do begin
        pos[*, ich, ihole] = centroid(subfields[*, *, ich])
;        strehl[ich, ihole] = max(subfields[*, *, ich])/mean(subfields[*, *, ich])
;        ;; Peak position
;        maxloc = (where(subfields[*,*,ich] eq max(subfields[*,*,ich])) )[0]
;        pos[*, ich, ihole] = [maxloc mod sz, maxloc/sz]
                                ;print, 'pos', [maxloc mod sz, maxloc/sz]
      endfor                    ; ich
    endfor                      ; ihole
    
    ;; Make the tilts in units of 1/100 pixel.
    xtilts = fltarr(Nch, Npinh)
    xtilts[1, *] = (pos[0, 1, *] - pos[0, 0, *])*100

    ytilts = fltarr(Nch, Npinh)
    ytilts[1, *] = (pos[1, 1, *] - pos[1, 0, *])*100


    ;; Make fits to the "tilts" in units of 1/100 pixel.
    red_pinh_make_fits, px, py, dims[0], dims[1] $
                        , xtilts = xtilts $ ; Input
                        , ytilts = ytilts $ ; Input
                        , xoffs = xoffs $   ; Output
                        , yoffs = yoffs     ; Output


    ;; File name templates for the subfield-size images and offsets.
    stat = refstates[iref].fullstate
    cams = [ref_cam, pd_cam]
    ptemplates = cams + '.' + stat + '.pinh.%03d' 
    xtemplates = cams + '.' + stat + '.xoffs.%03d'
    ytemplates = cams + '.' + stat + '.yoffs.%03d'
    

    ;; Wavelength in Å is passed to momfbd as a string.
    stat = string(states[iref].PF_WAVELENGTH*1e10,format='(i04)') ; or just states[iref].prefilter

    diversity_measured_thispair = fltarr(Npinh)

    
    ;; Set up manager and slaves.
    red_momfbd_setup, PORT = port, NSLAVES = nslaves, NMANAGERS=1, /FREE, NTHREADS=nthreads

;    diversity_values = self.diversity * [ 1 + [-10, -5, -3, -1, 0, 1, 3, 5, 10, 15]*0.01 ]
;    Ndiversities = n_elements(diversity_values)
;    metric_values = fltarr(Ndiversities, Npinh)

    ;; Start iterations
    maxit = 100
    xconvergence = fltarr(Nch,   MaxIt)
    yconvergence = fltarr(Nch,   MaxIt)
    mconvergence = fltarr(Npinh, MaxIt)
    dconvergence = fltarr(Nch,   MaxIt)
    dvalues      = fltarr(Nch,   MaxIt)

    diversity = float(self.diversity) ; Initial value
;    diversity *= 4.
    
    ;; First iterate tilts only, then also with diversity fitting
    
    for finddiversity = 0, 1 do begin

      for it = 0, MaxIt-1 do begin
        ;;for idiv = 0, Ndiversities-1 do begin

        ;; Run momfbd on the pinholes and read back the results.
        ;; Momfbd is run with the CALIBRATE keyword, so the output
        ;; tilts are in pixels. No need to convert from radians.
        red_pinh_run_momfbd, images $
                             , xoffs = xoffs, yoffs = yoffs $
                             , px, py, sz $
                             , xtemplates = xtemplates, ytemplates = ytemplates $
                             , ptemplates $
                             , stat $
                             , self.telescope_d $
                             , self.image_scale $
                             , strtrim(ref_caminfo.pixelsize, 2) $
                             , DIVERSITY = [0, diversity*1e3] $
;                           , DIVERSITY = [0, diversity_values[idiv]*1e3] $
                             , finddiversity = finddiversity $
                             , foc = foc $ ; Output focus in mm
                             , WORKDIR = workdir $
                             , PORT = port $
                             , NSLAVES = nslaves $
                             , NTHREADS = nthreads $
                             , XTILTS = dxtilts $  ; Output in pixels/100
                             , YTILTS = dytilts $  ; Output in pixels/100
                             , METRICS = metrics $ ; Output
                             , subfieldpadding = subfieldpadding $
                             , show_plots=0


        ;; Make weights based on outliers
        mx = biweight_mean(dytilts[1,*],ms,wy)
        my = biweight_mean(dxtilts[1,*],ms,wx)

        ww = wx/max(wx) * wy/max(wy)

        if finddiversity then begin

          dfoc = (foc[1, *] - foc[0, *])/2.
          mf = biweight_mean(dfoc,ms,wf)
          ww *= wf/max(wf)

          cgplot,(dfoc*wf/max(wf) )[where(wf ne 0)]
;stop          
          dfoc = total(dfoc*ww)/total(ww)
          olddiversity = diversity
          diversity += dfoc
          dvalues[1, it] = diversity
          dconvergence[1, it] = diversity-olddiversity

        endif

        ;; Updates, weighted to remove outliers
        xtilts[1, *] += ww*dxtilts[1, *]/2 ;*100
        ytilts[1, *] += ww*dytilts[1, *]/2 ;*100

        oldxoffs = xoffs
        oldyoffs = yoffs

        ;; Make fits to the tilts 
        red_pinh_make_fits, px, py, dims[0], dims[1] $
                            , XTILTS = xtilts $ ; Input
                            , YTILTS = ytilts $ ; Input
                            , xoffs = xoffs $   ; Output
                            , yoffs = yoffs     ; Output

                                ;     ;; Updates
                                ;     xoffs += dxoffs
                                ;     yoffs += dyoffs

        dxoffs = xoffs-oldxoffs
        dyoffs = yoffs-oldyoffs


        red_stats, oldxoffs, name = 'oldxoffs'
        red_stats, xoffs, name = 'xoffs'
        red_stats, dxoffs, name = 'dxoffs'

        red_stats, yoffs, name = 'yoffs'
        red_stats, dyoffs, name = 'dyoffs'

        ;; Convergence
        for ich = 1, Nch-1 do begin
          xconvergence[ich, it] = max(dxoffs[*, *, ich])-min(dxoffs[*, *, ich])
          yconvergence[ich, it] = max(dyoffs[*, *, ich])-min(dyoffs[*, *, ich])
          mconvergence[ich, it] = median(metrics[*])
        endfor
;      xconvergence[*, it] = max(max(dxoffs, dim=1), dim=1) - min(min(dxoffs, dim=1), dim=1)
;      yconvergence[*, it] = max(max(dyoffs, dim=1), dim=1) - min(min(dyoffs, dim=1), dim=1)
;        mconvergence[*, it] = metrics
        
        if keyword_set(show_plots) then begin 
          if finddiversity then begin
            red_pinh_plot_convergence, IT = it $ ; # of iterations so far
                                       , DTOL = dtol $
                                       , TTOL = ttol $
                                       , XCONVERGENCE = xconvergence $ ; X tilt convergence (out)
                                       , YCONVERGENCE = yconvergence $ ; Y tilt convergence (out)
                                       , MCONVERGENCE = mconvergence $ ; Metric convergence (out)     
                                       , dconvergence = dconvergence $ ; Diversity convergence (out)     
                                       , dvalues = dvalues
          endif else begin
            red_pinh_plot_convergence, IT = it $ ; # of iterations so far
                                       , DTOL = dtol $
                                       , TTOL = ttol $
                                       , XCONVERGENCE = xconvergence $ ; X tilt convergence (out)
                                       , YCONVERGENCE = yconvergence $ ; Y tilt convergence (out)
                                       , MCONVERGENCE = mconvergence   ; Metric convergence (out)     
          endelse
        endif


        wset, 0
        tvscl, congrid(im_ref, dims[0]/2, dims[1]/2)
;tvscl, congrid(im_pd_remapped, dims[0]/2, dims[1]/2)
        tvscl, congrid(red_applyoffsets(im_pd_remapped, xoffs[*, *, 1], yoffs[*, *, 1]), dims[0]/2, dims[1]/2)
;tvscl, congrid(red_applyoffsets(im_pd_remapped, oldxoffs[*, *, 1], oldyoffs[*, *, 1]), dims[0]/2, dims[1]/2)
        

        ;; If converged, then get out of the it loop.
        if finddiversity then begin
          
          if abs(dconvergence[1, it]) lt dtol then begin
            print, inam + ' : Converged due to small diversity update.'
            break
          endif
        endif else begin
          if max(abs(xconvergence[*, it])) lt ttol and $
             max(abs(yconvergence[*, it])) lt ttol then begin
            print, inam + ' : Converged due to small shifts.'
            break
          endif
        endelse

;      metric_values[idiv, *] = metrics

      endfor                    ; it
    endfor                      ; finddiversity

    red_pinh_plot_convergence, IT = it-1 $ ; # of iterations so far
                               , DTOL = dtol $
                               , TTOL = ttol $
                               , XCONVERGENCE = xconvergence $ ; X tilt convergence (out)
                               , YCONVERGENCE = yconvergence $ ; Y tilt convergence (out)
                               , MCONVERGENCE = mconvergence $ ; Metric convergence (out)     
                               , dconvergence = mconvergence   ; Diversity convergence (out)     


    red_momfbd_setup, PORT = port, NMANAGERS=0

;    cgplot,diversity_values*1e3,mean(metric_values,dim=2),psym=16,/yno $
;           , xtitle = 'diversity focus / 1 mm', ytitle = 'L'

    
  endfor                        ; istate

  print, string(diversity*1e3, format = '(f4.1)')+' mm'
  openw, llun, /get_lun, 'calib/diversity.txt'
  printf, llun, string(diversity*1e3, format = '(f4.1)')+' mm'
  free_lun, llun

stop

end
