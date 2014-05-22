; docformat = 'rst'

;+
; Run C++ momfbd on pinholes.
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
;
; :Params:
;
;   images : 
;
;      The full-FOV pinhole images.
;
;    xoffs : 
;
;
;
;    yoffs : 
;
;
;
;    simx : 
;
;
;
;    simy : 
;
;
;
;    sz : in, type=integer                     
;
;       Subfield size used for momfbding.
;
;    xtemplates : 
;
;
;
;    ytemplates : 
;
;
;
;    ptemplates :                           
;
;
;
;    stat :                           
;
;
;
;    telescope_d : 
;
;
;
;    image_scale : 
;
;
;
;    pixel_size : 
;
;
;
;
; :Keywords:
; 
;    diversity : in, optional
; 
; 
; 
;    finddiversity : in, optional
; 
; 
; 
;    workdir : in, optional
; 
; 
; 
;    port : in, optional
; 
; 
; 
;    nslaves : in, optional
; 
; 
; 
;    nthreads : in, optional
; 
; 
; 
;    xtilts : out, optional 
; 
; 
; 
;    ytilts : out, optional 
; 
; 
; 
;    foc : out, optional 
; 
; 
; 
;    metrics : out, optional 
;
;    margin : in, optional, type=integer
; 
;       Margin for momfbd to use for swapping tilts for image shifts.
; 
; 
; 
; :History:
; 
;   2013-09-04 : MGL. Use red_momfbd_*, not momfbd_*.
;
;   2013-09-10 : MGL. Moved writing of subfield images to disk (and
;                the associated dark level refinement) from
;                red::pinholecalib. Add a margin to the subfields
;                written to disk.  
;
;   2014-01-27 : MGL. New keyword: margin. Adapt to input simx, simy,
;                etc., being 1D arrays.
; 
; 
;   2014-04-08 : THI. New keyword show_plots, change default to not
;                display plots
; 
;-
pro red_pinh_run_momfbd, images, xoffs, yoffs, simx, simy, sz $
                         , xtemplates, ytemplates, ptemplates $
                         , stat $
                         , telescope_d, image_scale, pixel_size $
                         , DIVERSITY = diversity $
                         , FINDDIVERSITY = finddiversity $
                         , WORKDIR = workdir $
                         , PORT = port $
                         , NSLAVES = nslaves $
                         , NTHREADS = nthreads $
                         , XTILTS = xtilts $ 
                         , YTILTS = ytilts $ 
                         , FOC = foc $       
                         , METRICS = metrics  $ 
                         , margin = margin $
                         , show_plots = show_plots

  ana_output = 1                ; Otherwise use MOMFBD format output (for development/debugging).

  Nch = (size(xoffs, /dim))[2] 
  Npinh = (size(simx, /dim))[0]

  if n_elements(margin) eq 0 then margin = 0
  if n_elements(diversity) eq 0 then diversity = replicate(0.0, Nch)
  if n_elements(finddiversity) eq 0 then finddiversity = 0
  if n_elements(workdir) eq 0 then workdir = './'
  if n_elements(nslaves) eq 0 then nslaves = 6
  if n_elements(nthreads) eq 0 then nthreads = 1

  ;; Restart manager and slaves if needed.
  red_momfbd_setup, PORT = port, NSLAVES = nslaves, NMANAGERS=1, NTHREADS=nthreads
  
  ;; Remove old metric files if needed.
  spawn, 'cd '+workdir+' ; rm -f metric.*'

  if keyword_set(show_plots) then begin 
      window, 0
  endif
  
  ;; New subfield offsets, submit jobs
  d = sz/2 + margin

  for ihole=0,Npinh-1 do begin
     
     ;; Read out subfield from images and offsets
     xosubfields = xoffs[simx[ihole]-d:simx[ihole]+d-1,simy[ihole]-d:simy[ihole]+d-1, *]
     yosubfields = yoffs[simx[ihole]-d:simx[ihole]+d-1,simy[ihole]-d:simy[ihole]+d-1, *]
     pinhsubfields = images[simx[ihole]-d:simx[ihole]+d-1,simy[ihole]-d:simy[ihole]+d-1, *]
     
     Nbins = 1000
     for ich = 0, Nch-1 do begin
        
        pinhsubfields[*, *, ich] = pinhsubfields[*, *, ich]/max(pinhsubfields[*, *, ich])
        
        ;; Fine tune dark level and normalization for each subfield.
        ;; Get the bias by fitting a Gaussian to the histogram peak.
        
        ;; Calculate some statistics on the background
        sindx = where(pinhsubfields[*,*,ich] lt 0.01*max(pinhsubfields[*,*,ich]))
        mn = median((pinhsubfields[*,*,ich])[sindx])
        st = stdev((pinhsubfields[*,*,ich])[sindx])
        
        ;;mn = biweight_mean(pinhsubfields[*,*,ich])
        ;;st = robust_sigma(pinhsubfields[*,*,ich])
        
        hmax = mn + 30.*st
        hmin = mn - 30.*st
        
        hh = histogram(pinhsubfields[*, *, ich], min = hmin, max = hmax $
                       , Nbins = Nbins, locations = locations)
        binsize = (hmax - hmin) / (Nbins - 1)
        intensities = locations + binsize/2.
        
        dummy = max(smooth(hh,15),ml)
        
        Nfit = 51
        indx = ml + indgen(Nfit) - (Nfit-1)/2
        
        yfit=mpfitpeak(float(intensities[indx]), float(hh[indx]), a, Nterms = 5)
        if keyword_set(show_plots) then begin 
            plot,intensities[indx],hh[indx], /xstyle, /ystyle
            oplot,intensities[indx],yfit,color=fsc_color('yellow')
            oplot,[0,0]+a[1],[0,1]*max(hh),color=fsc_color('yellow')
            oplot,[0,0]+mn,[0,1]*max(hh),color=fsc_color('cyan')
        endif
        
        ;print,intensities[ml]
        
        ;; Subtract the fitted dark level and renormalize to unit max.
        ;; The momfbd program will then renormalize all channels to
        ;; the same total energy.
        pinhsubfields[*, *, ich] = pinhsubfields[*, *, ich] - a[1]
        pinhsubfields[*, *, ich] = pinhsubfields[*, *, ich]/max(pinhsubfields[*, *, ich])
        
        fname = workdir+string(format = '(%"'+ptemplates[ich]+'")', ihole)
        fzwrite, fix(round(pinhsubfields[*, *, ich]*15000.)), fname, ' '

        if ich gt 0 then begin
           fname = workdir+string(format = '(%"'+xtemplates[ich]+'")', ihole)
           fzwrite, fix(round(xosubfields[*, *, ich])), fname, ' '
           fname = workdir+string(format = '(%"'+ytemplates[ich]+'")', ihole)
           fzwrite, fix(round(yosubfields[*, *, ich])), fname, ' '
        endif

     endfor                     ; ich

     ;; Construct a config file for pinhole calibration. 
     cfglines = ['object{']
     cfglines = [cfglines, '  WAVELENGTH='+strtrim(string(long(stat)*1e-10),2)] 
     for ich = 0, Nch-1 do begin
        cfglines = [cfglines, '  channel{']
        cfglines = [cfglines, '    IMAGE_DATA_DIR=' + workdir]
        cfglines = [cfglines, '    FILENAME_TEMPLATE=' + ptemplates[ich]]
        cfglines = [cfglines, '    DIVERSITY='+strtrim(string(diversity[ich]), 2)+' mm']
        if ich gt 0 then begin
           cfglines = [cfglines, '    XOFFSET='+ workdir $
                       + string(format = '(%"'+xtemplates[ich]+'")', ihole)]
           cfglines = [cfglines, '    YOFFSET='+ workdir $
                       + string(format = '(%"'+ytemplates[ich]+'")', ihole)]
        endif
        cfglines = [cfglines, '  }']
     endfor
     cfglines = [cfglines, '}']
     cfglines = [cfglines, 'BASIS=Zernike']
     cfglines = [cfglines, 'GRADIENT=gradient_diff']
     cfglines = [cfglines, 'GETSTEP=getstep_steepest_descent']
     cfglines = [cfglines, 'PROG_DATA_DIR=./data/']
     if keyword_set(finddiversity) then begin
        cfglines = [cfglines, 'MODES=2-4']                               ;; Tilts + focus
     endif else begin
        cfglines = [cfglines, 'MODES=2-3']                               ;; Tilts 
     endelse
     cfglines = [cfglines, 'NUM_POINTS='+strtrim(string(round(sz)),2)]
     cfglines = [cfglines, 'TELESCOPE_D='+strtrim(string(telescope_d),2)]
     cfglines = [cfglines, 'ARCSECPERPIX='+strtrim(string(image_scale),2)]
     cfglines = [cfglines, 'PIXELSIZE='+strtrim(string(pixel_size),2)]
     if ana_output then begin
        cfglines = [cfglines, 'FILE_TYPE=ANA']
     endif else begin
        cfglines = [cfglines, 'FILE_TYPE=MOMFBD']
     endelse
     cfglines = [cfglines, 'CALIBRATE']
     ;; Want to process this pinhole only once, centered.
     cfglines = [cfglines, 'SIM_X='+red_stri(d)]
     cfglines = [cfglines, 'SIM_Y='+red_stri(d)]
;        cfglines = [cfglines, 'IMAGE_NUMS='+string(ihole, format = '(i0)')]
     cfglines = [cfglines, 'BORDER_CLIP=0']
     cfglines = [cfglines, 'GET_METRIC']
     
     ;; Submit to momfbd
     red_momfbd_submit, ihole, PORT=port, CFGLINES=cfglines, DIR=workdir, /FORCE

  endfor                        ; for ihole
  

  ;; Loop through pinholes and read results
  xtilts = fltarr(Nch, Npinh)
  ytilts = fltarr(Nch, Npinh)
  foc = fltarr(Nch, Npinh)
  metrics = fltarr(Npinh)
  imno = 0
  if ~ana_output then mrpointers = ptrarr(Npinh, /allocate_heap)

  ;; Pre-read the metric files. They are numbered with the momfbd
  ;; job numbers, not imno.
  repeat begin
     red_momfbd_check, PORT = port, NJOBS = njobs
     print, 'njobs=', njobs
     wait, 3
  endrep until njobs eq 0
  spawn, 'cd '+workdir+' ; cat metric.*', metriclist


  for ihole=0,Npinh-1 do begin
     
     metrics[ihole] = float(metriclist[imno])
        
     if arg_present(xtilts) or arg_present(ytilts) then begin

        ;; Read momfbd output
        ;print, 'Read results for ', ihole

        if ana_output then begin
           for ich = 0, Nch-1 do begin

              fname = 'momfbd_submit.'+strtrim(string(ihole), 2)+'.alpha.'+strtrim(string(ich+1), 2)+'.1'
              
              openr, flun, /get_lun, workdir+fname
              oneline = ''
              readf, flun, oneline
              readf, flun, oneline
              free_lun, flun
              
              tmp = strsplit(oneline, ' ', /extract)
              xtilts[ich, ihole] = float(tmp[2])
              ytilts[ich, ihole] = float(tmp[3])
              if keyword_set(finddiversity) then foc[ich, ihole] = float(tmp[4])
              
           endfor
        endif else begin        ; end ana_output, begin momfbd_output
           *mrpointers[ihole] = red_momfbd_results(ihole, PORT=port, DIR=workdir, /NOWAIT)
           
           xtilts[*, ihole] = (*mrpointers[ihole]).patch.alpha[0, *]
           ytilts[*, ihole] = (*mrpointers[ihole]).patch.alpha[1, *]
           if keyword_set(finddiversity) then foc[*, ihole] = (*mrpointers[ihole]).patch.alpha[2, *]
           
        endelse                 ; end momfbd_output
        
     endif                      ; arg_present
     
  endfor                        ; for ihole

  ;;  if ~ana_output then ptr_free, mrpointers
  
end

