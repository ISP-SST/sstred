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
;    sz :                           
;
;
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
; 
; 
; :History:
; 
;   2013-09-04 : MGL. Use red_momfbd_*, not momfbd_*.
; 
; 
; 
;-
pro red_pinh_run_momfbd, xoffs, yoffs, simx, simy, sz $
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
                         , METRICS = metrics 

  ana_output = 1                ; Otherwise use MOMFBD format output (for development/debugging).

  Nch = (size(xoffs, /dim))[2] 
  Npinhx = (size(simx, /dim))[0]
  Npinhy = (size(simy, /dim))[0]

  if n_elements(diversity) eq 0 then diversity = replicate(0.0, Nch)
  if n_elements(finddiversity) eq 0 then finddiversity = 0
  if n_elements(workdir) eq 0 then workdir = './'
  if n_elements(nslaves) eq 0 then nslaves = 6
  if n_elements(nthreads) eq 0 then nthreads = 1

  ;; Restart manager and slaves if needed.
  red_momfbd_setup, PORT = port, NSLAVES = nslaves, NMANAGERS=1, NTHREADS=nthreads
  
  ;; Remove old metric files if needed.
  spawn, 'cd '+workdir+' ; rm -f metric.*'

  
  ;; New subfield offsets, submit jobs
  imno = 0
  d=sz/2
  for ix=0,Npinhx-1 do begin
     for iy=0,Npinhy-1 do begin

        ;; Read out subfield (ix,iy) from offsets
        xosubfields = xoffs[simx[ix]-d:simx[ix]+d-1,simy[iy]-d:simy[iy]+d-1, *]
        yosubfields = yoffs[simx[ix]-d:simx[ix]+d-1,simy[iy]-d:simy[iy]+d-1, *]

        for ich = 0, Nch-1 do begin
           if ich gt 0 then begin
              fname = workdir+string(format = '(%"'+xtemplates[ich]+'")', imno)
              fzwrite, fix(round(xosubfields[*, *, ich])), fname, ''
              fname = workdir+string(format = '(%"'+ytemplates[ich]+'")', imno)
              fzwrite, fix(round(yosubfields[*, *, ich])), fname, ''
           endif
        endfor

        ;; Construct a config file for pinhole calibration. 
        cfglines = ['object{']
        cfglines = [cfglines, '  WAVELENGTH='+strtrim(string(long(stat)*1e-10),2)] 
        for ich = 0, Nch-1 do begin
           cfglines = [cfglines, '  channel{']
           cfglines = [cfglines, '    IMAGE_DATA_DIR=' + workdir]
           cfglines = [cfglines, '    FILENAME_TEMPLATE=' + ptemplates[ich]]
           cfglines = [cfglines, '    DIVERSITY='+strtrim(string(diversity[ich]), 2)+' mm']
;           cfglines = [cfglines, '    '+acl[ich]]
           if ich gt 0 then begin
              cfglines = [cfglines, '    XOFFSET='+ workdir $
                          + string(format = '(%"'+xtemplates[ich]+'")', imno)]
              cfglines = [cfglines, '    YOFFSET='+ workdir $
                          + string(format = '(%"'+ytemplates[ich]+'")', imno)]
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
        cfglines = [cfglines, 'IMAGE_NUMS='+string(imno, format = '(i0)')]
        cfglines = [cfglines, 'BORDER_CLIP=0']
        cfglines = [cfglines, 'GET_METRIC']

        ;; Submit to momfbd
        red_momfbd_submit, imno, PORT=port, CFGLINES=cfglines, DIR=workdir, /FORCE
        
        imno += 1
        
     endfor                     ; for iy
  endfor                        ; for ix
  

  ;; Loop through pinholes and read results
  xtilts = fltarr(Nch, Npinhx, Npinhy)
  ytilts = fltarr(Nch, Npinhx, Npinhy)
  foc = fltarr(Nch, Npinhx, Npinhy)
  metrics = fltarr(Npinhx, Npinhy)
  imno = 0
  if ~ana_output then mrpointers = ptrarr(Npinhx*Npinhy, /allocate_heap)

  ;; Pre-read the metric files. They are numbered with the momfbd
  ;; job numbers, not imno.
  repeat begin
     red_momfbd_check, PORT = port, NJOBS = njobs
     print, 'njobs=', njobs
     wait, 3
  endrep until njobs eq 0
  spawn, 'cd '+workdir+' ; cat metric.*', metriclist

  for ix=0, Npinhx-1 do begin
     for iy=0, Npinhy-1 do begin

        metrics[ix, iy] = float(metriclist[imno])
        
        if arg_present(xtilts) or arg_present(ytilts) then begin

           ;; Read momfbd output
           print, 'Read results for ', ix, iy

           if ana_output then begin
              for ich = 0, Nch-1 do begin

                 fname = 'momfbd_submit.'+strtrim(string(imno), 2)+'.alpha.'+strtrim(string(ich+1), 2)+'.1'

                 openr, flun, /get_lun, workdir+fname
                 oneline = ''
                 readf, flun, oneline
                 readf, flun, oneline
                 free_lun, flun

                 tmp = strsplit(oneline, ' ', /extract)
                 xtilts[ich, ix, iy] = float(tmp[2])
                 ytilts[ich, ix, iy] = float(tmp[3])
                 if keyword_set(finddiversity) then foc[ich, ix, iy] = float(tmp[4])

              endfor
           endif else begin     ; ana_output
              *mrpointers[imno] = red_momfbd_results(imno, PORT=port, DIR=workdir, /NOWAIT)
              
              xtilts[*, ix, iy] = (*mrpointers[imno]).patch.alpha[0, *]
              ytilts[*, ix, iy] = (*mrpointers[imno]).patch.alpha[1, *]
              if keyword_set(finddiversity) then foc[*, ix, iy] = (*mrpointers[imno]).patch.alpha[2, *]

           endelse              ; momfbd_output

        endif                   ; arg_present

        imno += 1

     endfor                     ; for iy
  endfor                        ; for ix

  ;; ptr_free, mrpointers
  
end

