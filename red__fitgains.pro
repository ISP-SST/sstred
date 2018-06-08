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
;    all : in, optional, type=boolean
;
;       If set, bypass selection and process all data.
;
;    extra_nodes : in, optional, type=array
;
;       Make extra spline nodes at the specified wavelengths, given in
;       Ångström from the core point.
;
;    pref : in, optional, type=string
;
;       Select prefilter to process.
;
;    npar  : 
;   
;   
;   
;    niter  : 
;   
;   
;   
;    rebin  : 
;   
;   
;   
;    xl  : 
;   
;   
;   
;    yl  : 
;   
;   
;   
;    densegrid  : 
;   
;   
;   
;    res  : 
;   
;   
;   
;    thres  : 
;   
;   
;   
;    initcmap  : 
;   
;   
;   
;    x0  : 
;   
;   
;   
;    x1  : 
;   
;   
;   
;    state  : 
;   
;   
;   
;    nosave  : 
;   
;   
;   
;    myg  : 
;   
;   
;   
;    w0  : 
;   
;   
;   
;    w1  : 
;   
;   
;   
; 
; 
; :History:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2013-08-27 : MGL. Added support for logging. Let the subprogram
;                find out its own name.
; 
;   2013-12-20 : PS join with _ng version: reflectivity as keyword,
;                npar is number of *additional* terms (default: 1 =
;                linear)
; 
;   2013-12-21 : MGL. Changed default npar to 2.
; 
;   2016-03-22 : JLF. Added support for .lc4 flat cubes (see 
;  	 	 red::prepflatcubes_lc4). Added a bit of error checking.
;
;   2016-10-04 : MGL. Adapted for CHROMIS data and split_classes. Set
;                up for red_get_imean to use cgplot so plots can
;                easily be saved. Individual npar defaults for
;                prefilters, so we can batch process. Get output file
;                names from filenames method.
;
;   2016-10-26 : MGL. Special case for prefilters with only a single
;                wavelength point.
;
;   2017-04-13 : MGL. Make SOLARNET FITS headers.
;
;   2017-04-20 : MGL. New keyword extra_nodes.
;
;   2018-06-08 : MGL. New keyword pref.
;
;-
pro red::fitgains, all = all $
                   , densegrid = densegrid $
                   , extra_nodes = extra_nodes $
                   , fit_reflectivity = fit_reflectivity $
                   , ifit = ifit $
                   , initcmap = initcmap $
                   , myg = myg $
                   , niter = niter $
                   , nosave = nosave $
                   , npar = npar $
                   , nthreads = nthreads $
                   , pref = pref $
                   , rebin = rebin $
                   , res = res $
                   , state = state $
                   , thres = thres $
                   , w0 = w0, w1 = w1 $
                   , x0 = x0, x1 = x1 $
                   , xl = xl, yl = yl 
  

  ;; Defaults
  if n_elements(niter) eq 0 then niter = 3L
  if n_elements(rebin) eq 0 then rebin = 100L
  
  ;; Prepare for logging (after setting of defaults). Set up a
  ;; dictionary with all parameters that are in use
  red_make_prpara, prpara, all
  red_make_prpara, prpara, densegrid
  red_make_prpara, prpara, extra_nodes
  red_make_prpara, prpara, fit_reflectivity
  red_make_prpara, prpara, ifit
  red_make_prpara, prpara, initcmap
  red_make_prpara, prpara, myg
  red_make_prpara, prpara, niter
  red_make_prpara, prpara, nosave
  red_make_prpara, prpara, npar
  red_make_prpara, prpara, nthreads
  red_make_prpara, prpara, pref
  red_make_prpara, prpara, rebin
  red_make_prpara, prpara, state
  red_make_prpara, prpara, thres
  red_make_prpara, prpara, w0
  red_make_prpara, prpara, w1
  red_make_prpara, prpara, x0
  red_make_prpara, prpara, x1
  red_make_prpara, prpara, xl
  red_make_prpara, prpara, yl

  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  outdir = self.out_dir + '/flats/'

  files = file_search(self.out_dir + 'flats/spectral_flats/*_filenames.txt' $
                      , count = Nfiles)

  print, inam + ' : Found '+strtrim(Nfiles, 2)+' files'

  ;; Select data
  selectionlist = red_strreplace(file_basename(files),'_filenames.txt','') ;strarr(Nfiles)

;  for ifile = 0L, Nfiles-1 do begin
;    spawn, 'cat ' + files[ifile], namelist
;    self -> extractstates, namelist, states
;    selectionlist[ifile] = states[0].prefilter + ' ' + states[0].cam_settings
;  endfor                        ; ifile
  
  case 1 of

    keyword_set(all) : begin
      Nselect = Nfiles
      sindx = indgen(Nselect)
    end

    n_elements(pref) ne 0 : begin
      sindx = where(strmatch(selectionlist, '*_'+pref), Nselect)
      if Nselect eq 0 then begin
        print, inam + ' : Keyword pref='+pref+' does not correspond to any output from prepflatcube.'
        print, files
        stop
      end
    end
    
    else : begin
      tmp = red_select_subset(selectionlist $
                              , qstring = 'Select data to be processed' $
                              , count = Nselect, indx = sindx)
    end

  endcase
  

;  if keyword_set(all) then begin
;    Nselect = Nfiles
;    sindx = indgen(Nselect)
;  endif else begin
;    tmp = red_select_subset(selectionlist $
;                            , qstring = 'Select data to be processed' $
;                            , count = Nselect, indx = sindx)
;  endelse

  for iselect = 0L, Nselect-1 do begin

    ifile = sindx[iselect]

    ;; Load data
    print
    print, file_basename(files[ifile], '_filenames.txt')
    spawn, 'cat ' + files[ifile], namelist
    Nwav = n_elements(namelist)
    
    ffile = red_strreplace(files[ifile], '_filenames.txt', '_flats_data.fits') ; Input flat intensities
    wfile = red_strreplace(files[ifile], '_filenames.txt', '_flats_wav.fits')  ; Input flat wavelengths
    sfile = red_strreplace(files[ifile], '_filenames.txt', '_fit_results.sav') ; Output save file
    pfile = red_strreplace(files[ifile], '_filenames.txt', '_fit_results.png') ; Plot file
    pfile0 = red_strreplace(files[ifile], '_filenames.txt', '_fit_results0.png') ; Plot file

    self -> extractstates, namelist, states
    outnames = self -> filenames('cavityflat', states)

    if Nwav eq 0 then continue
    if Nwav eq 1 then begin
      print, inam+' : Only one wavelength point.'
      print, inam+' : Will assume it is a continuum point and copy the ordinary flat.'

      if ~keyword_set(nosave) then begin

        print, inam + ' : Copying '+namelist[0]+' -> '+outnames[0]
        file_copy, namelist[0], outnames[0], /overwrite

        ;; Need to set res (the cavity error map?).
        dat = readfits(ffile)
        wav = float(readfits(wfile)) * 1e10 ; Must be in Å.
        dims = size(dat, /dim)
;        if keyword_set(fit_reflectivity) then npar_t = max([nparr,3]) else npar_t = max([nparr,2])
        npar_t = 1
        res = dblarr([npar_t, dims[1:2]])
        res[0,*,*] = dat

        ;; Need to set xl, yl (from red_get_imean):
;        imean = fltarr(dims[0])
;        yl = imean
        yl = [1.];median(dat)
        iwav = wav ;- median(res[1,*,*])
        xl = iwav
        
        ;res /= yl

        fit = {pars:res, yl:yl, xl:xl, oname:outnames}
        save, file = sfile, fit
        fit = 0B
      endif      
      
      continue
    endif

    ;; At this point, Nwav gt 1.

    pref = states[0].prefilter

    if n_elements(npar) eq 0 then begin
      ;; Add prefilters and nparr values here as we gain experience.
      case pref of
        '3934' : nparr = 3
        '3969' : nparr = 3
        '4862' : nparr = 3
        '6302' : nparr = 2
        '6563' : nparr = 5
        '8542' : nparr = 5
        else: nparr = 2         ; default
      endcase
    endif else nparr = npar
    ;; npar_t = npar + (keyword_set(fit_reflectivity) ? 3 : 2)
    if keyword_set(fit_reflectivity) then npar_t = max([nparr,3]) else npar_t = max([nparr,2])

    cub = readfits(ffile)
    wav = float(readfits(wfile)) * 1e10 ; Must be in Å.

    case 1 of 
      n_elements(w0) gt 0 and n_elements(w1) gt 0 : begin
        cub = (temporary(cub))[w0:w1,*,*]
        wav = wav[w0:w1]
        namelist = namelist[w0:w1]
      end
      n_elements(w0) gt 0 : begin
        cub = (temporary(cub))[w0:*,*,*]
        wav = wav[w0:*]
        namelist = namelist[w0:*]
      end
      n_elements(w1) gt 0 : begin
        cub = (temporary(cub))[0:w1,*,*]
        wav = wav[0:w1]
        namelist = namelist[0:w1]
      end
      else:
    endcase
    Nwav = n_elements(namelist) ; Nwav has to be adjusted if w0 or w1 were used.
      
    if n_elements(extra_nodes) gt 0 then begin
      if n_elements(myg) eq 0 then myg = wav
      myg = [extra_nodes, myg]
      myg = myg[sort(myg)]
    endif 

    dat = cub                   ; Why do we make a copy of cub?

    ;; Init output vars
    dims = size(dat, /dim)
    res = dblarr([npar_t, dims[1:2]])
    ratio = fltarr(dims)

    res[0,*,*] = total(dat,1) / nwav
    ;; Is this still needed?
    IF keyword_set(fit_reflectivity) THEN IF npar_t GT 3 THEN res[3:*,*,*] = 1.e-3

    ;; Init cavity map?   
    if(keyword_set(initcmap)) then begin
      print, inam + ' : Initializing cavity-errors with parabola-fit'
      res[1,*,*] = red_initcmap(wav, dat, x0 = x0, x1 = x1)
    endif

    cgwindow
    title = selectionlist[ifile] + ' npar='+strtrim(nparr, 2)

    ;; Loop niter
    for iiter = 0L, niter - 1 do begin

      ;;dum = 2
      ;;if(keyword_set(fit_reflectivity)) then dum += 1
      ;;  print, inam+ 'Normalizing polynomial coefs'
      ;;  for ii = dum, npar_t-1 do begin
      ;;     kk = median(reform(res[ii,*,*]))
      ;;     print, ' <C_'+string(ii-dum,format='(I0)')+'> = ', kk
      ;;     res[ii,*,*] -= kk 
      ;;  endfor
      
      ;; Get mean spectrum using Hermitian Spline
      yl = red_get_imean(wav, dat, res, npar_t, iiter $
                         , xl = xl, rebin = rebin, densegrid = densegrid, thres = thres $
                         , myg = myg, reflec = fit_reflectivity $
                         , title = title + ' iter='+strtrim(iiter, 2))

      if iiter eq 0 then cgcontrol, output = pfile0
      
      ;; Pixel-to-pixel fits using a C++ routine to speed-up things
      if iiter eq 0 then begin
        res1 = res[0:1,*,*]
        if(keyword_set(ifit)) then begin
          red_ifitgain, res1, wav, dat, xl, yl, ratio
        endif else red_cfitgain, res1, wav, dat, xl, yl, ratio, nthreads=nthreads
        res[0:1,*,*] = temporary(res1)
      endif else begin
        if keyword_set(fit_reflectivity) then $
           red_cfitgain2, res, wav, dat, xl, yl, ratio, pref, nthreads = nthreads $
        else begin
          if(keyword_set(ifit)) then begin
            red_ifitgain, res, wav, dat, xl, yl, ratio 
          endif else red_cfitgain, res, wav, dat, xl, yl, ratio, nthreads = nthreads
        endelse
      endelse

    endfor                      ; iiter

    yl = red_get_imean(wav, dat, res, npar_t, iiter $
                       , xl = xl, rebin = rebin, densegrid = densegrid $
                       , thres = thres, myg = myg, reflec = fit_reflectivity $
                       , title = title)

    ;; Create cavity-error-free flat (save in "ratio" variable)
    print, inam + ' : Recreating cavity-error-free flats ... ', FORMAT='(A,$)'

    for ii = 0L, nwav - 1 do ratio[ii,*,*] *= reform(res[0,*,*]) * $
       reform(red_get_linearcomp(wav[ii], res, npar_t, reflec = fit_reflectivity))

    print, 'done'
    
    ;; Save results
    if ~keyword_set(nosave) then begin

      for iwav = 0L, Nwav - 1 do begin

        output = reform(ratio[iwav,*,*])

        ;; Make FITS header
        head = red_readhead(namelist[iwav])
        check_fits, output, head, /UPDATE, /SILENT  
        fxaddpar, head, 'DATE', red_timestamp(/iso), 'UTC creation date of FITS header'
        fxaddpar, head, 'FILENAME', file_basename(outnames[iwav]), after = 'DATE'
        self -> headerinfo_addstep, head, prstep = 'Make cavity free flats' $
                                    , prproc = inam, prpara = prpara

        print, inam + ' : Saving file -> '+outnames[iwav]
        red_writedata, outnames[iwav], output, header = head, /overwrite

      endfor                    ; iwav
      
      fit = {pars:res, yl:yl, xl:xl, oname:outnames}
      save, file = sfile, fit
      fit = 0B

      ;; Save the plot made by red_get_imean
      cgcontrol, output = pfile

    endif

  endfor                        ; iselect

  print, inam + ' : Done!'
  print, inam + ' : Please check flats/spectral_flats/*fit_results.png.'

end
