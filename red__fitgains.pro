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
;    densegrid  : 
;   
;   
;   
;    extra_nodes : in, optional, type=array
;
;       Make extra spline nodes at the specified wavelengths, given in
;       Ångström from the core point.
;
;    initcmap  : 
;   
;   
;   
;    myg  : 
;   
;   
;   
;    niter  : 
;   
;   
;   
;    npar  : 
;   
;   
;   
;    nosave  : 
;   
;   
;   
;    pref : in, optional, type=string
;
;       Select prefilter to process.
;
;    rebin  : 
;   
;   
;   
;    res  : 
;   
;   
;   
;    state  : 
;   
;   
;   
;    thres  : 
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
;    xl  : 
;   
;   
;   
;    yl  : 
;   
;   
;   
;    w0 : in, type=integer, default=0
;   
;      The index of the lowest wavelength point to include in the fit.
;      Use this if necessary to deselect blue points that are too far
;      from each other.
;   
;    w1 : in, type=integer, default="Npoints-1"
;   
;      The index of the highest wavelength point to include in the
;      fit. Use this if necessary to deselect red points that are too
;      far from each other.
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
;   2022-09-01 : MGL. Automatically mask bad pixels and pixels without
;                light.
;
;   2023-11-02 : JdlCR. Modifications for multiple LC fitting with
;                a single cavity map value for the new demodulation
;                and flat-fielding scheme.
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
                   , xl = xl, yl = yl $
                   , sig = sig


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
  inam = red_subprogram(/low, calling = inam1)                              

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
    pfile = red_strreplace(files[ifile], '_filenames.txt', '_fit_results.pdf') ; Plot file
    pfile0 = red_strreplace(files[ifile], '_filenames.txt', '_fit_results0.pdf') ; Plot file

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
        dims = dims[-2:-1]
        
;        if keyword_set(fit_reflectivity) then npar_t = max([nparr,3]) else npar_t = max([nparr,2])
        npar_t = 1
        res = dblarr([npar_t, dims[0:1]])
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
    

    cub = readfits(ffile)
    wav = float(readfits(wfile)) * 1e10 ; Must be in Å.
    
    
    dim = size(cub, /dim)
    sig1 = dblarr(dim[0])+1.0
    
    if(n_elements(sig) ne 0) then begin
      nsig = n_elements(sig)
      nmax = min([dim[0],nsig])
      sig1[0:nmax-1] = sig[0:nmax-1]
    endif
    sig = sig1

    
    if(n_elements(dim) eq 3) then begin
      cub = reform(temporary(cub), [dim[0], 1, dim[1], dim[2]])
      dim = size(cub, /dim)
    endif
    Nlc = dim[1]

    if keyword_set(fit_reflectivity) then npar_t = max([nparr,3]) else npar_t = max([nparr,2])
    npar_t_save = npar_t
    npar_one = npar_t-1 ;; parameters per lc without the cavity error
    
    npar_t = Nlc*(npar_t-1)+1 ;; all LCs see the same cavity map and are fitted at once.
    
    case 1 of 
      n_elements(w0) gt 0 and n_elements(w1) gt 0 : begin
        cub = (temporary(cub))[w0:w1,*,*,*]
        wav = wav[w0:w1]
        namelist = namelist[w0:w1]
      end
      n_elements(w0) gt 0 : begin
        cub = (temporary(cub))[w0:*,*,*,*]
        wav = wav[w0:*]
        namelist = namelist[w0:*]
      end
      n_elements(w1) gt 0 : begin
        cub = (temporary(cub))[0:w1,*,*,*]
        wav = wav[0:w1]
        namelist = namelist[0:w1]
      end
      else:
    endcase
    Nwav = n_elements(namelist) ; Nwav has to be adjusted if w0 or w1 were used.

    red_message, 'Wavelength points wav : ['+strjoin(strtrim(wav,2),', ')+']'

    if n_elements(extra_nodes) gt 0 then begin
      if n_elements(myg) eq 0 then myg = wav
      myg = [extra_nodes, myg]
      myg = myg[sort(myg)]
      red_message, 'Wavelength points myg : ['+strjoin(strtrim(myg,2),',')+']'
    endif

    dat = cub                   ; Why do we make a copy of cub?
    totdat = total(dat,1) / nwav
    
    ;; Make a gaintable from the sum of the dat frames and use it to
    ;; define the pixels we will use to do the fitting.
    
    mask = red_flat2gain(total(totdat,1)/nlc) ne 0
    mindx = where(mask, Nmask)
    if Nmask eq 0 then stop

    dims = size(dat, /dim)    

    mdat = fltarr(dims[0], Nlc, 1, Nmask) ; Make a pseudo-2d array that will work with the fitting routines
    for i = 0, dims[0]-1 do for ss=0,Nlc-1 do mdat[i,ss, 0, *] = (reform(dat[i,ss, *, *]))[mindx]
    
    ;; Init output vars
    res = dblarr([npar_t, dims[2:3]])
    ratio = fltarr(dims)

    mdims = size(mdat, /dim)
    mres = dblarr([npar_t, mdims[2:3]])
    mratio = fltarr(mdims)
    
    off = indgen(Nlc)*npar_t_save
    for ss = 0, Nlc-1 do begin
      if(ss ge 2) then off[ss] -= (ss-1)
    endfor

    ;; Init the gains to the mean of each lc state
    for ss=0,Nlc-1 do mres[off[ss],*,*] = (reform(totdat[ss,*,*]))[mindx]
                                ;res[0,*,*] = totdat / nwav
    ;; Is this still needed?
;    IF keyword_set(fit_reflectivity) THEN IF npar_t GT 3 THEN res[3:*,*,*] = 1.e-3    
    ;;IF keyword_set(fit_reflectivity) THEN IF npar_t GT 3 THEN mres[3:*,*,*] = 1.e-3    ;; // JDLCR: fix this!!!!!

    
    ;; Init cavity map?   
    if(keyword_set(initcmap)) then begin
      print, inam + ' : Initializing cavity-errors with parabola-fit'
      mres[1,*,*] = red_initcmap(wav, total(mdat, 2)/nlc, x0 = x0, x1 = x1)      
      ;;res[1,*,*] = red_initcmap(wav, dat, x0 = x0, x1 = x1)      
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
                                ;yl = red_get_imean(wav, dat, res, npar_t, iiter $
                                ;                   , xl = xl, rebin = rebin, densegrid = densegrid, thres = thres $
                                ;                   , myg = myg, reflec = fit_reflectivity $
                                ;                   , title = title + ' iter='+strtrim(iiter, 2))
      yl = red_get_imean(wav, mdat, mres, npar_t, iiter $
                         , xl = xl, rebin = rebin, densegrid = densegrid, thres = thres $
                         , myg = myg, reflec = fit_reflectivity $
                         , title = title + ' iter='+strtrim(iiter, 2))

      if iiter eq 0 then cgcontrol, output = pfile0
      ddim = size(mres, /dim)

      ;; Pixel-to-pixel fits using a C++ routine to speed-up things
      if iiter eq 0  and ~keyword_set(initcmap) then begin
        ;;mres1 = mres[0:1,*,*]
        ;;ddim = size(mres, /dim)
        mres1 = dblarr(Nlc+1,ddim[1], ddim[2])
        mres1[0:1,*,*] = mres[0:1,*,*]
        for dd=1,Nlc-1 do mres1[dd+1,*,*] = mres[off[dd],*,*]
               
        red_cfitgain, mres1, wav, mdat, xl, yl, mratio, sig, nthreads=nthreads          
        ;;res[0:1,*,*] = temporary(res1)        
        mres[0:1,*,*] = mres1[0:1,*,*]
        for ss=1,nlc-1 do mres[off[ss],*,*] = mres1[1+ss,*,*]
        
      endif else begin
        if keyword_set(fit_reflectivity) then $
           ;;red_cfitgain2, res, wav, dat, xl, yl, ratio, pref, nthreads = nthreads $           
           red_cfitgain2, mres, wav, mdat, xl, yl, mratio, pref, nthreads = nthreads $           
        else begin
          red_cfitgain, mres, wav, mdat, xl, yl, mratio, sig, nthreads = nthreads            
        endelse
      endelse

    endfor                      ; iiter

                                ;yl = red_get_imean(wav, dat, res, npar_t, iiter $
                                ;                   , xl = xl, rebin = rebin, densegrid = densegrid $
                                ;                   , thres = thres, myg = myg, reflec = fit_reflectivity $
                                ;                   , title = title)

    yl = red_get_imean(wav, mdat, mres, npar_t, iiter $
                       , xl = xl, rebin = rebin, densegrid = densegrid $
                       , thres = thres, myg = myg, reflec = fit_reflectivity $
                       , title = title)

    ;; Create cavity-error-free flat (save in "ratio" variable)
    print, inam + ' : Recreating cavity-error-free flats ... ', FORMAT='(A,$)'

    
    tmp_par = dblarr(npar_one+1, ddim[1], ddim[2])
    for ss=0,Nlc-1 do begin
      
      if(ss eq 0) then begin
        tmp_par[*,*,*] = mres[0:npar_one,*,*]
      endif else begin
        tmp_par[0,*,*] = mres[off[ss],*,*]
        tmp_par[1,*,*] = mres[1,*,*]

        if(npar_one+1 gt 2) then begin
          tmp_par[2:*,*,*] = mres[off[ss]+1:off[ss]+npar_one-1,*,*]
        endif
      endelse
      
      for ii = 0L, nwav - 1 do begin
        
        mratio[ii,ss,*,*] *= reform(tmp_par[0,*,*]) * $    
                             reform(red_get_linearcomp(wav[ii], tmp_par, npar_one+1, reflec = fit_reflectivity))
      endfor
    endfor
      ;;for ii = 0L, nwav - 1 do ratio[ii,*,*] *= reform(res[0,*,*]) * $    
    ;; reform(red_get_linearcomp(wav[ii], res, npar_t, reflec = fit_reflectivity))    

    print, 'done'
    
    ;; We want to save ratio and res so we need to fill them with data
    ;; from mratio and mres.
    tmp = fltarr(dims[2:3])    
    for ss=0,Nlc-1 do for iwav = 0L, nwav - 1 do begin
      tmp[mindx] = reform(mratio[iwav, ss, 0, *])
      ratio[iwav, ss, *, *] = tmp
    endfor
    
    
    tmp = dblarr(dims[2:3])        
    for ii = 0, npar_t-1 do begin
      tmp[mindx] = mres[ii, 0, *]
      res[ii, *, *] = tmp
    endfor
    
    ;; Save results
    if ~keyword_set(nosave) then begin

      for iwav = 0L, Nwav - 1 do begin
        for ss=0, Nlc-1 do begin
          output = reform(ratio[iwav,ss,*,*,*])
          
          if(Nlc eq 1) then begin
         
            ;; Make FITS header
            head = red_readhead(namelist[iwav])
            check_fits, output, head, /UPDATE, /SILENT  
            fxaddpar, head, 'DATE', red_timestamp(/iso), 'UTC creation date of FITS header'
            fxaddpar, head, 'FILENAME', file_basename(outnames[iwav]), after = 'DATE'
            self -> headerinfo_addstep, head, prstep = 'SPECTRAL-COORDINATE-DISTORTION-CORRECTION' $
                                        , prproc = inam, prpara = prpara
            
            print, inam + ' : Saving file -> '+outnames[iwav]
            red_writedata, outnames[iwav], output, header = head, /overwrite
          endif else begin

            lc = 'lc'+string(ss,format='(I1)')
            
            iname = red_strreplace(namelist[iwav], 'lc0', lc)
            oname = red_strreplace(iname, '.flat', '_cavityfree.flat')

            head = red_readhead(iname)
     
            check_fits, output, head, /UPDATE, /SILENT  
            fxaddpar, head, 'DATE', red_timestamp(/iso), 'UTC creation date of FITS header'
            fxaddpar, head, 'FILENAME', file_basename(outnames[iwav]), after = 'DATE'
            self -> headerinfo_addstep, head, prstep = 'SPECTRAL-COORDINATE-DISTORTION-CORRECTION' $
                                        , prproc = inam, prpara = prpara
            
            print, inam + ' : Saving file -> '+oname
            red_writedata, oname, output, header = head, /overwrite
            
            
          endelse
        endfor                  ; iwav
      endfor
      
      fit = {pars:res, yl:yl, xl:xl, oname:outnames}
      save, file = sfile, fit
      fit = 0B

      ;; Save the plot made by red_get_imean
      cgcontrol, output = pfile

    endif
    
    scrollwindow, /free, xs = dims[2], ys = dims[3] $
                  , title = 'Cavity map for '+selectionlist[ifile] + ' npar='+strtrim(nparr, 2)
    tvscl, red_histo_opt(reform(res[1,*,*]), 5e-3)

  endfor                        ; iselect

  print, inam + ' : Done!'
  print, inam + ' : Please check flats/spectral_flats/*fit_results.png and the displayed cavity map(s).'

end
