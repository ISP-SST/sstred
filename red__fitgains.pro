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
;    densegrid : in 
;   
;   
;   
;    extra_nodes : in, optional, type=array
;
;       Make extra spline nodes at the specified wavelengths, given in
;       Ångström from the core point.
;   
;    fit_reflectivity : in, optional, type=boolean
;   
;       Fit reflectivity in addition to cavity map.
;   
;    initcmap : in, optional, type=boolean
;   
;       Initialize the cavity map.
;   
;    myg :  : in, optional, type=array
;   
;       Specify the wavelength nodes, alternative to using the data
;       tuning points.
;   
;    niter : in, optional, type=integer, default=3
;   
;       Number of iterations.
;   
;    npar : in, optional, type=integer, default=1
;   
;       Number of additional terms (default: 1 = linear)
;   
;    nosave : in, optional, type=boolean 
;   
;       Do not save the output.
;   
;    nthreads : in, optional, type=integer
;   
;       Number of threads used by red_cfitgain and red_cfitgain2.
;   
;    pref : in, optional, type=string
;
;       Process data for this prefilter.
;
;    rebin : in, optional, type=integer
;   
;   
;   
;    res : out, optional, type=array
;   
;       The results of the fit, including the cavity map.
;   
;    sig : 
;   
;   
;   
;    thres : in
;   
;   
;   
;    x0 : in 
;   
;   
;   
;    x1 : in 
;   
;   
;   
;    xl : in
;   
;   
;   
;    yl : in
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
pro red::fitgains, all = all $                             ;
                   , densegrid = densegrid $               ; Used in red_get_imean
                   , extra_nodes = extra_nodes $           ; 
                   , fit_reflectivity = fit_reflectivity $ ; Boolean
;                   , ifit = ifit $                         ; Not used
                   , initcmap = initcmap $                 ; Boolean, initialize cavity map
                   , myg = myg $                           ; Specify the wavelength nodes instead of using the data tuning points
                   , niter = niter $                       ; 
                   , nosave = nosave $                     ; 
                   , npar = npar $                         ; 
                   , nthreads = nthreads $                 ; 
                   , pref = pref $                         ; 
                   , rebin = rebin $                       ; Used in red_get_imean
                   , res = res $                           ; Result, output
;                   , state = state $                       ; Not used
                   , thres = thres $                       ; Used in red_get_imean
                   , w0 = w0, w1 = w1 $                    ; 
                   , x0 = x0, x1 = x1 $                    ; Used in red_initcmap
                   , xl = xl, yl = yl $                    ; 
                   , sig = sig                             ; Used like pref (by red_cfitgain) but is an array?


  ;; Defaults
  if n_elements(niter) eq 0 then niter = 3L
  if n_elements(rebin) eq 0 then rebin = 100L
  
  ;; Prepare for logging (after setting of defaults). Set up a
  ;; dictionary with all parameters that are in use
  red_make_prpara, prpara, all
  red_make_prpara, prpara, densegrid
  red_make_prpara, prpara, extra_nodes
  red_make_prpara, prpara, fit_reflectivity
;  red_make_prpara, prpara, ifit
  red_make_prpara, prpara, initcmap
  red_make_prpara, prpara, myg
  red_make_prpara, prpara, niter
  red_make_prpara, prpara, nosave
  red_make_prpara, prpara, npar
  red_make_prpara, prpara, nthreads
  red_make_prpara, prpara, pref
  red_make_prpara, prpara, rebin
;  red_make_prpara, prpara, state
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

  ;; Loop over selected prefilters
  for iselect = 0L, Nselect-1 do begin 

    ifile = sindx[iselect]

    ;; Load data
    print
    print, file_basename(files[ifile], '_filenames.txt')
    spawn, 'cat ' + files[ifile], namelist
    Nwav = n_elements(namelist)

    if Nwav eq 0 then continue  ; No data
    
    ffile = red_strreplace(files[ifile], '_filenames.txt', '_flats_data.fits')   ; Input flat intensities
    wfile = red_strreplace(files[ifile], '_filenames.txt', '_flats_wav.fits')    ; Input flat wavelengths
    sfile = red_strreplace(files[ifile], '_filenames.txt', '_fit_results.sav')   ; Output save file
    pfile = red_strreplace(files[ifile], '_filenames.txt', '_fit_results.pdf')   ; Plot file
    pfile0 = red_strreplace(files[ifile], '_filenames.txt', '_fit_results0.pdf') ; Plot file

    self -> extractstates, namelist, states
    outnames = self -> filenames('cavityflat', states)

    pref = states[0].prefilter

    if Nwav eq 1 then begin     ; Most likely CHROMIS Ca II continuum
      
      red_message, 'Only one wavelength point: '+states[0].tuning
      red_message, 'Will assume it is a continuum point so the  ordinary flat is free from cavity errors.'

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
        yl = [1.]               ;median(dat)
        iwav = wav              ;- median(res[1,*,*])  ; iwav is used as a loop variable later....
        xl = iwav
        
        ;res /= yl

        fit = {pars:res, yl:yl, xl:xl, oname:outnames}
        save, file = sfile, fit
        fit = 0B
      endif      
      
      continue
    endif

    ;; At this point, Nwav gt 1, so we have a scan.

    ;; Continued processing in three steps: 1. Fit cavity error and
    ;; reflectivity model; 2. Use the fit to make flats without cavity
    ;; errors. 3. Save the cavity error free flats 4. Save the fit;
    ;; Some wavelength points could be excluded from 1. but should
    ;; still be included in 2 and 3.
    
    if n_elements(npar) eq 0 then begin
      ;; Add prefilters and nparr values here as we gain experience.
      case pref of
        '3934' : nparr = 3
        '3969' : nparr = 3
        '4862' : nparr = 3
        '6302' : nparr = 2
        '6563' : nparr = 5
        '8542' : nparr = 5
        else : nparr = 2        ; default
      endcase
    endif else nparr = npar
    ;; npar_t = npar + (keyword_set(fit_reflectivity) ? 3 : 2)
    
    ;; Read the flats cube and the wavelength tunings
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

    orig_cub = cub              ; Preserve original flats cube for steps 2 and 3.
    orig_wav = wav
    orig_namelist = namelist
    orig_outnames = outnames

    ;; Step 1: fitting

    findx = indgen(Nwav)        ; Wavelength index

    ;; Take w0 and w1 keywords into account
    if n_elements(w0) eq 0 then w0 = 0
    if n_elements(w1) eq 0 then w1 = Nwav-1
    findx = where(findx ge w0 and findx le (w1 lt 0 ? Nwav+w1 : w1), Nfit $
                  , complement = lindx, Ncomplement = Nlinks)
    ;; findx for tunings to include in the fit and make cavityfree
    ;; flats for, lindx for tunings for which to make links to the
    ;; ordinary flats.
    if Nfit eq 0 then stop
    
;    findx = findx[w0:w1]

    ;; Info to keep for the linking
    lcub = cub[lindx,*,*,*]
    lwav = wav[lindx]
    lnamelist = namelist[lindx]
    loutnames = outnames[lindx]

    ;; Info restricted to the tunings used in the fit
    cub = (temporary(cub))[findx,*,*,*]
    wav = wav[findx]
    namelist = namelist[findx]
    outnames = outnames[findx]
;    cub = (temporary(cub))[w0:w1,*,*,*]
;    wav = wav[w0:w1]
;    namelist = namelist[w0:w1]
;    outnames = outnames[w0:w1]

    Nwav = n_elements(namelist) ; Nwav has to be adjusted if w0 or w1 were used.



    
    if keyword_set(fit_reflectivity) then npar_t = max([nparr,3]) else npar_t = max([nparr,2])
    npar_t_save = npar_t
    npar_one = npar_t-1 ;; parameters per lc without the cavity error
    
    npar_t = Nlc*(npar_t-1)+1 ;; all LCs see the same cavity map and are fitted at once.


;    ;; Do we have polcal data for this prefilter? If we do, we
;    ;; probably want to exclude polcal flats at that wavelength from
;    ;; the fit, as it is usually far outside the line.
;    psum_files = file_search('polcal_sums/*-R/cam*_'+pref+'_*pols.fits', count = Npolcal)
;    if Npolcal gt 0 then begin
;
;      self -> extractstates, psum_files[0], psum_state
;
;      stop
;      ;; Check if the polcal tuning is in the tunings. If it is, offer
;      ;; to remove it. Adjust findx (and the arrays modified for w0,w1
;      ;; above).
;      
;      red_message, 'Wavelength points for fitting (wav) : ['+strjoin(strtrim(wav,2),', ')+']'
;      diff = red_differential(wav)
;      diff = diff[1:*]
;      diff = diff/median(diff)
;      
;      if max(diff, maxloc) gt 5 then begin
;        case maxloc of
;          0 : begin
;            tmp = 'first'
;            default = '1-'+red_stri(n_elements(wav)-1)
;          end
;          n_elements(wav)-2 : begin
;            tmp = 'last'
;            default = '0-'+red_stri(n_elements(wav)-2)
;          end
;          else : tmp = ''
;        endcase
;        
;        if tmp ne '' then begin
;          print
;          red_message, ['The '+tmp+' wavelength point might be the polcal wavelengh' $
;                        , ' and could then be skipped.']
;          selectionlist = string(1000*wav, format = '(I5)')
;          tmp = red_select_subset(selectionlist, default = default, indx = sel $
;                                  , qstring = 'Which points do you want to keep?')
;
;          wav = wav[sel]
;          Nwav = n_elements(wav)
;          cub = cub[sel, *, *, *]
;        endif
;      endif
;
;    endif                       ; Npolcal

        
    if n_elements(extra_nodes) gt 0 then begin
          if n_elements(myg) eq 0 then myg = wav
          myg = [extra_nodes, myg]
          myg = myg[sort(myg)]
          red_message, 'Wavelength points myg : ['+strjoin(strtrim(myg,2),',')+']'
    endif

    dat = cub
    totdat = total(dat,1) / nwav
    
    ;; Make a gaintable from the sum of the dat frames and use it to
    ;; define the pixels we will use to do the fitting.
    
    mask = red_flat2gain(total(totdat,1)/Nlc) ne 0
    mindx = where(mask, Nmask)
    if Nmask eq 0 then stop

    dims = size(dat, /dim)    

    mdat = fltarr(dims[0], Nlc, 1, Nmask) ; Make a pseudo-2d array that will work with the fitting routines
    for i = 0, dims[0]-1 do for ilc=0,Nlc-1 do mdat[i,ilc, 0, *] = (reform(dat[i,ilc, *, *]))[mindx]
    
    ;; Init output vars
    res = dblarr([npar_t, dims[2:3]])
    ratio = fltarr(dims)

    mdims = size(mdat, /dim)
    mres = dblarr([npar_t, mdims[2:3]])
    mratio = fltarr(mdims)
    
    off = indgen(Nlc)*npar_t_save
    for ilc = 0, Nlc-1 do begin
      if(ilc ge 2) then off[ilc] -= (ilc-1)
    endfor

    ;; Init the gains to the mean of each lc state
    for ilc=0,Nlc-1 do mres[off[ilc],*,*] = (reform(totdat[ilc,*,*]))[mindx]
                                ;res[0,*,*] = totdat / nwav
    ;; Is this still needed?
;    IF keyword_set(fit_reflectivity) THEN IF npar_t GT 3 THEN res[3:*,*,*] = 1.e-3    
    ;;IF keyword_set(fit_reflectivity) THEN IF npar_t GT 3 THEN mres[3:*,*,*] = 1.e-3    ;; // JDLCR: fix this!!!!!

    
    ;; Init cavity map?   
    if(keyword_set(initcmap)) then begin
      print, inam + ' : Initializing cavity-errors with parabola-fit'
      mres[1,*,*] = red_initcmap(wav, total(mdat, 2)/Nlc, x0 = x0, x1 = x1)      
      ;;res[1,*,*] = red_initcmap(wav, dat, x0 = x0, x1 = x1)      
    endif

    cgwindow
    ;;  title = selectionlist[ifile] + ' npar='+strtrim(nparr, 2)
    title = pref + ' npar='+strtrim(nparr, 2)

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
        for ilc=1,Nlc-1 do mres[off[ilc],*,*] = mres1[1+ilc,*,*]
        
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

    ;; Step 2: Create cavity-error-free flat (save in "ratio" variable)

;    ;; For the following steps we need the original data
;    dat = orig_cub
;    wav = orig_wav
;    namelist = orig_namelist  
;    nwav = n_elements(wav)

    Nprogress = Nlc*Nwav
    iprogress = 0L
    tmp_par = dblarr(npar_one+1, ddim[1], ddim[2])
    for ilc=0,Nlc-1 do begin
      
      if ilc eq 0 then begin
        tmp_par[*,*,*] = mres[0:npar_one,*,*]
      endif else begin
        tmp_par[0,*,*] = mres[off[ilc],*,*]
        tmp_par[1,*,*] = mres[1,*,*]
        if(npar_one+1 gt 2) then tmp_par[2:*,*,*] = mres[off[ilc]+1:off[ilc]+npar_one-1,*,*]
      endelse
      
      for iwav = 0L,Nwav-1 do begin
        red_progressbar, iprogress,  Nprogress, 'Recreating cavity-error-free flats'
        iprogress++
        mratio[iwav,ilc,*,*] *= reform(tmp_par[0,*,*]) * $    
                                reform(red_get_linearcomp(wav[iwav], tmp_par, npar_one+1 $
                                                          , reflec = fit_reflectivity))
      endfor                    ; iwav
      
    endfor                      ; ilc
    ;;for ii = 0L, nwav - 1 do ratio[ii,*,*] *= reform(res[0,*,*]) * $    
    ;; reform(red_get_linearcomp(wav[ii], res, npar_t, reflec = fit_reflectivity))    

    ;; We want to save ratio and res so we need to fill them with data
    ;; from mratio and mres.
    tmp = fltarr(dims[2:3])    
    for ilc=0,Nlc-1 do for iwav = 0L,Nwav-1 do begin
      tmp[mindx] = reform(mratio[iwav, ilc, 0, *])
      ratio[iwav, ilc, *, *] = tmp
    endfor                      ; ilc
    
    tmp = dblarr(dims[2:3])        
    for ii = 0, npar_t-1 do begin
      tmp[mindx] = mres[ii, 0, *]
      res[ii, *, *] = tmp
    endfor
    
    ;; 3. Save cavity error free flats
    
    if ~keyword_set(nosave) then begin

      for iwav = 0L, Nwav - 1 do begin
        for ilc=0, Nlc-1 do begin
          output = reform(ratio[iwav,ilc,*,*,*])
          
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

            lc = 'lc'+string(ilc,format='(I1)')
            
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

        endfor                  ; ilc
      endfor                    ; iwav

      ;; 4. Save the fit
      fit = {pars:res, yl:yl, xl:xl, oname:outnames}
      save, file = sfile, fit
      fit = 0B

      ;; Save the plot made by red_get_imean
      cgcontrol, output = pfile

      ;; 5. Links to regular flats (for tunings not included in the
      ;; fit)
      if Nlinks gt 0 then begin
        red_message, 'Creating soft links to regular flats for tunings not included in the fit.'
        foreach iwav, lindx do begin ; Loop link indices
          for ilc=0, Nlc-1 do begin
            if Nlc eq 1 then begin
              oname = orig_outnames[iwav]
              fname = orig_namelist[iwav]
            endif else begin
              lc = 'lc'+string(ilc,format='(I1)')
              iname = red_strreplace(orig_namelist[iwav], 'lc0', lc)
              oname = red_strreplace(iname, '.flat', '_cavityfree.flat')
              fname = iname
            endelse
            ;; Make the link
            print, inam + ' : Linking file -> '+oname
            file_delete, oname, /allow_nonexistent ; file_link will not overwrite
            file_link, fname, oname
          endfor                ; ilc
        endforeach              ; iwav
      endif                     ; Nlinks
      
      
    endif                       ; nosave
    
    scrollwindow, /free, xs = dims[2], ys = dims[3] $
                  , title = 'Cavity map for '+selectionlist[ifile] + ' npar='+strtrim(nparr, 2)
    tvscl, red_histo_opt(reform(res[1,*,*]), 5e-3)

  endfor                        ; iselect

  print, inam + ' : Done!'
  print, inam + ' : Please check flats/spectral_flats/*fit_results.png and the displayed cavity map(s).'

end
