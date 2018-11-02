; docformat = 'rst'

;+
; Demodulate restored images for a particular scan and tuning state.
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, Institute for Solar Physics
; 
; 
; :Returns:
; 
; 
; :Params:
; 
;    outname
; 
; :Keywords:
; 
;     clips : in, optional, type=array
;
;       Used to compute stretch vectors for the wideband alignment.
;
;    nthreads : in, optional, type=integer
;
;       The number of threads to use.
; 
;    nbrfac : in, optional, type=float
;
;       A scaling factor to apply to NBR data to get the correct units.
;
;    nbrstates : in, type=structarr
;   
;       An array of structs with the states of the input NBR images.
; 
;    nbtfac : in, optional, type=float
;
;       A scaling factor to apply to NBT data to get the correct units.
;
;    nbtstates : in, type=structarr
;   
;       An array of structs with the states of the input NBT images.
;
;     tiles : in, optional, type=array
;
;       Used to compute stretch vectors for the wideband alignment.
; 
;    wbg : in, optional, type=image
;
;      The global WB image to use as destretching reference. If this
;      keyword is not used, no destretching will be done.
; 
;    wbstates : in, optional, type=structarr
;   
;       An array of structs with the states of the input WB images. If
;       provided, will use them to stretch the NB images when
;       combining them.
; 
; 
; :History:
; 
;   2018-10-09 : MGL. First version, based on Jaime's
;                red::polarim and pol::demodulate.
; 
;-
pro crisp::demodulate, outname, immr, immt $
                       , clips = clips $
                       , nbrfac = nbrfac $
                       , nbrstates = nbrstates $
                       , nbtfac = nbtfac $
                       , nbtstates = nbtstates $
                       , newflats = newflats $
                       , no_ccdtabs = no_ccdtabs $
                       , nthreads = nthreads $
                       , overwrite = overwrite $
                       , smooth_by_subfield = smooth_by_subfield $
                       , tiles = tiles $
                       , units = units $
                       , wbg = wbg $
                       , wcs = wcs $
                       , wbstates = wbstates

  ;; Feature requested by Vasco: Support polarimetry with 5173.
  ;; Problem: there is no telmat for 5173 but it seems to work well
  ;; with 5250 telmat. Possible solution: add an optional keyword in
  ;; red_telmat allowing selection of the filter for which telmat is
  ;; to be used (default the actual filter of the data). This keyword
  ;; would have to be propagated from make_nb_cube, through
  ;; demodulate, to red_telmat.

  ;; What to do with the newflats keyword? The polarim method uses it
  ;; in call to red_getstates_polarim.
  
  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  if file_test(outname) and ~keyword_set(overwrite) then begin
    print, inam + ' : Already exists: '+outname
    return
  endif
  
  red_make_prpara, prpara, clips         
  red_make_prpara, prpara, nbrfac 
  red_make_prpara, prpara, nbtfac
  red_make_prpara, prpara, newflats
  red_make_prpara, prpara, no_ccdtabs 
  red_make_prpara, prpara, smooth_by_subfield
  red_make_prpara, prpara, tiles 
  
  ;; Defaults of keywords
  if n_elements(Nthreads) eq 0 then Nthreads = 1
  if n_elements(units) eq 0 then units = 'D.N.'
  
  outdir = file_dirname(outname)
  
  ;; Should be four images from each camera (four LC states)
  Nstokes = 4
  Nlc = 4
  if n_elements(nbrstates) ne Nlc then stop
  if n_elements(nbtstates) ne Nlc then stop
  if n_elements(wbstates)  ne Nlc then stop

;  Nelements = 16
  
  rfiles = nbrstates.filename
  tfiles = nbtstates.filename
  wfiles = wbstates.filename

  for ilc = 0, Nlc-1 do begin
    red_append, rimg, momfbd_read(rfiles[ilc])
    red_append, timg, momfbd_read(tfiles[ilc])
    red_append, wimg, momfbd_read(wfiles[ilc])
  endfor

  ;; FOV is given in the momfbd output, use it to crop the modulation matrices.
;  mr = momfbd_read(tfiles[0], /names) ; Use /names to avoid reading
;  the data parts
  mr = rimg[0]
  x0 = mr.roi[0] + mr.margin
  x1 = mr.roi[1] - mr.margin
  y0 = mr.roi[2] + mr.margin
  y1 = mr.roi[3] - mr.margin
  Nx = x1 - x0 + 1
  Ny = y1 - y0 + 1
  
  prefilter = wbstates[0].prefilter
;  nbtdetector = nbtstates[0].detector
;  nbrdetector = nbrstates[0].detector
  camw = wbstates[0].camera
  camt = nbtstates[0].camera
  camr = nbrstates[0].camera

  ;; Geometrical distortion maps for the two NB cameras
  self -> getalignment, align=align, prefilters=prefilter
  indx = where(align.state2.camera eq camt, Nalign)
  case Nalign of
    0    : stop                 ; Should not happen!
    1    : amapt = invert(      align[indx].map           )
    else : amapt = invert( mean(align[indx].map, dim = 3) )
  endcase
  indx = where(align.state2.camera eq camr, Nalign)
  case Nalign of
    0    : stop                 ; Should not happen!
    1    : amapr = invert(      align[indx].map           )
    else : amapr = invert( mean(align[indx].map, dim = 3) )
  endcase

  ;; Get some header info
  
  for ilc = 0, Nlc-1 do begin
    wbhead = red_readhead(wfiles[ilc])
    red_fitspar_getdates, wbhead $
                          , date_beg = date_beg $
                          , date_avg = date_avg $
                          , date_end = date_end 
    red_append, tbegs, red_time2double((strsplit(date_beg,'T',/extract))[1])
    red_append, tavgs, red_time2double((strsplit(date_avg,'T',/extract))[1])
    red_append, tends, red_time2double((strsplit(date_end,'T',/extract))[1])
    red_append, xps, fxpar(wbhead, 'XPOSURE')
    red_append, texps, fxpar(wbhead, 'TEXPOSUR')
    red_append, nsums, fxpar(wbhead, 'NSUMEXP')
  endfor                        ; ilc
  tbeg = min(tbegs)
  tavg = mean(tavgs)
  tend = max(tends)
  
  if n_elements(wcs) gt 0 then wcs.time = tavg

;  ;;tmp = red_mozaic(*self.timg[0])
;  ;;dim = red_getborder(tmp, x0, x1, y0, y1)
;  ;;dim = size(tmp, /dim)
;  x0 = self.x0
;  x1 = self.x1
;  y0 = self.y0
;  y1 = self.y1
;  if keyword_set(noclip) then begin
;    tmp = size(red_mozaic((*self.timg[1])), /dim)
;    x0 = 0L
;    x1 = tmp[0] - 1
;    y0 = 0L
;    y1 = tmp[1] - 1
;  endif
;  
;  dim = [x1 - x0 + 1, y1 - y0 + 1]
                                ;
;; Create arrays for the image mozaics
  img_t = fltarr(Nx, Ny, 4)
  img_r = fltarr(Nx, Ny, 4)

  if keyword_set(smooth_by_subfield) then begin

    ;; Divide modulation matrices into momfbd-subfields, smooth with
    ;; PSF from each subfield, make mosaics of the result.
    
    immts = red_matrix2momfbd(timg, immt)
    immrs = red_matrix2momfbd(rimg, immr)

    mymt = fltarr(Nlc, Nstokes, Nx, Ny)
    mymr = fltarr(Nlc, Nstokes, Nx, Ny)
  
    ;; Mozaic images and demodulation matrix
    for ilc = 0L, Nlc-1 do begin
      img_r[*,*,ilc] = red_mozaic(rimg[ilc], /crop) * nbrfac
      img_t[*,*,ilc] = red_mozaic(timg[ilc], /crop) * nbtfac
      for istokes = 0L, Nstokes-1 do mymr[ilc,istokes,*,*] = red_mozaic(immrs[ilc,istokes], /crop) 
      for istokes = 0L, Nstokes-1 do mymt[ilc,istokes,*,*] = red_mozaic(immts[ilc,istokes], /crop) 
    endfor

;    if keyword_set(cmap) then cmap = (red_mozaic(red_conv_cmap(cmap,*self.timg[1])))[x0:x1,y0:y1]

  endif else begin

    ;; Smooth the modulation matrices with average PSF.

  endelse


  ;; Load the simultaneous WB images.  
  ;;img_wb = fltarr(dim[0], dim[1], Nlc)
  img_wb = fltarr(Nx, Ny, Nlc)
  for ilc = 0L, Nlc-1 do img_wb[*,*,ilc] = red_readdata(wfiles[ilc]) 

;  if keyword_set(cmap) then cmap = red_stretch(temporary(cmap), grid1)

  rest = fltarr(Nx, Ny, 4)
  resr = fltarr(Nx, Ny, 4)
  
  
  if n_elements(wbg) gt 0 then begin

    ;; Demodulate with destretching
;    stop
;    grid = fltarr(2, 2, Nlc)
;    for ilc = 0, Nlc-1 do begin
;      grid[0, 0, ilc] = 
;    end
;    grid0 = red_dsgridnest(wbg, img_wb[*,*,0], tiles, clip)
;    grid1 = red_dsgridnest(wbg, img_wb[*,*,1], tiles, clip)
;    grid2 = red_dsgridnest(wbg, img_wb[*,*,2], tiles, clip)
;    grid3 = red_dsgridnest(wbg, img_wb[*,*,3], tiles, clip)
    
    for ilc = 0L, Nlc-1 do begin
      grid = red_dsgridnest(wbg, img_wb[*,*,ilc], tiles, clips)
      for istokes = 0L, 3 do begin
        
        rest[*,*,istokes] += red_stretch(reform(mymt[ilc,istokes,*,*]) * img_t[*,*,ilc], grid)
        resr[*,*,istokes] += red_stretch(reform(mymr[ilc,istokes,*,*]) * img_r[*,*,ilc], grid)

;        rest[*,*,istokes] = red_stretch(reform(mymt[0,istokes,*,*]) * img_t[*,*,0], grid0) + $
;                            red_stretch(reform(mymt[1,istokes,*,*]) * img_t[*,*,1], grid1) + $
;                            red_stretch(reform(mymt[2,istokes,*,*]) * img_t[*,*,2], grid2) + $
;                            red_stretch(reform(mymt[3,istokes,*,*]) * img_t[*,*,3], grid3)
;        
;        resr[*,*,istokes] = red_stretch(reform(mymr[0,istokes,*,*]) * img_r[*,*,0], grid0) + $
;                            red_stretch(reform(mymr[1,istokes,*,*]) * img_r[*,*,1], grid1) + $
;                            red_stretch(reform(mymr[2,istokes,*,*]) * img_r[*,*,2], grid2) + $
;                            red_stretch(reform(mymr[3,istokes,*,*]) * img_r[*,*,3], grid3) 

      endfor                    ; istokes
    endfor                      ; ilc
    
  endif else begin

    ;; No WBG image given, demodulate without destretching.

    for ilc = 0L, Nlc-1 do begin
      for istokes = 0L, 3 do begin

        rest[*,*,istokes] += reform(mymt[ilc,istokes,*,*]) * img_t[*,*,ilc]
        resr[*,*,istokes] += reform(mymr[ilc,istokes,*,*]) * img_r[*,*,ilc]

;        rest[*,*,istokes] = reform(mymt[0,istokes,*,*]) * img_t[*,*,0] + $
;                            reform(mymt[1,istokes,*,*]) * img_t[*,*,1] + $
;                            reform(mymt[2,istokes,*,*]) * img_t[*,*,2] + $
;                            reform(mymt[3,istokes,*,*]) * img_t[*,*,3]
;      
;        resr[*,*,istokes] = reform(mymr[0,istokes,*,*]) * img_r[*,*,0] + $
;                            reform(mymr[1,istokes,*,*]) * img_r[*,*,1] + $
;                            reform(mymr[2,istokes,*,*]) * img_r[*,*,2] + $
;                            reform(mymr[3,istokes,*,*]) * img_r[*,*,3] 

      endfor                    ; istokes
    endfor                      ; ilc

  endelse
  

  ;; These are now Stokes cubes!
  img_t = temporary(rest)
  img_r = temporary(resr)



  dum = size(img_t,/dim)
  drm = dum / 8.0 
  xx0 = drm[0] - 1
  xx1 = dum[0] - drm[0] - 1
  yy0 = drm[1] - 1
  yy1 = dum[1] - drm[1] - 1
  
  aver = (mean(img_t[xx0:xx1,yy0:yy1,0]) + mean(img_r[xx0:xx1,yy0:yy1,0])) / 2.
  sct = aver / mean(img_t[xx0:xx1,yy0:yy1,0])
  scr = aver / mean(img_r[xx0:xx1,yy0:yy1,0])
  
  res = (sct * (img_t) + scr * (img_r)) / 2.
  
  print, inam + ' : Combining data from transmitted and reflected camera'
  print, '   -> Average Intensity = '  + red_stri(aver)
  print, '   -> Tcam scale factor -> ' + red_stri(sct)
  print, '   -> Rcam scale factor -> ' + red_stri(scr)
  
  ;; Average obs. time (used for demodulation) 
  time_obs = 0.d0
  for ii = 0L, 3 do time_obs += red_time2double((timg[ii]).time)
  time_obs = red_time2double(time_obs * 0.25d0, /dir)
  
;  if n_elements(ext_time) gt 0 then begin
;    iscan = long(self.scan)
;    if n_elements(ext_time) gt iscan + 1 then time_obs = ext_time[iscan]
;  endif
  
  ;; telescope model
  
  line = (strsplit(wbstates[0].fpi_state,'_',/extract))[0]
  print, inam+' : Detected spectral line -> '+line



  ;; self.telog is the name of the turret log file? Replace accessing
  ;; this with a call to red_logdata?
  if file_test(self.telog) then begin

    red_logdata, self.isodate, tavg, azel = azel,  tilt = tilt
    year = (strsplit(self.isodate, '-', /extract))[0]
    mtel = red_telmat(line, {TIME:tavg, AZ:azel[0], ELEV:azel[1], TILT:tilt} $
                      , /no_zero, year=year)

;    ;; mtel =
;    ;; sst_mueller_all((*self.timg[0]).date, time_obs, self.telog, line)
;    if n_elements(mdate) eq 0 then mdate = strjoin(strsplit((timg[0]).date,'-/.',/extra),'/')
;    year = long((strsplit(mdate,'/',/extract))[0])
;    
;    telpos = red_read_azel( self.telog,mdate)
;    print, inam + ' : time_obs = '+time_obs
;    mtel = red_telmat(line, telpos, time_obs, /no_zero, year=year)
  
    imtel = invert(mtel) 
    imtel /= imtel[0]

    ;; Apply the matrix

    res1 = fltarr(Nx, Ny, Nstokes)
    for j=0, Nstokes-1 do begin
      for i=0, Nstokes-1 do begin
        res1[*, *, j] += res[*, *, i] * imtel[i, j]
      endfor                    ; i
    endfor                      ; j
    res = temporary(res1)
  endif else print, inam + ' : WARNING, SST position LOG not found -> telescope polarization not corrected!!!!'

  

  
  ;; Save result
  
;  if keyword_set(savecams) then begin
;    odir = file_dirname(self.tfiles[0]) + '/demodulated_cameras/'
;    file_mkdir, odir
;    
;    tcam = (strsplit(file_basename(self.tfiles[0]), '.',/extract))[0]
;    rcam = (strsplit(file_basename(self.rfiles[0]), '.',/extract))[0]
;
;    ofil = tcam + '.'+self.state + '.f0'
;    res1 = fltarr(dim[0], dim[1],4)
;    for j=0, 3 do for i=0, 3 do res1[*, *, j] += img_t[*, *, i] * imtel[i, j]
;    fzwrite,  temporary(res1) * sct, odir + ofil,' '
;    
;    ofil = rcam + '.'+self.state + '.f0'
;    res1 = fltarr(dim[0], dim[1],4)
;    for j=0, 3 do for i=0, 3 do res1[*, *, j] += img_r[*, *, i] * imtel[i, j]
;    fzwrite,  temporary(res1) * scr, odir + ofil,' '
;
;    return
;    
;  endif


  if keyword_set(nosave) then return
  
  file_mkdir, outdir
  
  ;; Save as a fitscube with all the usual metadata. Could be read
  ;; into crispex?
  
;  print, inam + ' : saving file -> '+ outdir + outname
;  head = 'TIME_OBS='+time_obs+' DATE_OBS='+timg[0].date
;  fzwrite, res, outdir + outname, head+''

  dims = [Nx, Ny, 1, Nstokes, 1] 
  res = reform(res, dims)

  ;; Make header
  hdr = red_readhead(nbrstates[0].filename)
  check_fits, res, hdr, /UPDATE, /SILENT
  red_fitsaddkeyword, anchor = anchor, hdr, 'FILENAME', file_basename(outname)
  red_fitsaddkeyword, anchor = anchor, hdr, 'OBS_HDU', 1

  
  ;; Add info about this step
  self -> headerinfo_addstep, hdr $
                              , prstep = 'Demodulate' $
                              , prpara = prpara $
                              , prproc = inam

  
  ;; New date, at better position
  red_fitsaddkeyword, anchor = anchor, after = 'SOLARNET', /force, hdr $
                      , 'DATE', red_timestamp(/iso), 'Creation UTC date of FITS header '

  ;; Remove some irrelevant keywords
  red_fitsdelkeyword, hdr, 'TAB_HDUS' ; No tabulated headers in summed file
  red_fitsdelkeyword, hdr, 'FRAMENUM' ; No particular frame number
  red_fitsdelkeyword, hdr, 'SCANNUM'  ; No particular scan number
  red_fitsdelkeyword, hdr, 'CADENCE'  ; Makes no sense to keep?
  
  
  ;; Set various keywords in hdr. DATE* and statistics are useful.
  ;; Also make sure we implement the Stokes dimension properly!
  red_fitsaddkeyword, anchor = anchor, hdr, 'BUNIT', units, 'Units in array'
  red_fitsaddkeyword, anchor = anchor, hdr, 'BTYPE', 'Intensity', 'Type of data in array'
  red_fitsaddkeyword, anchor = anchor, hdr $
                      , 'DATE-BEG', self.isodate + 'T' + red_timestring(tbeg) $
                      , 'Start time of combined observation'
  red_fitsaddkeyword, anchor = anchor, hdr $
                      , 'DATE-AVG', self.isodate + 'T' + red_timestring(tavg) $
                      , 'Average time of combined observation'
  red_fitsaddkeyword, anchor = anchor, hdr $
                      , 'DATE-END', self.isodate + 'T' + red_timestring(tend) $
                      , 'End time of combined observation'
  red_fitsaddkeyword, anchor = anchor, hdr, 'XPOSURE', total(xps)
  red_fitsaddkeyword, anchor = anchor, hdr, 'TEXPOSUR', mean(texps)
  red_fitsaddkeyword, anchor = anchor, hdr, 'NSUMEXP', round(total(nsums))

  ;; Get the frame numbers involved
  for ilc = 0, Nlc-1 do begin
    nhdr = red_readhead(nbrstates[ilc].filename)
    red_append, fnumsum, red_expandrange(fxpar(nhdr, 'FNUMSUM'))
  endfor                        ; ilc
  fnumsum = fnumsum[sort(fnumsum)]
  red_fitsaddkeyword, anchor = anchor, hdr $
                      , 'FNUMSUM', red_collapserange(fnumsum,ld='',rd='') $
                      , 'List of frame numbers used'
  red_fitsaddkeyword, anchor = anchor, hdr $
                      , 'CAMERA', nbrstates[0].camera+','+nbtstates[0].camera
  red_fitsaddkeyword, anchor = anchor, hdr $
                      , 'DETECTOR', nbrstates[0].detector+','+nbtstates[0].detector



  ;; Add "global" metadata
  red_metadata_restore, anchor = anchor, hdr
  
  ;; Statistics
  for istokes = 0, Nstokes-1 do begin
    red_append, statistics, red_image_statistics_calculate(res[*, *, 0, istokes, 0])
  endfor                        ; istokes
  Nbins = 2L^16                 ; Use many bins!
  binsize = (max(double(statistics.datamax)) - min(double(statistics.datamin))) / (Nbins - 1.)
  hist = lonarr(Nbins)
  for istokes = 0, Nstokes-1 do begin
    ;; Accumulate the histogram
    hist += histogram(res[*, *, 0, istokes, 0] $
                      , min = min(statistics.datamin) $
                      , max = max(statistics.datamax) $
                      , Nbins = Nbins, /nan)
  endfor                        ; istokes
  cubestats = red_image_statistics_combine(statistics $
                                           , hist = hist $
                                           , binsize = binsize)
  
  self -> fitscube_initialize, outname, hdr, lun, fileassoc, dims 
  for istokes = 0, Nstokes-1 do begin
    self -> fitscube_addframe, fileassoc, res[*, *, 0, istokes, 0] $
                               , istokes = istokes
    
  endfor                        ; istokes
  

  ;; Close the file
  self -> fitscube_finish, lun, wcs = wcs
  
  anchor = 'DATE-END'
  for itag = n_tags(statistics[0])-1, 0, -1 do begin
    itagc = where((tag_names(statistics[0]))[itag] eq tag_names(cubestats), Nmatch)
    if Nmatch eq 1 then $
       self -> fitscube_addvarkeyword, outname $
                                       , (tag_names(statistics[0]))[itag] $
                                       , statistics.(itag) $
                                       , anchor = anchor $
                                       , keyword_value = cubestats.(itagc) $
                                       , axis_numbers = 4
  end
  
  
  
;    if keyword_set(cmap) then begin
;      odir = file_dirname(self.tfiles[0]) + '/cavity_map/'
;      file_mkdir, odir
;      ofile = 'cmap.' + self.state + '.f0'
;      print, inam + 'saving cavity map -> '+odir+ofile
;      fzwrite, float(cmap), odir+ofile, head + ''
;    endif
  
  
end
