; docformat = 'rst'

;+
; Demodulate restored images for a particular scan and tuning state.
; 
; This method does (or at least is meant to do) the same as
; pol__demodulate.pro in the master (old CRISP) branch, except that it
; reads the inverse modulation matrices from disk rather than having
; pointers to them from within the pol class object. (We don't use a
; pol class in the new code base.) We rely on red_mozaik to handle
; both the inverse mod matrices and the image data the same way when
; it comes to clipping.
; 
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
;   2018-10-09 : MGL. First version, based on Jaime's pol::demodulate.
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
  camw = wbstates[0].camera
  camt = nbtstates[0].camera
  camr = nbrstates[0].camera

  ;; Geometrical distortion maps for the two NB cameras
  self -> getalignment, align=align, prefilters=prefilter
  indx = where(align.state2.camera eq camt, Nalign)
  case Nalign of
    0    : stop                 ; Should not happen!
    1    : amapt =      align[indx].map
    else : amapt = mean(align[indx].map, dim = 3)
  endcase
  amapt /= amapt[2, 2]          ; Normalize
  indx = where(align.state2.camera eq camr, Nalign)
  case Nalign of
    0    : stop                 ; Should not happen!
    1    : amapr =      align[indx].map
    else : amapr = mean(align[indx].map, dim = 3)
  endcase
  amapr /= amapr[2, 2]          ; Normalize

  
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
    red_append, xps,   fxpar(wbhead, 'XPOSURE')
    red_append, texps, fxpar(wbhead, 'TEXPOSUR')
    red_append, nsums, fxpar(wbhead, 'NSUMEXP')
  endfor                        ; ilc
  tbeg = min(tbegs)
  tavg = mean(tavgs)
  tend = max(tends)
  
  if n_elements(wcs) gt 0 then wcs.time = tavg

                                ;
  ;; Create arrays for the image mozaics
  img_t = fltarr(Nx, Ny, 4)
  img_r = fltarr(Nx, Ny, 4)

  ;; load files
  if ~keyword_set(noflat) then begin

    ;; Adapt this block to new code base!

    dims = size(immt, /dim)
    Nxx = dims[1]               ; Detector size
    Nyy = dims[2]
    utflat = fltarr(Nxx, Nyy)
    urflat = fltarr(Nxx, Nyy)
    tgain = fltarr([Nxx,Nyy,Nlc])
    rgain = fltarr([Nxx,Nyy,Nlc])
    
    ;; utflat and urflat are the "unpolarized" ("cavityfree") flats
    ;; for the current state (ingoring the LC state)
    self -> get_calib, nbtstates[0], cflatdata = utflat
    self -> get_calib, nbrstates[0], cflatdata = urflat
    
    for ilc = 0, Nlc-1 do begin

      ;; tgain and rgain are the gaintables for the current state,
      ;; including the LC states.

      self -> get_calib, nbtstates[ilc], gaindata = gaindata ; Read the gain
      tgain[0, 0, ilc] = gaindata
      
      self -> get_calib, nbrstates[ilc], gaindata = gaindata ; Read the gain
      rgain[0, 0, ilc] = gaindata

    endfor                      ; ilc
    
    print, inam + ' : Applying flat ratios to the inverse modulation matrix.'
    immt_dm = red_demat(tgain, utflat, immt)
    immr_dm = red_demat(rgain, urflat, immr)
    print, 'done'

  endif else begin
    immt_dm = immt
    immr_dm = immr
  endelse

  
  if keyword_set(smooth_by_subfield) then begin

    ;; Divide geometrically distortion-corrected modulation matrices
    ;; into momfbd-subfields, smooth with PSF from each subfield, make
    ;; mosaics of the result.
    
    immts = red_matrix2momfbd(timg, immt_dm, amap = amapt)
    immrs = red_matrix2momfbd(rimg, immr_dm, amap = amapr)

    mymt = fltarr(Nlc, Nstokes, Nx, Ny)
    mymr = fltarr(Nlc, Nstokes, Nx, Ny)

    ;; The red_mozaic calls below return matrices in the cropped size.
    ;; Do we need to correct their geometry as well at this point?

    ;; Mozaic images and demodulation matrix
    for ilc = 0L, Nlc-1 do begin
      img_r[*,*,ilc] = red_mozaic(rimg[ilc], /crop) * nbrfac
      img_t[*,*,ilc] = red_mozaic(timg[ilc], /crop) * nbtfac
      for istokes = 0L, Nstokes-1 do mymr[ilc,istokes,*,*] = red_mozaic(immrs[ilc,istokes], /crop) 
      for istokes = 0L, Nstokes-1 do mymt[ilc,istokes,*,*] = red_mozaic(immts[ilc,istokes], /crop) 
    endfor

    ;; if keyword_set(cmap) then cmap = (red_mozaic(red_conv_cmap(cmap,*self.timg[1])))[x0:x1,y0:y1]

  endif else begin

    stop

    ;; Smooth the modulation matrices with average PSF.


  endelse


  ;; Load the simultaneous WB images.  
  ;;img_wb = fltarr(dim[0], dim[1], Nlc)
  img_wb = fltarr(Nx, Ny, Nlc)
  for ilc = 0L, Nlc-1 do img_wb[*,*,ilc] = red_readdata(wfiles[ilc]) 

  ;; if keyword_set(cmap) then cmap = red_stretch(temporary(cmap), grid1)

  rest = fltarr(Nx, Ny, 4)
  resr = fltarr(Nx, Ny, 4)
  
  
  if n_elements(wbg) gt 0 then begin

    ;; Demodulate with destretching
    for ilc = 0L, Nlc-1 do begin
      grid = red_dsgridnest(wbg, img_wb[*,*,ilc], tiles, clips)
      for istokes = 0L, 3 do begin
        rest[*,*,istokes] += red_stretch(reform(mymt[ilc,istokes,*,*]) $
                                         * img_t[*,*,ilc], grid)
        resr[*,*,istokes] += red_stretch(reform(mymr[ilc,istokes,*,*]) $
                                         * img_r[*,*,ilc], grid)
      endfor                    ; istokes
    endfor                      ; ilc

  endif else begin

    ;; No global WB image, demodulate without destretching.
    for ilc = 0L, Nlc-1 do begin
      for istokes = 0L, 3 do begin
        rest[*,*,istokes] += reform(mymt[ilc,istokes,*,*]) * img_t[*,*,ilc]
        resr[*,*,istokes] += reform(mymr[ilc,istokes,*,*]) * img_r[*,*,ilc]
      endfor                    ; istokes
    endfor                      ; ilc

  endelse

  
  ;; These are now Stokes cubes!
  img_t = temporary(rest)
  img_r = temporary(resr)



  dum = size(img_t,/dim)
  drm = dum / 8.0 
  xx0 = round(drm[0] - 1)
  xx1 = round(dum[0] - drm[0] - 1)
  yy0 = round(drm[1] - 1)
  yy1 = round(dum[1] - drm[1] - 1)
  
  mediant = median(img_t[xx0:xx1,yy0:yy1,0])
  medianr = median(img_r[xx0:xx1,yy0:yy1,0])
  aver = (mediant + medianr) / 2.
  sct = aver / mediant
  scr = aver / medianr
  
  res = (sct * (img_t) + scr * (img_r)) / 2.
  
  print, inam + ' : Combining data from transmitted and reflected camera'
  print, '   -> Average Intensity = '  + red_stri(aver)
  print, '   -> Tcam scale factor -> ' + red_stri(sct) + ' (after '+red_stri(nbtfac)+')'
  print, '   -> Rcam scale factor -> ' + red_stri(scr) + ' (after '+red_stri(nbrfac)+')'


  if 0 then begin

    meant=mean(img_t[xx0:xx1,yy0:yy1,0])
    meanr=mean(img_r[xx0:xx1,yy0:yy1,0])
    mediant=median(img_t[xx0:xx1,yy0:yy1,0])
    medianr=median(img_r[xx0:xx1,yy0:yy1,0])

    cghistoplot,img_r[xx0:xx1,yy0:yy1,0],xrange=[2.6,4.3]*1e-8
    cghistoplot,img_t[xx0:xx1,yy0:yy1,0],/oplot,color='blue' 
    cgoplot,[1,1]*medianr,[0,20000],color='red',/over,line=2
    cgoplot,[1,1]*mediant,[0,20000],color='blue',/over,line=2
    cgoplot,[1,1]*meanr,[0,20000],color='red',/over          
    cgoplot,[1,1]*meant,[0,20000],color='blue',/over         
    
    aver = (meant + meanr) / 2.
    sct = aver / meant
    scr = aver / meanr
    
    print, inam + ' : Combining data from transmitted and reflected camera'
    print, '   -> Average Intensity = '  + red_stri(aver)
    print, '   -> Tcam scale factor -> ' + red_stri(sct) + ' (after '+red_stri(nbtfac)+')'
    print, '   -> Rcam scale factor -> ' + red_stri(scr) + ' (after '+red_stri(nbrfac)+')'

    aver2 = (mediant + medianr) / 2.
    sct2 = aver2 / mediant
    scr2 = aver2 / medianr
    
    print, inam + ' : Combining data from transmitted and reflected camera'
    print, '   -> Average Intensity = '  + red_stri(aver2)
    print, '   -> Tcam scale factor -> ' + red_stri(sct2) + ' (after '+red_stri(nbtfac)+')'
    print, '   -> Rcam scale factor -> ' + red_stri(scr2) + ' (after '+red_stri(nbrfac)+')'
  
  endif
  
  ;; Telescope model 
  line = (strsplit(wbstates[0].fpi_state,'_',/extract))[0]
  print, inam+' : Detected spectral line -> '+line



  red_logdata, self.isodate, tavg, azel = azel,  tilt = tilt
  year = (strsplit(self.isodate, '-', /extract))[0]
  mtel = red_telmat(line, {TIME:tavg, AZ:azel[0], ELEV:azel[1], TILT:tilt} $
                    , /no_zero, year=year)
  imtel = invert(mtel) 
  imtel /= imtel[0]

  ;; Apply the telescope Muller matrix
  res1 = fltarr(Nx, Ny, Nstokes)
  for j=0, Nstokes-1 do for i=0, Nstokes-1 do $
     res1[*, *, j] += res[*, *, i] * imtel[i, j]
  res = temporary(res1)
  
  
  
  if keyword_set(nosave) then return
  
  ;; Save result as a fitscube with all the usual metadata.
  
  file_mkdir, outdir

  dims = [Nx, Ny, 1, Nstokes, 1] 
  res = reform(res, dims)

  ;; Make header
  hdr = red_readhead(nbrstates[0].filename)
  check_fits, res, hdr, /update, /silent
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
  red_fitsdelkeyword, hdr, 'CADENCE' ; Makes no sense to keep?
  
  
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

  ;; Add statistics metadata
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
