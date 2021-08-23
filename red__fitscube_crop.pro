; docformat = 'rst'

;+
; Crop a fitscube while maintaining meta data.
; 
; :Categories:
;
;    SST pipeline
; 
; :Author:
; 
;    Mats Löfdahl, Institute for Solar Physics
; 
; :Params:
; 
;    infile : in, type=string
; 
;      The name of the file containing the fitscube to be cropped.
; 
; 
; :Keywords:
;
;    centered : in, optional, type=boolean
;
;      Set this to get a centered FOV. Specify the size with the size
;      keyword or 
;
;    roi : in, out, optional, type=array
;
;      Give the ROI corner coordinates [llx,urx,lly,ury] to bypass
;      interactive FOV selection. If the FOV is selected interactively
;      and/or specified with the size and centered keywords the roi is
;      returned in this keyword.
;
;    iscan : in, optional, type=integer
;
;      Select image to display in xroi GUI by setting iscan, ituning,
;      and/or iscan. Indices not specified signals a dimension to sum
;      over. 
;
;    istokes : in, optional, type=integer
;
;      Select image to display in xroi GUI by setting iscan, ituning,
;      and/or iscan. Indices not specified signals a dimension to sum
;      over. 
;
;    ituning : in, optional, type=integer
;
;      Select image to display in xroi GUI by setting iscan, ituning,
;      and/or iscan. Indices not specified signals a dimension to sum
;      over. 
;
;    overwrite : in, optional, type=boolean
;
;      Set this to overwrite older files.
; 
;    nospectral : in, optional, type=boolean
;   
;      If set, don't make any spectral, "flipped" data cubes.
; 
;    outfile : in, optional, type=string, default="Based on infile"
; 
;      The name of the file in which to store the cropped cube.
; 
;    size : in, optional, type="integer or array"
; 
;      Set this to a FOV size to only select center of FOV
;      interactively. If an array, it is interpreted as sx=size[0],
;      sy=size[1]. If a single integer, sx=sy=size.
; 
;    xc : in, out, optional, type=integer
;   
;      X pixel coordinate of center of FOV. 
;
;    yc : in, out, optional, type=integer
;   
;      Y pixel coordinate of center of FOV. 
; 
; 
; :History:
; 
;    2017-11-08 : MGL. First version.
; 
;    2017-11-10 : MGL. Copy (and modify when needed) extensions.
;                 Renamed keyword corners to roi and changed the order
;                 of elements.
; 
;    2018-01-12 : MGL. Use subroutine red_fitscube_getframe rather
;                 than method fitscube_getframe.
; 
; 
;-
pro red::fitscube_crop, infile $
                        , centered = centered $
                        , iscan = iscan $
                        , istokes = istokes $
                        , ituning = ituning $
                        , nospectral = nospectral $
                        , outfile = outfile $
                        , overwrite = overwrite $
                        , roi = roi $
                        , size = size $
                        , xc = xc $
                        , yc = yc 

  ;; Do we want to "crop" in time, tuning, stokes as well?
  
  inam = red_subprogram(/low, calling = inam1)

  ;; Check the input file
  if ~file_test(infile) then begin
    print, inam + ' : Input file does not exist: '+infile
    retall
  endif
  if strmatch(infile, '*_sp.fits') then begin
    print, inam + ' : Do not crop a spectral fitscube, do the image file instead and have it flipped.'
    print, infile
    retall
  endif
  
  ;; Make default output file name if needed
  if n_elements(outfile) eq 0 then begin
    if strmatch(infile, '*_im.fits') then begin
      old_ending = '_im.fits'
      new_ending = '_cropped_im.fits'
    endif else begin
      old_ending = '.fits'
      new_ending = '_cropped.fits'
    endelse
    outfile = file_dirname(infile) + '/' + file_basename(infile, old_ending) + new_ending
  endif
  
  ;; Check the output file
  if file_test(outfile) then begin
    if file_same(outfile, infile) then begin
      print, inam + ' : Input and output file names refer to the same file.'
      print, infile
      print, outfile
      retall
    endif
    if ~keyword_set(overwrite) then begin
      print, inam + ' : Output file exists, use /overwrite to reuse it.'
      print, outfile
      retall
    endif
  endif

  ;; Make prpara
  red_make_prpara, prpara, roi 
  red_make_prpara, prpara, infile 
  red_make_prpara, prpara, nospectral 
  red_make_prpara, prpara, outfile 
  red_make_prpara, prpara, overwrite 
  red_make_prpara, prpara, size

  ;; Read header,
  in_hdr = headfits(infile)
  naxis = fxpar(in_hdr, 'NAXIS*')  ; Original dimensions

  if naxis[2] eq 1 and naxis[3] eq 1 then nospectral = 1
  
  
  ;; Make sure we have the FOV information in terms of the roi
  ;; array.
  case 1 of
    n_elements(iscan) ne 0 and n_elements(istokes) ne 0 and n_elements(ituning) ne 0 : begin
      red_fitscube_getframe, infile, dispim $
                             , ituning = ituning $
                             , istokes = istokes $
                             , iscan = iscan
    end
    n_elements(iscan) ne 0 and n_elements(istokes) ne 0 : begin
      dispim = 0.0
      for ituning = 0, naxis[2]-1 do begin
        red_fitscube_getframe, infile, thisframe $
                               , ituning = ituning $
                               , istokes = istokes $
                               , iscan = iscan
        dispim += thisframe
      endfor                    ; ituning
    end
    n_elements(iscan) ne 0 and n_elements(ituning) ne 0 : begin
      dispim = 0.0
      for istokes = 0, naxis[3]-1 do begin
        red_fitscube_getframe, infile, thisframe $
                               , ituning = ituning $
                               , istokes = istokes $
                               , iscan = iscan
        dispim += thisframe
      endfor                    ; istokes
    end
    n_elements(istokes) ne 0 and n_elements(ituning) ne 0 : begin
      dispim = 0.0
      for iscan = 0, naxis[4]-1 do begin
        red_fitscube_getframe, infile, thisframe $
                               , ituning = ituning $
                               , istokes = istokes $
                               , iscan = iscan
        dispim += thisframe
      endfor                    ; iscan
    end
    keyword_set(nospectral) : begin
      dispim = 0.0
      for iscan = 0, naxis[4]-1 do begin
        red_fitscube_getframe, infile, thisframe $
                               , ituning = ituning $
                               , istokes = istokes $
                               , iscan = iscan
        dispim += thisframe
      endfor                    ; iscan
    end
    else : begin
      ;; If the keywords gave no clue, default to the first frame. 
      red_fitscube_getframe, infile, dispim
    end
  endcase
  tmp = red_crop(dispim $
                 , centered = centered $
                 , roi = roi $
                 , /nocropping $ ; Don't crop, just get the roi
                 , size = size $
                 , xc = xc $
                 , yc = yc $
                )

  ;; FOV size
  Sx = roi[1] - roi[0] + 1
  Sy = roi[3] - roi[2] + 1


  ;; Adapt header to new FOV.
  hdr = in_hdr
  self -> headerinfo_addstep, hdr $
                              , prstep = 'CROPPING' $
                              , prpara = prpara $
                              , prproc = inam

  ;; Initialize outfile.
  dims = naxis
  dims[0:1] = [Sx, Sy]
  self -> fitscube_initialize, outfile, hdr, lun, fileassoc, dims 

  ;; Loop through data frames: read, crop, fitscube_addframe.
  Nframes = round(product(dims[2:*]))
  for iframe = 0, Nframes-1 do begin
    red_progressbar, iframe, Nframes, /predict, 'Copying cropped frames'
    red_fitscube_getframe, infile, frame, iframe = iframe
    frame = frame[roi[0]:roi[1], roi[2]:roi[3]]
    self -> fitscube_addframe, fileassoc, temporary(frame), iframe = iframe
  endfor                        ; iframe

  ;; Read WCS extension, adapt roi coordinates to new FOV.
  red_fitscube_getwcs, infile, coordinates = old_wcs_coord, distortions = wcs_dist
  wcs_coord = old_wcs_coord
  hpl_dims = size(wcs_coord.hplt, /dim)
  xx = roi[[[0,1],[0,1]]]/float(Naxis[0]-1)
  yy = roi[[[2,2],[3,3]]]/float(Naxis[1]-1)
  Ntun = hpl_dims[2]
  if n_elements(hpl_dims) le 3 then Nscans = 0 else Nscans = hpl_dims[3]
  for i_wave = 0, Ntun-1 do for i_time = 0, Nscans-1 do begin
    wcs_coord[i_wave, i_time].hpln = interpolate(old_wcs_coord[i_wave, i_time].hpln, xx, yy)
    wcs_coord[i_wave, i_time].hplt = interpolate(old_wcs_coord[i_wave, i_time].hplt, xx, yy)
  endfor                        ; iwave, itime

  
  ;; Use fitscube_finish to close the cube and possibly make a flipped
  ;; version.
;  if keyword_set(nospectral) then begin
  self -> fitscube_finish, lun, wcs = wcs_coord
;  endif else begin
;    self -> fitscube_finish, lun, flipfile = flipfile, wcs = wcs_coord
;  endelse
  
  if size(wcs_dist,/n_dim) gt 0 then begin
    cmap = reform(red_crop(wcs_dist.wave, roi = roi))
    red_fitscube_addcmap, outfile, reform(cmap, roi[1]-roi[0]+1, roi[3]-roi[2]+1, 1, 1, Nscans)
  endif

  
  ;; Find all extensions and copy them to the new file(s), except the
  ;; already copied WCS coordinate extensions.

  ;; First variable keywords. In addition to copying the extension,
  ;; also the main header is updated so do this with fitscube_addvar.
  var_keys = red_fits_var_keys(in_hdr)
  for ikey = 0, n_elements(var_keys)-1 do begin
    self -> fitscube_addvarkeyword, outfile, var_keys[ikey], old_file = infile
;    if n_elements(flipfile) ne 0 then $
;       self -> fitscube_addvarkeyword, flipfile, var_keys[ikey], old_file = infile, /flip
  endfor
  
  fits_open, infile, fcb
  free_lun, fcb.unit

  for iext = 1, fcb.nextend do begin ; 

    ;; Some extensions are already taken care of, continue with the next one.
    if fcb.extname[iext] eq 'Main' then continue              ; Main
    if fcb.extname[iext] eq 'WCS-TAB' then continue           ; WCS coordinates
    if fcb.extname[iext] eq 'WCSDVARR' then continue          ; WCS distortions
    if strmatch(fcb.extname[iext], 'VAR-EXT-*') then continue ; Variable keywords

    if fcb.extname[iext] eq 'MWCINFO' then begin
      print, inam + ' : Copying and cropping extension '+ fcb.extname[iext]

      ;; Read the extension
      fxbopen, bunit, infile, 'MWCINFO', bbhdr
      fxbreadm, bunit, row = 1 $
                , ['ANG', 'CROP', 'FF', 'GRID', 'ND', 'SHIFT', 'TMEAN', 'X01Y01'] $
                ,  rANG,   rCROP,   rFF,   rGRID,   rND,  rSHIFT,   rTMEAN,   rX01Y01
      ;; Note that the strarr wfiles cannot be read by fxbreadm!
      fxbread, bunit, rWFILES, 'WFILES', 1
      fxbclose, bunit
      wfiles_col = where(strtrim(fxpar(bbhdr,'TTYPE*'),2) eq 'WFILES')+1
      
      ;; Change rX01Y01 = [x0,x1,y0,y1] due to cropping
      cX01Y01 = rX01Y01[[0, 0, 2, 2]] + roi

      ;; Write the extension
      fxbcreate, bunit, outfile, bbhdr
      fxbwritm, bunit $
                , ['ANG', 'CROP', 'FF', 'GRID', 'ND', 'SHIFT', 'TMEAN', 'X01Y01'] $
                ,  rANG,  rCROP,  rFF,  rGRID,  rND,  rSHIFT,  rTMEAN,  cX01Y01
      fxbwrite, bunit, rWFILES, wfiles_col, 1
      fxbfinish, bunit
      
    endif
    
    ;; Copy other extensions
    case fcb.xtension[iext] of
      'BINTABLE' : begin
        print, inam + ' : Copying extension '+ fcb.extname[iext]
        red_fits_copybinext, infile, outfile, fcb.extname[iext]
;        if n_elements(flipfile) ne 0 then $
;           red_fits_copybinext, infile, flipfile, fcb.extname[iext]
      end
      'IMAGE' : begin
        if fcb.extname[iext] eq 'WBIMAGE' then begin
          extim = readfits(infile, ehdr, ext = iext)
          ;; Crop also the image extension
          extim = extim[roi[0]:roi[1], roi[2]:roi[3]]
          check_fits, extim, ehdr, /update
          fxaddpar, ehdr, 'DATE', red_timestamp(/utc, /iso)
          ;; Should add a PRPARA processing step for cropping in ehdr?
          ;; Or is it enough that the main header shows the cropping?
          writefits, outfile, extim, ehdr, /append
        endif else stop
      end
      else : begin
        ;; Extensions that we don't know how to handle:
        print, inam + ' : Can only copy binary extensions.'
        print, fcb.xtension[iext]
        stop
      end
    endcase
    
  endfor                        ; iext

  if ~keyword_set(noflip) then begin
    ;; Make a flipped version
    print, 'Flip it!'
    red_fitscube_flip, outfile $
                       , flipfile = flipfile $
                       , overwrite = overwrite
  endif

  print, 'Cropped cube written to '+outfile
  if n_elements(flipfile) ne 0 then $
     print, 'Flipped version written to '+flipfile
  
end
