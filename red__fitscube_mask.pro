; docformat = 'rst'

;+
; Define and apply a mask on the FOV of a fitscube and then crop to fit.
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
; :Params:
; 
;    infile : in, type=string
; 
;      The path to the input fitscube.
; 
; 
; :Keywords:
;
;    ituning : in, optional, type=integer
;
;      The tuning index of the Stokes image displayed when defining
;      the mask. If not given, the core tuning will be used.
;
;    noflip : in, optional, type=boolean
;   
;      If set, don't make any spectral, "flipped" data cubes.
; 
;    overwrite : in, type=boolean
;   
;      Overwrite an existing file.
; 
;    no_qhull : in, type=boolean
; 
;      Instead of computing the mask from the convex hull of several
;      ROIs, just use the union of the ROIs.
; 
;    tag : in, type=string, default='masked'
; 
;      If not overwriting the original file, add this tag to the
;      original filename.
;
;    no_crop : in, type=boolean
;
;      If set, don't crop the FOV.
;
;    include_I : in, type=boolean
;
;      If set, use mask also for stokes I and crop NB and WB cubes.
;
;    mask : in, out, type=array
;
;      FOV mask to be reused with other cubes (particularly WB cube). 
; 
; :History:
; 
;    2025-11-11 . MGL. First version.
;
;    2026-03-24 . OA. Added keywords, 'no_crop', 'include_I', 'mask'.
;                     Added recursive call for WB cube.
;
;-
pro red::fitscube_mask, infile $
                        , ituning = ituning $
                        , noflip = noflip $
                        , no_qhull = no_qhull $
                        , include_I = include_I $
                        , no_crop = no_crop $
                        , no_wb_mask = no_wb_mask $
                        , mask = mask $
                        , nomissing_nans = nomissing_nans $
                        , overwrite = overwrite $
                        , tag = tag

  if n_elements(tag) eq 0 then tag = 'masked'
  outfile = red_strreplace(infile, '_im.fits', '_' + tag + '_im.fits')

  if ~keyword_set(overwrite) && file_test(outfile) then begin
    red_message, 'Output file already exists. Use other tag or use /overwrite.'
    return
  endif
  
  indir = file_dirname(infile)
  tmpfile = indir + '/tmp'+cgTimeStamp(Random_Digits=5, /UTC)+'.fits'

  hdr = headfits(infile)
  steps = fxpar(hdr, 'PRSTEP*')
  naxis = fxpar(hdr, 'NAXIS*')

  Nx = naxis[0]
  Ny = naxis[1]
  Ntuning = naxis[2]
  Nstokes = naxis[3]
  Nscans = naxis[4]

  croostalk_done = 0B
  padding_done = 0B
  for istep=0, n_elements(steps)-1 do begin
    step = steps[istep]
    if strmatch(step,'*CROSSTALK*') then crosstalk_done = 1B
    if strmatch(step,'*PADDING*') then begin
      padding_done = 1B
      nomissing_nans = 1B
    endif
  endfor

  if Nstokes gt 1 then begin
    istokes = Nstokes-1
    if crosstalk_done and padding_done then $
      backgr = !Values.F_NaN $
    else $
      backgr = 0.
  endif else begin
    istokes = 0
    include_I = 1B
    if padding_done then $
      backgr = !Values.F_NaN $
    else $
      backgr = 0.
  endelse
  if keyword_set(include_I) then begin
    no_crop = 0B
    if Ntuning gt 1 then no_wb_mask = 0B
  endif

  if n_elements(ituning) eq 0 then begin
    red_message, 'Find core tuning point'
    red_fitscube_statistics, infile, frame_statistics, cube_statistics
    mn = min(frame_statistics[*,0,0].datamedn, ituning)
  endif
  
  ;; Base mask on total Stokes V for the selected tuning point 
  dispim = 0.0
  for iscan = 0, Nscans-1 do begin
    red_fitscube_getframe, infile, thisframe $
                           , ituning = ituning $
                           , istokes = istokes $
                           , iscan = iscan
    dispim += thisframe
  endfor                        ; iscan
  origmask = finite(dispim)
  dispim = bytscl(red_histo_opt(dispim))

  if n_elements(mask) eq 0 then begin
    red_message, 'Define the mask'

    if keyword_set(no_qhull) then begin
      XROI, dispim, Regions_Out=ROIs, /Block, title = 'Select areas.'
      mask = 0*origmask
      for iroi = 0, n_elements(ROIs)-1 do $
        mask or= ROIs[iroi] -> ComputeMask(Dimensions=[Nx, Ny], Mask_Rule=2) 
      mask and= origmask
      Obj_Destroy, ROIs
    endif else begin
      red_message, ['The XROI ellipse tool sometimes "remembers" points from when you expand' $
                      , 'and then contract a region. You will therefore have the opportunity' $
                      , 'to check the region and rerun XROI if you are not satisfied.' $
                     ]
      repeat begin
        XROI, dispim, Regions_Out=ROIs, /Block, title = 'Select areas.'
        undefine, roix, roiy
        for iroi = 0, n_elements(ROIs)-1 do begin
          indx = where((*rois[iroi].data)[0,*] ne 0 and (*rois[iroi].data)[1,*] ne 0)
          help, roix, roiy, indx, reform((*rois[iroi].data)[0,indx]), reform((*rois[iroi].data)[1,indx])
          red_append, roix, reform((*rois[iroi].data)[0,indx])
          red_append, roiy, reform((*rois[iroi].data)[1,indx])
        endfor               ; iroi
        cgimage, dispim, /keep, /save
        cgplot, roix, roiy, color = 'blue', /over, thick = 3
        qhull, roix, roiy, Tr, bounds = bounds
        ;; Determine the vertices order from the t output variable from
        ;; QHULL, and print the vertices order
        verts=INTARR(N_ELEMENTS(tr[0,*])+1)
        verts[0]=tr[0,0]
        verts[1]=tr[1,0]
        for i=1,N_ELEMENTS(tr[0,*])-1 do begin
          id=WHERE(tr[0,*] eq verts[i])
          verts[i+1]=tr[1,id]
        endfor
        ;;print," Polygon vertices order: ", verts
        ;; Order the points in x1 and y1 vectors based on vertices order
        roix2=roix[verts]
        roiy2=roiy[verts]
        cgplot,roix2, roiy2,color='red',/over, thick = 3
        ROI=IDLanROI()
        ROI.SetProperty,Data=TRANSPOSE([[roix2],[roiy2]])
        mask = ROI.ComputeMask(Dimensions=[Nx, Ny], Mask_Rule=2)
        Obj_Destroy, ROI
        Obj_Destroy, ROIs
        s = ''
        read, 'Mask ok [Y/n]?', s
        done = s eq '' || strupcase(strmid(s, 0, 1)) ne 'N'
           
      endrep until done
    endelse
  endif

  ;; Copy the original cube to a temporary file name and then
  ;; overwrite the frames with the masked frames. This way we maintain
  ;; all the metadata.
  file_copy, infile, tmpfile
  Nframes = long(Nscans)*long(NStokes)*long(Ntuning)
  indx = where(mask eq 0)
  for iscan=0, Nscans-1 do begin
    for ituning=0,Ntuning-1 do begin
      for istokes = 0, NStokes-1 do begin
        iframe = long(ituning) + long(istokes)*Ntuning $
          + long(iscan)*Ntuning*Nstokes
        red_progressbar, iframe, Nframes, 'Copying frames'
        red_fitscube_getframe, infile, thisframe, iframe = iframe
        if ~keyword_set(include_I) then begin
          if istokes gt 0 then $
            thisframe[indx] = backgr
        endif else $
          thisframe[indx] = backgr
        red_fitscube_addframe, tmpfile, thisframe, iframe = iframe 
      endfor
    endfor
  endfor                 

  ;; Add info about this step
  hdr = headfits(tmpfile)
  inam = red_subprogram(/low, calling = inam1)
  red_make_prpara, prpara, roix2
  red_make_prpara, prpara, roiy2
  red_make_prpara, prpara, no_qhull
  red_make_prpara, prpara, include_I
  red_make_prpara, prpara, no_crop
  red_make_prpara, prpara, nomissing_nans
  red_make_prpara, prpara, noflip
  red_make_prpara, prpara, infile
  red_make_prpara, prpara, outfile
  red_make_prpara, prpara, overwrite
  self -> headerinfo_addstep, hdr $
     , prstep = 'MASKING' $
     , prpara = prpara $
     , prproc = inam
  red_fitscube_newheader, tmpfile, hdr

  if ~keyword_set(nomissing_nans) then begin
    ;; Set padding pixels to missing-data, i.e., NaN.
    self -> fitscube_missing, tmpfile $
       , /noflip $
       , missing_type = 'nan' 
  endif

  ;; Then run fitscube_crop on the temporary file with a sensible
  ;; output file name. Either with a standard tag or possibly
  ;; overwriting the original file. Use the bb as input to
  ;; fitscube_crop.
  if ~keyword_set(no_crop) then begin
    bb = red_boundingbox(mask)
    self -> fitscube_crop, tmpfile $
                            , nospectral = noflip $
                            , outfile = outfile $
                            , /overwrite $
                           , roi = bb[[0, 2, 1, 3]]
    file_delete, tmpfile
  endif else file_move, tmpfile, outfile

  ;; Do we need to generate also a matching WB file? Could use the
  ;; same bb but not care about the masking --> more context. Use same
  ;; tag in the new wb file name, or overwrite if that's what you want.
  
  print, outfile

  if ~keyword_set(noflip) then begin
    red_fitscube_flip, outfile, flipfile = flipfile, overwrite = overwrite
    print, flipfile
  endif

  if keyword_set(include_I) and ~keyword_set(no_wb_mask) then begin
    indx = where(strmatch(hdr, '*Align reference:*'), Nwhere)
    key = strtrim((strsplit(hdr[indx[-1]], '=', /extract))[0], 2)
    wcfile = red_strreplace(fxpar(hdr, key), 'Align reference: ', '')
    if ~file_test(wcfile) then begin
      print, 'WB file ',wcfile,' does not exist.'
      stop
    endif
    self->fitscube_mask, wcfile, /noflip, /no_wb_mask $
                          , ituning = 0 $
                          , mask = mask $
                          , nomissing_nans = nomissing_nans $
                          , overwrite = overwrite $
                          , tag = tag
  endif
     
end

dir = '/scratch_beegfs/mats/NEW/2025-08-18/CRISP2/'
cd, dir
a = crisp2red(/dev, /no)

infile = 'cubes_nb/nb_6173_2025-08-18T08:52:59_08:52:59=1-4_stokes_corrected_im.fits'

a -> fitscube_mask, infile, ituning = 6, /overwrite ;, /no_qhull

end
