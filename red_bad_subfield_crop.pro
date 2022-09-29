; docformat = 'rst'

;+
; Find cropping based on bad subfields in momfbd output files, either
; automatically or interactively - or both.
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
;    files : in, type=strarr
; 
;      The momfbd (wideband) files to use to detect the bad subfields.
; 
;    crop : in, out, type=fltarr(4) 
;
;      See documentation in make_wb_cube method. Note, this cropping
;      is by default specified in the momfbd-output orientation. Use
;      the direction keyword to change this - on output only.
; 
; :Keywords:
;
;    autocrop : in, optional, type=booean
;
;      Try to determine the largest FOV that avoids any bad momfbd
;      subfields along the edges. If this keyword is set, the input
;      value of the crop keyword is ignored and is set to the
;      auto-detected crop parameters.
;
;    direction : in, optional, type=integer, default=0
;
;      The relative orientation of reference cameras of different
;      instruments. The crop parameter will be re-arranged if this
;      keyword is used.
;
;    interactive : in, optional, type=boolean
;
;      Set this keyword to define the FOV by use of the XROI GUI. If
;      autocrop is set, then use the so defined FOV as an
;      initialization of the FOV in the GUI. Otherwise use the crop
;      keyword (or its default).
; 
; 
; :History:
; 
;    2018-01-12 : MGL. First version, based on code taken from
;                 chromis::make_wb_cube.
; 
;    2020-12-23 : MGL. New keyword direction.
; 
;-
function red_bad_subfield_crop, files, crop $
                                , autocrop = autocrop $
                                , direction =  direction $
                                , interactive = interactive

  Nfiles = n_elements(files)

  if n_elements(direction) eq 0 then direction = 0
  
  if keyword_set(autocrop) then begin
    
    dispim = 0.0

    for ifile = 0L, Nfiles -1 do begin

      red_progressbar, ifile, Nfiles, 'Analyze momfbd subfields for autocrop', /predict

      mr = momfbd_read(files[ifile], /img)

      if keyword_set(interactive) then dispim += red_mozaic(mr, /crop)

      if ifile eq 0 then begin
        subf_dims = size(mr.patch.img, /dim)
        Ssubf_y = subf_dims[0]
        Ssubf_x = subf_dims[1]
        Nsubf_y = subf_dims[2]
        Nsubf_x = subf_dims[3]
        subf_mean   = fltarr(Nsubf_x, Nsubf_y, Nfiles)
        subf_median = fltarr(Nsubf_x, Nsubf_y, Nfiles)
        subf_stddev = fltarr(Nsubf_x, Nsubf_y, Nfiles)
      endif

      for isubf_x = 0, Nsubf_x-1 do for isubf_y = 0, Nsubf_y-1 do begin
        subf_mean[isubf_x, isubf_y, ifile]   = mean(mr.patch[isubf_y, isubf_x].img)
        subf_median[isubf_x, isubf_y, ifile] = median(mr.patch[isubf_y, isubf_x].img)
        subf_stddev[isubf_x, isubf_y, ifile] = stddev(mr.patch[isubf_y, isubf_x].img)
      endfor                    ; isubf_x, isubf_y
    endfor                      ; ifile

    if size(subf_stddev, /n_dim) eq 2 then begin
      subf_detect = subf_stddev*subf_median^2
    endif else begin
      subf_detect = total(subf_stddev*subf_median^2, 3)
    endelse
    
    ;; Try removing subfields in X first
    subf_detect_x = max(subf_detect/median(subf_detect),dim=2)
    i0 = -1
    i1 = Nsubf_x
    repeat i0++ until subf_detect_x[i0] lt 2 or i0 eq Nsubf_x-1
    repeat i1-- until subf_detect_x[i1] lt 2 or i1 eq 0
    if i0 ge i1 then begin
      ;; Failed autodetection, use entire FOV
      i0 = 0
      i1 = Nsubf_x-1      
    endif
    
    ;; Then Y
    subf_detect_y = max(subf_detect[i0:i1, *]/median(subf_detect[i0:i1]),dim=1)
    j0 = -1
    j1 = Nsubf_y
    repeat j0++ until subf_detect_y[j0] lt 2 or j0 eq Nsubf_y-1
    repeat j1-- until subf_detect_y[j1] lt 2 or j1 eq 0
    if j0 ge j1 then begin
      ;; Failed autodetection, use entire FOV
      j0 = 0
      j1 = Nsubf_y-1    
    endif

    hdr = red_readhead(files[0])
    im_dim = fxpar(hdr, 'NAXIS*')

;    im = red_mozaic(mr, /clip)
;    im_dim = size(im, /dim)

    overlapfacs = [Ssubf_x*Nsubf_x/float(im_dim[0])*[1, 1], Ssubf_y*Nsubf_y/float(im_dim[1])*[1, 1]] * 0.8

    ;; If user selected autocrop, then do use the detected cropping.
    ;; If user selected interactive, then use the detected cropping
    ;; only if crop keyword was not provided.
    
    roi_name = 'Autocrop'
    crop = ceil( [[i0, (Nsubf_x-1-i1)]*Ssubf_x, [j0, (Nsubf_y-1-j1)]*Ssubf_y] / overlapfacs )
    
  endif else if keyword_set(interactive) then begin

    dispim = 0.0

    for ifile = 0L, Nfiles -1 do begin
      red_progressbar, ifile, Nfiles, 'Read momfbd output for interactive crop', /predict
      mr = momfbd_read(files[ifile], /img)
      dispim += red_mozaic(mr, /crop)
    endfor                      ; ifile

  endif
  
  hdr = red_readhead(files[0])
  im_dim = fxpar(hdr, 'NAXIS*')

  ;; Default cropping
  case n_elements(crop) of
    1 : begin                   ; Same crop from all sides
      roi_name = 'From crop keyword'
      crop = replicate(crop, 4)
    end
    2 : begin                   ; One crop in X, another in Y
      roi_name = 'From crop keyword'
      crop = [crop[0], crop[0], crop[1], crop[1]]
    end
    4 : begin                   ; Leave it alone
      if n_elements(roi_name) eq 0 then roi_name = 'From crop keyword'
    end
    else : begin
      roi_name = 'Default'
      crop = [0, 0, 0, 0]
    end
  endcase

  if crop[0] lt 0 || crop[2] lt 0 $
     || crop[1] ge im_dim[0] || crop[3] ge im_dim[1] then begin
    roi_name = 'Default'
    crop = [0, 0, 0, 0]
  endif
  
  x0 = crop[0]
  x1 = im_dim[0]-1 - crop[1]
  y0 = crop[2]
  y1 = im_dim[1]-1 - crop[3]
  
  if keyword_set(interactive) then begin

    ;; Use XROI GUI to select a rectangular area. 

    red_missing, dispim, missing_type_wanted='nan', /inplace
    
    ;; Initialize the FOV
    X_in = [x0, x1, x1, x0]
    Y_in = [y0, y0, y1, y1]
    roi_in = OBJ_NEW('IDLgrROI', X_in, Y_in)
    roi_in -> setproperty, name = roi_name

    print
    print, 'Use the XROI GUI to either modify an initial ROI/FOV or define a new'
    print, 'one from scratch. Look for bad (hightened contrast) subfields.'
    print, 'Select Quit in the File menu. The last ROI is used.'
    print
    
    ;; Fire up the XROI GUI.
    xroi, bytscl(red_histo_opt(dispim)), regions_in = [roi_in], regions_out = roi, /block $
          , tools = ['Translate-Scale', 'Rectangle'] $
          , title = 'Modify or define FOV based on summed image'

    roi[-1] -> getproperty, roi_xrange = roi_x
    roi[-1] -> getproperty, roi_yrange = roi_y

    obj_destroy, roi_in
    obj_destroy, roi

    ;; We don't want out-of-bounds coordinates here
    x0 = round(roi_x[0]) >0
    y0 = round(roi_y[0]) >0
    x1 = round(roi_x[1]) <(im_dim[0]-1)
    y1 = round(roi_y[1]) <(im_dim[1]-1)


    crop = [x0, im_dim[0]-1-x1, y0, im_dim[1]-1-y1]
    
  endif
  
  case direction of
    0 :                           ;  X,  Y (no change)
    1 : crop = crop[[3, 2, 0, 1]] ; -Y,  X
    2 : crop = crop[[1, 0, 3, 2]] ; -X, -Y
    3 : crop = crop[[2, 3, 1, 0]] ;  Y, -X
    4 : crop = crop[[2, 3, 0, 1]] ;  Y,  X
    5 : crop = crop[[1, 0, 2, 3]] ; -X,  Y
    6 : crop = crop[[3, 2, 1, 0]] ; -Y, -X
    7 : crop = crop[[0, 1, 3, 2]] ;  X, -Y
    else : stop
  endcase

  if max(direction eq [1, 3, 4, 6]) eq 1 then begin
    ;; X and Y switched
    y0 = crop[2]
    y1 = im_dim[0]-1 - crop[3]
    x0 = crop[0]
    x1 = im_dim[1]-1 - crop[1]
  endif else begin
    x0 = crop[0]
    x1 = im_dim[0]-1 - crop[1]
    y0 = crop[2]
    y1 = im_dim[1]-1 - crop[3]
  endelse

  return, [x0, x1, y0, y1]

end







