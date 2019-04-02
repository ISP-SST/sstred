; docformat = 'rst'

;+
; Correct a polarimetric fitscube for crosstalk from I to Q, U, and V. 
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
;     filename : in, type=string
; 
;       The name of the file containing the polarimetric fitscube.
; 
; 
; :Keywords:
; 
;     flip : in, optional, type=boolean
;   
;       Produce a flipped version if this keyword is set. 
; 
;    force : in, optional, tyope=boolean
;   
;       Do not care if the correction is done already.
; 
; :History:
; 
;   2019-04-01 : MGL. First version.
; 
;-

pro red::fitscube_crosstalk, filename  $
                             , flip = flip $
                             , force = force $
                             , tuning_selection = tuning_selection

  inam = red_subprogram(/low, calling = inam1)

  hdr = headfits(filename)

  ;; Information about processing steps in the formation of the input
  ;; file. 
  prprocs = fxpar(hdr, 'PRPROC*')
  prparas = fxpar(hdr, 'PRPARA*')

  ;; Cube dimensions
  naxis = fxpar(hdr, 'NAXIS*')
  Nx      = naxis[0]
  Ny      = naxis[1]
  Ntuning = naxis[2]
  Nstokes = naxis[3]
  Nscans  = naxis[4]

  ;; Check that it is in fact a polarimetric cube
  if Nstokes eq 1 then begin
    print, inam + ' : This is not a polarimetric cube:'
    print, filename
    return
  endif

  if ~keyword_set(force) then begin
    ;; Check that it is not already corrected for crosstalk.
    pos = where(strmatch(prprocs, inam), Nmatch)
    if Nmatch gt 0 then begin
      print, inam + ' : This file is already corrected for crosstalk:'
      print, filename
      print, inam + " : Don't do it again."
      return
    endif
  endif
  
  ;; Get the wavelength coordinate
  red_fitscube_getwcs, filename, coordinates = coordinates
  wav = coordinates[*,0].wave[0,0]

  ;; Get scan numbers
  tmp = red_fitsgetkeyword(filename, 'SCANNUM', variable_values = variable_values)
  scannumbers = reform(variable_values.values)
  
  ;; Get medians of the I component of the first scan, to be used for
  ;; selecting the wavelength points.
  medi = fltarr(Ntuning)
  for ituning = 0, Ntuning-1 do begin
    red_fitscube_getframe, filename, frame, istokes = 0, iscan = 0, ituning = ituning
    medi[ituning] = median(frame)
  endfor                        ; ituning

  if keyword_set(tuning_selection) then begin
    ppc = tuning_selection
    ;; Translate negative indices to positive ones
    negindx = where(ppc lt 0, Nwhere)
    if Nwhere gt 0 then begin
      ppc[negindx] = Ntuning + ppc[negindx]
    endif
  endif else begin
    ;; Choose spectral points to use. We want as little signal as
    ;; possible so continuum points are good. For wide lines we
    ;; might not have them so pick end points if similar intensity,
    ;; or just one endpoint if one has significantly higher
    ;; intensity than the other.
    print, 'Select spectral points to calculate cross-talk from. Select with left mouse, end with right mouse.'
    ppc = red_select_spoints(wav, medi)
  endelse

  
  if n_elements(ppc) gt 1 then begin
    im = 0.
    for i = 0, n_elements(ppc)-1 do begin
      red_fitscube_getframe, filename, frame, istokes = 3, iscan = 0, ituning = ppc[i]
      im += abs(frame)
    endfor
  endif else begin
    red_fitscube_getframe, filename, im, istokes = 3, iscan = 0, ituning = ppc[0]
  endelse

  print, 'Deselect areas with magnetic structures and/or artifacts. End with File-->Quit.'
  ;;mask = red_select_area(red_histo_opt(im,2.e-3), /noedge, /xroi)
  mask = red_select_area(red_histo_opt(im,2.e-3), /xroi)

  
  ;; Get name of WB cube from the NB cube-making parameters
  pos = where(strmatch(prprocs, '*make_nb_cube'), Nmatch)
  if Nmatch eq 0 then stop
  make_nb_cube_paras = prparas[pos]
  wcfile = (json_parse(make_nb_cube_paras, /tostruct)).wcfile

  ;; Get the rotation and alignment parameters from the wideband cube
  if ~file_test(wcfile) then stop
  fxbopen, bunit, wcfile, 'MWCINFO', bbhdr
  fxbreadm, bunit, row = 1 $
            , ['ANG', 'CROP', 'FF', 'GRID', 'ND', 'SHIFT', 'TMEAN', 'X01Y01'] $
            ,   ANG, wcCROP, wcFF, wcGRID, wcND, wcSHIFT, wcTMEAN, wcX01Y01
  ;; Note that the strarr wfiles cannot be read by fxbreadm! Put it in
  ;; wbgfiles (WideBand Global).
  fxbread, bunit, wbgfiles, 'WFILES', 1
  fxbclose, bunit

  ;; Get dimensions of non-rotated images from the momfbd-restored WB
  ;; images.
  whdr = red_readhead(wbgfiles[0])
  Nxx = fxpar(whdr, 'NAXIS1')
  Nyy = fxpar(whdr, 'NAXIS2')

  for iscan = 0, Nscans-1 do begin

    ;; Construct a mask for the padding
    pad_mask = make_array(Nxx, Nyy, /float, value = 1.) 
    pad_mask = red_rotation(pad_mask, ang[iscan], wcshift[0,iscan], wcshift[1,iscan], background = 0, full = wcFF)
    pindx = where(pad_mask le 0.99) ; Pixels that are padding

    ;; Include the padding mask just in case it rotates into the
    ;; selected mask.
    this_mask = mask * pad_mask
    mindx = where(this_mask)
        
    ;;crt = red_get_ctalk(d, idx=ppc, mask=pixmask)
    crt = dblarr(Nstokes)
    numerator   = dblarr(Nstokes)
    denominator = 0d
    for i = 0, n_elements(ppc)-1 do begin
      red_fitscube_getframe, filename, im0, istokes = 0, iscan = iscan, ituning = ppc[i] ; Stokes I

      ;;denominator += median(im0[where(this_mask)] *
      ;;im0[where(this_mask)], /double)
      
      ;; Find the centroid by fitting a Gaussian
      a = red_histo_gaussfit(im0[mindx] * im0[mindx], FWlevel = 0.25)
      denominator += a[1]
      
      for istokes=1, Nstokes-1 do begin
        red_fitscube_getframe, filename, im, istokes = istokes, iscan = iscan, ituning = ppc[i]
        ;;numerator[istokes] += median(im0[where(this_mask)] * im[where(this_mask)], /double)
        a = red_histo_gaussfit(im0[mindx] * im[mindx], FWlevel = 0.25)
        numerator[istokes] += a[1]
      endfor                    ; istokes
    endfor
    crt = numerator/denominator 
    
    print, 'Scan '+strtrim(scannumbers[iscan], 2)+' : crosstalk from I -> Q,U,V =' $
           , crt[1], ', ', crt[2], ', ', crt[3], format='(A,F8.5,A,F8.5,A,F8.5)'
    
    for ituning = 0, Ntuning-1 do begin
      red_fitscube_getframe, filename, im0, istokes = 0, iscan = 0, ituning = ituning ; Stokes I
      im0[pindx] = median(im0[pindx])                                                 ; Set the padding to median
      red_fitscube_addframe, filename, im0, istokes = 0, iscan = 0, ituning = ituning    ; Write with updated padding
      for istokes=1, Nstokes-1 do begin
        ;;d[*,*,tt,ww] -= crt[tt]*d[*,*,0,ww]
        red_fitscube_getframe, filename, im, istokes = istokes, iscan = iscan, ituning = ituning
        im -= float(crt[istokes] * im0)
        im[pindx] = median(im[pindx]) 
        red_fitscube_addframe, filename, im, istokes = istokes, iscan = iscan, ituning = ituning
      endfor                    ; istokes
    endfor                      ; ituning
  endfor                        ; iscan

  
  ;; Add info about this step
  self -> headerinfo_addstep, hdr $
                              , prstep = 'Crosstalk correction' $
                              , prpara = prpara $
                              , prproc = inam

  ;; Write the updated header (may be time consuming!)
  modfits, filename, 0, hdr

  if keyword_set(flip) then begin
    self -> fitscube_flip, filename, flipfile = flipfile
  endif
  
end
