; docformat = 'rst'

pro fitscube_crosstalk_event, event

  ;; Here all events are processed. As a first step, retrieve the
  ;; pointer from the widget.
  widget_control, event.top, get_uvalue = info

  ;;Check which widget generated the event
  case event.id of

    (*info).wSlider : begin
      ;; Process the info from the widget_slider.

      ;; {WIDGET_SLIDER, ID:0L, TOP:0L, HANDLER:0L, VALUE:0L, DRAG:0}
      ;; What you only really need now is "VALUE"
      ;; This is the position of the slider. 
      (*info).cutoff = event.value

      ;; Change the image based on the weight and the selected cutoff
      indx = where(finite((*info).image))
      im2 = (*info).image + !values.f_nan
      indx2 = indx[where((*info).weight gt (*info).cutoff)]
      im2[indx2] = (*info).image[indx2]

      indx_data = where(finite(im2), Ndata, complement = indx_missing, ncomplement = Nmissing)
      mn = median(im2[indx_data])-2.5*robust_sigma(im2[indx_data])
      mx = median(im2[indx_data])+2.5*robust_sigma(im2[indx_data])

      imshow = bytscl(im2, mn, mx)
      rg = imshow
      b = imshow
      if Nmissing gt 0 then begin
        ;; Yellow for missing data
        rg[indx_missing] = 255
        b[indx_missing] = 0
      endif
      imshow = [[[rg]],[[rg]],[[b]]]
                                ;im2 = bytscl(im2, mn, mx)
      
      ;;  im2 = cgGmaScl((*info).image, gamma=(*info).gamma)
      ;; Display it
      ;;  (*info).oImg -> setData, im2
      (*info).oImg -> setData, imshow
    end

    (*info).wClose : begin
      ;; Destroy the widget.
      widget_control, (*info).tlb, /destroy
    end
    
    else :                      ; Nothing...
    
  endcase

  
end

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
;     force : in, optional, tyope=boolean
;   
;       Do not care if the correction is done already.
;
;     mag_mask : in, out,  optional, type=array
;
;       A mask that deselects magnetic signal. Should have the same
;       dimensions as a frame. Use keyword primarily to preserve a
;       mask from one call to the next so user only has to to use the
;       GUI once.
; 
;     margin : in, optional, type=integer, default=15
; 
;       Remove this many pixels from all edges when making the
;       magnetic mask. 
; 
;     tuning_selection : in, out, optional, type="integer or array"
;
;       The index or indices in the tuning dimension to use for
;       calculating the correction. Should correspond to continuum (or
;       as close to continuum as possible), where the polarimetric
;       signal is minimal. Negative indices are allowed. Keyword can
;       be used to preserve a mask from one call to the next so user
;       only has to to use the GUI once.
;
; :History:
; 
;   2019-04-01 : MGL. First version.
; 
;   2019-04-04 : MGL. New keywords nostatistics and tuning_selection. 
; 
;   2019-04-05 : MGL. Use red_fitscube_open and red_fitscube_close. 
;
;   2020-03-23 : OA. Changed calls of red_fitscube_getframe &
;                *_addframe with fileassoc instead of filename. Added
;                keyword 'fitscube_info' in those calls.
; 
;   2020-10-28 : MGL. Remove statistics calculations.
; 
;   2021-04-15 : MGL. Improved magnetic mask construction. New keyword
;                margin.
; 
;   2021-05-05 : MGL. New GUI for constructing the magnetic mask.
;
;-
pro red::fitscube_crosstalk, filename  $
                             , flip = flip $
                             , force = force $
                             , mag_mask = mag_mask $
                             , margin = margin $
                             , test = test $
                             , tuning_selection = tuning_selection 

  inam = red_subprogram(/low, calling = inam1)

  red_fitscube_open, filename, fileassoc, fitscube_info, /update

  if n_elements(margin) eq 0 then margin = 50
  
  hdr = fitscube_info.header
  
  ;; Information about processing steps in the formation of the input
  ;; file. 
  prprocs = fxpar(hdr, 'PRPROC*')
  prparas = fxpar(hdr, 'PRPARA*')

  ;; Cube dimensions
  Nx      = fitscube_info.dimensions[0]
  Ny      = fitscube_info.dimensions[1]
  Ntuning = fitscube_info.dimensions[2]
  Nstokes = fitscube_info.dimensions[3]
  Nscans  = fitscube_info.dimensions[4]

  ;; Check that it is in fact a polarimetric cube
  if Nstokes eq 1 then begin
    print
    print, inam + ' : This is not a polarimetric cube:'
    print, filename
    red_fitscube_close, fileassoc, fitscube_info
    return
  endif

  if ~keyword_set(force) then begin
    ;; Check that it is not already corrected for crosstalk.
    pos = where(strmatch(prprocs, inam), Nmatch)
    if Nmatch gt 0 then begin
      print
      print, inam + ' : This file is already corrected for crosstalk:'
      print, filename
      print, inam + " : Use /force to do it again."
      red_fitscube_close, fileassoc, fitscube_info
      return
    endif
  endif
  
  ;; Get the wavelength coordinate
  red_fitscube_getwcs, filename, coordinates = coordinates
  wav = coordinates[*,0].wave[0,0]

  ;; Get scan numbers
  scannum = red_fitsgetkeyword(filename, 'SCANNUM', variable_values = variable_values)
  if n_elements(variable_values) gt 0 then begin
    scannumbers = reform(variable_values.values)
    if Nscans ne n_elements(scannumbers) then begin
      print, 'n_elements(scannumbers) ne Nscans (i.e. NAXIS5).'
      print,'Something is wrong with the cube header. Please, check.'
      stop
    endif
  endif else begin
    scannumbers = [scannum]
    if Nscans ne n_elements(scannumbers) then begin
      print, 'n_elements(scannumbers) ne Nscans (i.e. NAXIS5).'
      print,'Something is wrong with the cube header. Please, check.'
      stop
    endif
  endelse
  
  ;; Get medians of the I component of the first scan, to be used for
  ;; selecting the wavelength points.
  medi = fltarr(Ntuning)
  for ituning = 0, Ntuning-1 do begin
    red_fitscube_getframe, fileassoc, frame, istokes = 0, iscan = 0, ituning = ituning,fitscube_info = fitscube_info
    medi[ituning] = median(frame)
  endfor                        ; ituning

  if n_elements(tuning_selection) gt 0 then begin
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
    tuning_selection = ppc
    wait, .1                    ; Make sure the plot has time to color the selected spots
  endelse

  
  if n_elements(mag_mask) eq 0 then begin

    ;; Make a magnetic mask based on Stokes V images of the selected
    ;; tunings, for all scans.
    lows  = fltarr(Nx, Ny)
    highs = fltarr(Nx, Ny)
    Vcube = fltarr(Nx, Ny, Nscans, n_elements(ppc))
    for iscan = 0, Nscans-1 do begin
      red_progressbar, iscan, Nscans, 'Reading and processing image data'
      for ippc = 0, n_elements(ppc)-1 do begin
        red_fitscube_getframe, fileassoc, frame, istokes = 3, iscan = iscan, ituning = ppc[ippc] $
                               , fitscube_info = fitscube_info
        red_missing, frame, missing_type_wanted = 'nan', /inplace
        
        frame_mask = finite(frame)
        frame_mask = erode(frame_mask,replicate(1,2*margin,2*margin))
        frame_indx = where(~frame_mask, Nmissing)
        if Nmissing gt 0 then frame[frame_indx] = !values.f_nan
        
        Vcube[0, 0, iscan, ippc] = frame
      endfor                    ; ippc
    endfor                      ; iscan
    Vcube = reform(vcube,Nx,Ny,Nscans*n_elements(ppc))
    
;    lows  = min(Vcube,dim=3, /nan)
;    highs = max(Vcube,dim=3, /nan)


    vdims = size(vcube, /dim)
    if n_elements(vdims) eq 2 || vdims[2] eq 1 then begin
      im = reform(vcube)
      ma = abs(im)
    endif else begin
      im = mean(vcube, dim=3, /nan)
      ma = max(abs(vcube), dim=3, /nan, ml)    
    endelse
;    indx = where(finite(im))

    window, 10, xs = Nx, ys = Ny
    indx_data = where(finite(im), Ndata, complement = indx_missing, ncomplement = Nmissing)
    mn = median(im[indx_data])-2.5*robust_sigma(im[indx_data])
    mx = median(im[indx_data])+2.5*robust_sigma(im[indx_data])
    cgimage, red_histo_opt(im) $
             , missing_color = 'black' $             
             ;, missing_color = 'yellow' $             
             , missing_index = 0 $
             , stretch = 1  
    
    
;    sgp = long(abs(highs) gt abs(lows))
;    sgm = long(abs(highs) lt abs(lows))
;    sg = sgp - sgm
;
;    
;    stop

;    rs = robust_sigma(ma[indx], goodvec = q)
;    im3 = im + !values.f_nan
;    im3[indx[q]] = im[indx[q]]

    mn = biweight_mean(ma[indx_data], sigma, w)

;    indx1 = where(finite(vcube[*, *, 0]))
;    mn1 = biweight_mean(abs((vcube[*, *, 0])[indx1]), sigma1, w1)
    
;    ma1 = ma + !values.f_nan
;    ma1[indx[q]] = ma[indx[q]]

;    ma2 = ma + !values.f_nan
;    indx2 = indx[where(w/max(w) gt .8)]
;    ma2[indx2] = ma[indx2]

    cutoff = 0.85
    
 
    print
    print, inam + ' : Please inspect the image in windows 10 and 11.'
    print
    print, 'Window 10:'
    red_strflow, ['This image is the average image from the selected tunings and all' $
                  , 'the scans. You can use this image in XROI and deselect regions that' $
                  , 'show magnetic activity. Note that all frames have been cropped by' $
                  , 'a margin of '+strtrim(margin, 2)+' pixels before averaging. You' $
                  , 'can change this with the "margin" keyword.']
    print
    print, 'Window 11:'
    red_strflow, ['This image is the same as the one in Window 10 but an attempt has been' $
                  , 'made to identify magnetic activity and mask those pixels out. The' $
                  , 'masked pixels are yellow. You can use this image in XROI and deselect' $
                  , 'additional regions that show magnetic activity. Or you can accept its' $
                  , 'masking as is.']
    print
    red_strflow, ['You can also choose to change the cutoff value of the automatic mask, either in a GUI or with the keyboard.' $
                  , 'A larger value will mask more pixels, valid interval is [0,1].' $
                  , 'You will then be shown the new automatic mask and can choose again.']
    print
    red_strflow, ['When is the mask acceptable? When it masks enough magnetic' $
                  , 'activity that there are many more non-magnetic pixels than magnetic' $
                  , 'ones. But it must not mask noise peaks and granulation pattern from' $
                  , 'the non-magnetic areas.']
    print

    repeat begin

      im2 = im + !values.f_nan
      indx2 = indx_data[where(w/max(w) gt cutoff)]
      im2[indx2] = im[indx2]
      
      
      window, 11, xs = Nx, ys = Ny
      indx_data = where(finite(im2), Ndata, complement = indx_missing, ncomplement = Nmissing)
      mn = median(im2[indx_data])-2.5*robust_sigma(im2[indx_data])
      mx = median(im2[indx_data])+2.5*robust_sigma(im2[indx_data])
      cgimage, bytscl(im2, mn, mx)$
               , missing_color = 'yellow' $
               , missing_index = 0 $
               , stretch = 1  

      print, 'Current cutoff: ', cutoff
      tmp = red_select_subset( ['XROI with image from window 10' $
                                , 'XROI with image from window 11' $
                                , 'Change cutoff in GUI' $
                                , 'Change cutoff with keyboard' $
                                , 'Accept mask in window 11' $
                               ] $
                               , count = count $
                               , default = 4 $
                               , indx = sel $
                               , maxcount = 1  $
                               , qstring = 'What is your preference?' $
                             )
      print

      case sel of

        0 : mag_mask = red_select_area(im,  /xroi) ; Average V image

        1 : mag_mask = red_select_area(im2, /xroi) ; Average V image with automatic mask

        2 : begin               ; GUI to modify automatic mask
          ;; Create a top level widget
          tlb = widget_base(title='Mask cutoff selector',/column)
          ;; Add a window widget to display the image
          screen_dims = get_screen_size()
          x_scroll_size = screen_dims[0]-200
          wWindow = widget_window(tlb, xsize=Nx, ysize=Ny $
                                  , x_scroll_size=x_scroll_size, y_scroll_size=screen_dims[0]-50 $
                                 )
          ;; Add slider widget. cw_fslider is like window_slider but handles floats.
          wSlider = cw_fslider(tlb, title='Mask cutoff' $
                               , minimum=0.00, value=cutoff, maximum=1.00 $
                               , format = '(f5.3)' $
                               , xSize=x_scroll_size <Nx, /drag)
          ;; Add a button to close the widget
          wClose = widget_button(tlb, value='Close')
          ;; Realize the widget (make them, but do nothing more yet!)
          widget_control, tlb, /realize
          ;; Get the windog object reference
          widget_control, wWindow, get_value=oWin
          ;; Put an image in the window
          imshow = bytscl(im2, mn, mx)
          rg = imshow
          b = imshow
          if Nmissing gt 0 then begin
            ;; Yellow for missing data
            rg[indx_missing] = 255
            b[indx_missing] = 0
          endif
          imshow = [[[rg]],[[rg]],[[b]]]
          oImg = image(imshow, margin=0, current=oWin)
          ;;gotta store all the info somewhere...
          info = ptr_new({tlb         : tlb      $
                          , wClose    : wClose   $
                          , wWindow   : wWindow  $
                          , wSlider   : wSlider  $
                          , image     : im       $
                          , weight    : w/max(w) $
                          , cutoff    : cutoff   $
                          , oImg      : oImg     $
                         })
          ;; And be able to retrieve it
          widget_control, tlb, set_uvalue=info
          ;; Now activate the widget!
          xmanager, 'red__fitscube_crosstalk', tlb $
;                      , cleanup = 'mats_sliderWidget_cleanup' $
                    , event_handler = 'fitscube_crosstalk_event'
          
          cutoff = (*info).cutoff
          print, 'Cutoff=', cutoff
        end 

        3 : begin
          cutoff_answer = ''
          read, 'Enter new cutoff (current value is '+string(cutoff, format = '(f4.2)')+') : ', cutoff_answer
          print
          if cutoff_answer ne '' then cutoff = float(cutoff_answer) < 1
        end
        
        4: mag_mask = finite(im2) ; Accept current mask

        else : stop

      endcase
      
    endrep until sel ne 2 and sel ne 3
    
;    tighttv, ma, 0
;    tighttv, ma2, 1
;    tighttv, im, 2
;    tighttv, im2, 3
  
;    a = red_histo_gaussfit(vcube, fw = .05, /nan)
;    
;    loindx = where(abs(lows) gt abs(highs), complement = hiindx)
;    im = fltarr(Nx, Ny)
;    im[loindx] = lows[loindx]
;    im[hiindx] = highs[hiindx]

;    stop

;    if n_elements(ppc) gt 1 then begin
;      im = 0.
;      for i = 0, n_elements(ppc)-1 do begin
;        red_fitscube_getframe, fileassoc, frame, istokes = 3, iscan = 0, ituning = ppc[i],fitscube_info = fitscube_info
;        im += abs(frame)
;      endfor
;    endif else begin
;      red_fitscube_getframe, fileassoc, im, istokes = 3, iscan = 0, ituning = ppc[0], fitscube_info = fitscube_info
;    endelse

;    print, 'Deselect areas with magnetic structures and/or artifacts.'
;    print, 'Yellow areas are already deselected'
;    print, 'End with File-->Quit.'
;;;mag_mask = red_select_area(red_histo_opt(im,2.e-3), /noedge, /xroi)
;    ;;   mag_mask = red_select_area(red_histo_opt(im,2.e-3), /xroi)
;    mag_mask = red_select_area(im2, /xroi)
;
  
  endif

  if keyword_set(test) then return
  
  ;; Get name of WB cube from the NB cube-making parameters, used to
  ;; make a mask that removes rotational padding.
  pos = where(strmatch(prprocs, '*make_nb_cube'), Nmatch)
  makemask = Nmatch ne 0
  
  if makemask then begin
    make_nb_cube_paras = prparas[pos[0]]
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

    ;; Dimensions of non-rotated images .
    Nxx = wcX01Y01[1] - wcX01Y01[0] + 1 
    Nyy = wcX01Y01[3] - wcX01Y01[2] + 1 
    
  endif

  for iscan = 0, Nscans-1 do begin

    if makemask then begin
      ;; Construct a mask for the padding
      pad_mask = make_array(Nxx, Nyy, /float, value = 1.) 
      pad_mask = red_rotation(pad_mask, ang[iscan], wcshift[0,iscan], wcshift[1,iscan] $
                              , background = 0, full = wcFF)
      pindx = where(pad_mask le 0.99) ; Pixels that are padding
      
      ;; Combine with the magnetic mask
      this_mask = mag_mask * pad_mask
      mindx = where(this_mask)
    endif else begin
      mindx = where(mag_mask)
    endelse 
      
    ;; Calculate the amount of crosstalk
    numerator   = dblarr(Nstokes)
    denominator = 0d
    for i = 0, n_elements(ppc)-1 do begin
      red_fitscube_getframe, fileassoc, im0, istokes = 0, iscan = iscan, ituning = ppc[i] $
                             , fitscube_info = fitscube_info ; Stokes I

      ;;denominator += median(im0[where(this_mask)] *
      ;;im0[where(this_mask)], /double)
      
      ;; Find the centroid by fitting a Gaussian
      a = red_histo_gaussfit(im0[mindx] * im0[mindx], FWlevel = 0.25)
      denominator += a[1]
      
      for istokes=1, Nstokes-1 do begin
        red_fitscube_getframe, fileassoc, im, istokes = istokes, iscan = iscan, ituning = ppc[i] $
                               , fitscube_info = fitscube_info

        ;;numerator[istokes] += median(im0[where(this_mask)] * im[where(this_mask)], /double)
        a = red_histo_gaussfit(im0[mindx] * im[mindx], FWlevel = 0.25)
        numerator[istokes] += a[1]
      endfor                    ; istokes
    endfor
    crt = numerator/denominator 
    
    print, 'Scan '+strtrim(scannumbers[iscan], 2)+' : crosstalk from I -> Q,U,V =' $
           , crt[1], ', ', crt[2], ', ', crt[3], format='(A,F8.5,A,F8.5,A,F8.5)'

    ;; Apply the crosstalk correction
    for ituning = 0, Ntuning-1 do begin
      ;; Read Stokes I
      red_fitscube_getframe, fileassoc, im0, istokes = 0, iscan = iscan, ituning = ituning $
                             , fitscube_info = fitscube_info 
      if makemask then begin
        im0[pindx] = median(im0[pindx]) ; padding = median(padding)
        red_fitscube_addframe, fileassoc, im0, istokes = 0, iscan = iscan, ituning = ituning $
                             , fitscube_info = fitscube_info ; Write with updated padding
      endif
      ;; Read and apply to Q, U, V
      for istokes=1, Nstokes-1 do begin
        red_fitscube_getframe, fileassoc, im, istokes = istokes, iscan = iscan, ituning = ituning $
                               , fitscube_info = fitscube_info
        im -= float(crt[istokes] * im0)
        if makemask then im[pindx] = median(im[pindx]) ; padding = median(padding)
        ;; The padding manipulation above should perhaps be done by red_missing.
        red_fitscube_addframe, fileassoc, im, istokes = istokes, iscan = iscan, ituning = ituning $
                               , fitscube_info = fitscube_info
      endfor                    ; istokes
    endfor                      ; ituning
  endfor                        ; iscan

  ;; Add info about this step
  red_make_prpara, prpara, margin
  red_make_prpara, prpara, force
  red_make_prpara, prpara, tuning_selection
  self -> headerinfo_addstep, hdr $
                              , prstep = 'STOKES-CROSSTALK-CORRECTION' $
                              , prpara = prpara $
                              , prproc = inam

  ;; Close the file and write the updated header
  red_fitscube_close, fileassoc, fitscube_info, newheader = hdr
  
  if keyword_set(flip) then begin
    red_fitscube_flip, filename, flipfile = flipfile, /overwrite
  endif
  
end

a = crispred(/dev)

cd, '/scratch/mats/2016.09.19/CRISP-aftersummer/'
filename = 'cubes_nb/nb_6302_2016-09-19T09:30:20_scans=0-43_stokes_corrected_im.fits'
filename = 'cubes_nb/nb_6302_2016-09-19T09:30:20_scans=2-20_stokes_corrected_im.fits'
filename = 'cubes_nb/nb_6302_2016-09-19T09:30:20_scans=2-8_stokes_corrected_im.fits'
force = !true
a -> fitscube_crosstalk, filename, /test $
                         , force = force $
                         , mag_mask = mag_mask $
                         , tuning_selection = tuning_selection



end
