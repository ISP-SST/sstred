; docformat = 'rst'

;+
; Add keywords required by GaVO DaCHs to a fitscube header.
;
; Add RESOLVPW, CSPMINn/MAXn, CBBMINn/MAXn and CROTA.
; Also OBS_MODE and POINT_ID (if missing)
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
;    filename : in, type=string
; 
;      The path of the fitscube file. 
; 
; 
; :Keywords:
; 
;    anchor : in, out, optional, type=string
; 
;      See red_fitsaddkeyword.pro
; 
;    test_checksum : in, optional, type=boolean
;   
;      Test the CHECKSUM keyword before adding new keywords and update
;      it afterwards.
; 
;    crota : in, optional, type=float
; 
;      Counterclockwise rotation of the FOV in degrees. If not given,
;      don't add the CROTA keyword.
; 
;    hdr : in, out, optional, type=strarr
; 
;      A header to be updated. If present, don't read the header from
;      the file and don't write it either. Also implies checksum=0.
; 
; :History:
; 
;    2026-06-24 : MGL. First version.
; 
;-
pro red_fitscube_addgavokeywords, filename $
                                  , anchor = anchor $
                                  , checksum = checksum $
                                  , crota = crota $
                                  , hdr = hdr

  if n_elements(anchor) eq 0 then anchor = 'FILENAME'
  
  hdr_present = arg_present(hdr)

  if ~hdr_present then begin
    ;; A header was not provided, so read it from the file
    hdr = headfits(filename)
  endif else begin
    ;; A header was provided so we will not write it to the file, it
    ;; is up to the calling program to do it. Also, leave checksum
    ;; testing to the calling program because a new header will change
    ;; it anyway.
    checksum = 0
  endelse
  
  if keyword_set(test_checksum) then begin
    checksum_status = fits_test_checksum(hdr, /trust_datasum)
    if checksum_status eq -1 then begin
      red_message, ['CHECKSUM keyword does not have the' $
                    , 'correct value, indicating possible data' $
                    , 'corruption.']
      hgrep,hdr,'CHECK'
      red_message, 'Will not add keywords.'
      return
    endif
  endif
  
  instrument = strtrim(fxpar(hdr, 'INSTRUME'), 2)
  naxis = fxpar(hdr, 'NAXIS*')


  ;; * OBS_MODE -----------
  obs_mode = fxpar(hdr, 'OBS_MODE', count = cnt)
  if cnt eq 0 then begin
    ;; Temporal dimension, NAXIS[4], does not matter for this.
    case !true of

      ;; No Stokes, no tuning --> imaging : [x,y] or [x,y,t]
      naxis[4] eq 1 and naxis[3] eq 1 : obs_mode = "imaging"

      ;; Stokes but no tuning --> polarimetry : [x,y,s] or [x,y,s,t]
      naxis[4] gt 1 and naxis[3] eq 1 : obs_mode = "polarimetry"

      ;; Tuning but no Stokes -->, spectrometry : [x,y,lambda] or [x,y,lambda,t]
      naxis[4] eq 1 and naxis[3] gt 1 : obs_mode = "spectrometry"

      ;; Stokes and tuning --> spectropolarimetry : [x,y,s,lambda] or [x,y,lambda,s,t]
      naxis[4] gt 1 and naxis[3] gt 1 : obs_mode = "spectropolarimetry"

                                ; This should not happen!
      else : stop

    endcase
    red_fitsaddkeyword, hdr, 'OBS_MODE', obs_mode, anchor = anchor
  endif
  
  ;; * POINT_ID -----------

  point_it = fxpar(hdr, 'POINT_ID', count = cnt)
  if cnt eq 0 then begin
    ;; Add POINT_ID = DATE-OBS if not already in the header
    point_it = fxpar(hdr, 'DATE-OBS', count = cnt)
    red_fitsaddkeyword, hdr, 'POINT_ID', point_id, anchor = anchor
    print, 'POINT_ID', point_id
  endif
  
  
  ;; * CROTA ------------
  
  if n_elements(crota) eq 1 then begin
    ;; Add the CROTA keyword if the value was provided
    red_fitsaddkeyword, hdr, 'CROTA', crota, '[deg] Counterclockwise rotation of the FOV', anchor = anchor
    print, 'CROTA', crota
  endif

  ;;* CBBMINn/MAXn ------------------------- 
  ;;
  ;; BoundingBox keywords, extreme values of the WCS coordinates in
  ;; the file (handled by red_fitscube_addboundingbox).

  red_fitscube_addboundingbox, filename, hdr = hdr, anchor = anchor

  
  ;; The following keywords do not make sense for WB data

  if strmid(file_basename(filename), 0, 2) ne 'wb' then begin

    ;; Read WCS cordinates from the file
    red_fitscube_getwcs, filename, coordinates = coordinates
    tags = tag_names(coordinates)

    ;; * CSPMINn/MAXn ------------- 
    ;;
    ;; Coordinate SamPling/SPacing/SteP min and max, i.e., Min and max
    ;; tuning steps, to be read from the wavelength WCS coordinate.
    ;;
    ;; Will be used to set spectral_sampling_step_min/max.
  
    ;; Read WCS cordinates from the file
    red_fitscube_getwcs, filename, coordinates = coordinates
    tags = tag_names(coordinates)
    
    ;; Get the WCS coordinates in order
    lambda = coordinates.wave
    lambda = reform(lambda[0,0,*,0])
    dlambda = (red_differential(lambda))[1:*]

    ctype = strtrim(fxpar(hdr, 'CTYPE*'), 2)
    icoord = (where(ctype eq 'WAVE-TAB'))[0]

    red_fitsaddkeyword, hdr, 'CSPMIN'+strtrim(icoord+1, 2), min(dlambda) $
                        , '[nm] Min tuning step', anchor = anchor
    red_fitsaddkeyword, hdr, 'CSPMAX'+strtrim(icoord+1, 2), max(dlambda) $
                        , '[nm] Max tuning step', anchor = anchor

    print, 'CSPMIN'+strtrim(icoord+1, 2), min(dlambda)
    print, 'CSPMAX'+strtrim(icoord+1, 2), max(dlambda)


    ;; * RESOLVPW ------------- 
    ;; 
    ;; Resolving power of spectrograph, lambda/dlambda. Will be used
    ;; to set spectral_resolution_min/max
    ;;
    ;; dlambda is the FWHM of the etalon transmission at the current
    ;; wavelength. We get this by interpolation in data provided by
    ;; Jaime.
    
    case instrument of

      'CRISP' : begin
        ;; Profile data from Jaime in June 2026. fwhm_lambda is
        ;; wavelength, fwhm_crisp, fwhm_crisp2 are the profile
        ;; fwhms, all in Å.
        fwhmfile = cgFindPathTo('FWHM_CRISP-1_2.txt', SUCCESS=success)
        if ~success then stop
        READCOL, fwhmfile, fwhm_lambda, fwhm_crisp, fwhm_crisp2
        ;; We need nm, not Å, for both wavelengths and FWHMs.
        fwhm = interpol(fwhm_crisp/10., fwhm_lambda/10., mean(lambda)) ; [nm]
        if 0 then begin
          cgwindow
          ;; Plot mÅ vs nm as in the CRISPRED paper
          
          cgplot, /add,  fwhm_lambda/10.,fwhm_crisp2*1e3,color='blue' $
                  , xtitle = '$\lambda$ / 1 nm', ytitle = 'transmission FWHM / 1 m$\Angstrom$' $
                  , title = 'Based on FWHM_CRISP-1_2.txt'
          cgplot, /add, /over,fwhm_lambda/10.,fwhm_crisp*1e3,color='red'
          cglegend, /add, title = ['CRISP2','CRISP'], color = ['blue', 'red'] $
                    , location = [0.15, 0.85], align = 0 $
                    , vspace = 2
          cgcontrol, out = 'FWHM_CRISP-1_2.pdf'
          ;; Plot also RESOLVPW
          cgwindow
          cgplot, /add,  fwhm_lambda/10.,fwhm_lambda/fwhm_crisp,color='red' $
                  , xtitle = '$\lambda$ / 1 nm', ytitle = 'RESOLVPW' $
                  , title = 'Based on FWHM_CRISP-1_2.txt'
          cgplot, /add, /over,fwhm_lambda/10., fwhm_lambda/fwhm_crisp2,color='blue'
          cglegend, /add, title = ['CRISP2','CRISP'], color = ['blue', 'red'] $
                    , location = [0.17, 0.12], align = 3 $
                    , vspace = 2
          cgcontrol, out = 'RESOLVPW_CRISP-1_2.pdf'
        endif
      end

      'CRISP2' : begin
        ;; Profile data from Jaime in June 2026. fwhm_lambda is
        ;; wavelength, fwhm_crisp, fwhm_crisp2 are the profile
        ;; fwhms, all in Å.
        fwhmfile = cgFindPathTo('FWHM_CRISP-1_2.txt', SUCCESS=success)
        if ~success then stop
        READCOL, fwhmfile, fwhm_lambda, fwhm_crisp, fwhm_crisp2
        ;; We need nm, not Å, for both wavelengths and FWHMs.
        fwhm = interpol(fwhm_crisp2/10., fwhm_lambda/10., mean(lambda)) ; [nm]
      end

      'CHROMIS' : begin
        ;; We don't have the data for this yet.
      end

      else : stop

    endcase

    if n_elements(fwhm) gt 0 then begin
      ;; Give with 2 digits, FWHM has that at best
      RESOLVPW = red_round(mean(lambda)/fwhm, n_digits = 2) 
      red_fitsaddkeyword, hdr, 'RESOLVPW', RESOLVPW $
                          , 'Resolving power, lambda/FWHM(transmission)', anchor = anchor
      print, 'RESOLVPW', RESOLVPW
    endif
    
  endif

  if ~hdr_present then begin

    ;; Write the header to the file
    red_fitscube_newheader, filename, hdr

    ;; Test checkum
    if keyword_set(test_checksum) then begin
      hdr = headfits(destfile)
      fits_add_checksum, hdr    ; Not sure why this needs to be done again...
      red_fitscube_newheader, destfile, hdr
      checksum_status = fits_test_checksum(hdr, /trust_datasum)
      if checksum_status eq 1 then begin
        red_message, 'The new CHECKSUM is present and correct'
      endif else begin
        red_message, 'CHECKSUM incorrect or missing.'
        stop
      endelse
    endif
    
  endif
 
end
