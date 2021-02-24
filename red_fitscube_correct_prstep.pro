; docformat = 'rst'

;+
; Change old PRSTEP labels to SOLARNET-conforming labels and fix some
; other PRSTEP info.
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
;      The path to the file, the header of which we want to correct.
; 
; 
; :Keywords:
; 
;   nowrite : in, optional, type=boolean
;   
;      Don't write the header back to the file.
; 
; 
; :History:
; 
;    2020-05-12 : MGL. First version.
; 
;    2021-02-23 : MGL. Make sure there is a PRREF for make_nb_cube
;                 step with the WB cube file name.
; 
;-
pro red_fitscube_correct_prstep, filename $
                                 , header = header $
                                 , nowrite = nowrite $
                                 , verbose = verbose

  ;; Name of this subprogram
  inam = red_subprogram(/low, calling = inam1)
  
  ;; Read the header
  header = headfits(filename)

  prsteps = fxpar(header, 'PRSTEP*', count = Nsteps, comment = prstep_comments)
  prmodes = fxpar(header, 'PRMODE*', count = Nmodes, comment = prmode_comments)

  ;; This subroutine should not be run on exported files. Maybe check
  ;; for that here and return if it is?
  
  
  ;; Loop through steps (starting from 1)
  for istep = 0, Nsteps-1 do begin

    oldprstep = prsteps[istep]

    ;; We can take PRPROC into account if we want to
    prproc = fxpar(header, 'PRPROC'+strtrim(istep+1, 2), count = Nproc)
    if Nproc eq 0 then prproc = ''
    
        
    case strtrim(prsteps[istep], 2) of

      'Gain making' : begin 
        if prproc eq 'red::make_intdif_gains' then begin
          ;; red__make_intdif_gains.pro
          prsteps[istep] = 'RECIPROCAL' ; Plus other operations?
        endif else begin
          ;; chromis__makegains.pro, crisp__makegains.pro
          prsteps[istep] = 'RECIPROCAL' ; Plus other operations?
        endelse
      end 
      
      'Prepare WB science data cube' : $    ; red__make_wb_cube.pro
         prsteps[istep] = strjoin(['CONCATENATION' $ ; Combining multiple files 
                                   , 'SPATIAL-ALIGNMENT' $
                                   , 'DESTRETCHING'] $
                                  , ',')

      'Prepare NB science data cube' : begin
        if strmatch(prproc, 'make_nb_cube') then begin
          ;; chromis__make_nb_cube.pro, crisp__make_nb_cube.pro
          prsteps[istep] = strjoin(['CONCATENATION' $ ; Combining multiple files
                                    , 'SPATIAL-ALIGNMENT' $
                                    , 'DESTRETCHING' $
                                    , 'INTENSITY-CALIBRATION'] $
                                   , ',')
        endif else begin
          ;; chromis__make_scan_cube.pro, crisp__make_scan_cube.pro
          prsteps[istep] = strjoin(['CONCATENATION' $ ; Combining  multiple files
                                    , 'SPATIAL-ALIGNMENT' $
                                    , 'DESTRETCHING' $
                                    , 'CALIBRATION-INTENSITY-SPECTRAL'] $
                                   , ',')
        endelse
      end

      'Flat cubes' : $                    ; chromis__prepflatcubes.pro, crisp__prepflatcubes.pro
         prsteps[istep] = 'CONCATENATION' ; Combining multiple files 

      'Summing' : $             ; red__sumpolcal.pro, crisp__copy_oldsums.pro ;
         prsteps[istep] = 'SUMMING'
      
      'Dark summing' : $        ; red__sumdark.pro
         prsteps[istep] = 'SUMMING'
      
      'Flat summing' : $        ; red__sumflat.pro
         prsteps[istep] = 'SUMMING'
      
      'Pinhole summing' : $     ; red__sumpinh.pro      
         prsteps[istep] = 'SUMMING'

      'Demodulate' : $          ; crisp__demodulate.pro    
         prsteps[istep] = 'DEMODULATION'

      'Fourier-filtering' : $   ; crisp__demodulate.pro        
         prsteps[istep] = 'FIXED-PATTERN-REMOVAL' 

      'Make cavity free flats'    : $ ; red__fitgains.pro  
         prsteps[istep] = 'SPECTRAL-COORDINATE-DISTORTION-CORRECTION'
      
      'Interpolate' : $         ; red__fitscube_cmapcorr.pro
         prsteps[istep] = 'SPECTRAL-DISTORTION-CORRECTION'                                 
                                                    
      'Convert science data cube from LP format' : $ ; red__fitscube_convertlp.pro       
         prsteps[istep] = 'DATA-CURATION'
      
      'Cropping' : $            ; red__fitscube_crop.pro
         prsteps[istep] = 'CROPPING'
      
      'Crosstalk correction' : $ ; red__fitscube_crosstalk.pro       
         prsteps[istep] = 'STOKES-CROSSTALK-CORRECTION'
      
      'INTEGERIZATION'  : $     ; red__fitscube_integer.pro    
         prsteps[istep] = 'BZERO-BSCALE-TRUNCATION'                            

      'MOMFBD image restoration'  : begin $ ; red__prepmomfbd_fitsheaders.pro      
         prsteps[istep] = 'MOMFBD'              
        if strtrim(fxpar(header, 'INSTRUME'), 2) eq 'CRISP' then begin
          ;; Also remove "Phase Diversity" from prmodes[istep] if this
          ;; is a CRISP cube.
          pm = strsplit(prmodes[istep], ',', /extract)
          pos = where('Phase Diversity' eq pm, Nmatch)
          if Nmatch gt 0 then begin
            case 1 of
              n_elements(pm) eq 1     : prmodes[istep] = ''
              pos eq 0                : prmodes[istep] = strjoin( pm[1:*    ],               ',')
              pos eq n_elements(pm)-1 : prmodes[istep] = strjoin( pm[0:pos-1],               ',')
              else                    : prmodes[istep] = strjoin([pm[0:pos-1], pm[pos+1:*]], ',')
            endcase
            fxaddpar, header, 'PRMODE'+strtrim(istep+1, 2), prmodes[istep], strtrim(prmode_comments[istep], 2)
          endif
        endif

      end 

      else :

    endcase

    if keyword_set(verbose) then begin
      print
      if prsteps[istep] eq oldprstep then begin
        message, /info, 'No change.'
        message, /info, 'PRSTEP'+strtrim(istep+1, 2)+' : "'+prsteps[istep]+ '"'
        stop
      endif else begin
        message, /info, 'New PRSTEP'+strtrim(istep+1, 2)+' : "'+prsteps[istep]+ '"'
        message, /info, 'Old PRSTEP'+strtrim(istep+1, 2)+' : "'+oldprstep+'"'
      endelse
      message, /info, 'PRPROC'+strtrim(istep+1, 2)+' : ' + prproc
    endif

    if strmatch(prproc, '*::make_nb_cube') then begin
      ;; We want to make sure there is a PRREF keyword with the WB
      ;; cube name. We can get that from PRPARA.
      chars = ['', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j']
      for ichar = 0, n_elements(chars)-1 do begin
        tmp = fxpar(header, 'PRREF'+strtrim(istep+1, 2)+chars[ichar], count = Nref)
        if Nref eq 0 then break
        red_append, prrefs, tmp
      endfor                    ; ichar

      if n_elements(prrefs) eq 0 || max(strmatch(prrefs, 'Align reference:*')) eq 0 then begin
        ;; No match, we need to generate an appropriate keyword
        make_prref = 1
      endif else begin
        ;; Does the (last) align reference include a directory? It
        ;; should, unless this is an exported file (in which case this
        ;; subroutine should not have been run on it. Maybe check for
        ;; that at entry?).
        indx = where(strmatch(header, '*Align reference:*'), Nwhere)
        key = strtrim((strsplit(header[indx[-1]], '=', /extract))[0], 2)
        wcfile = red_strreplace(fxpar(header, key), 'Align reference: ', '')
        if file_dirname(wcfile) eq '.' then begin
          ;; Remove it and make a new one below
          red_fitsdelkeyword, header, key
          make_prref = 1
        endif else begin
          make_prref = 0
        endelse
      endelse

      if make_prref then begin
        ;; Make the PRREF if needed
        prpara = fxpar(header, 'PRPARA'+strtrim(istep+1, 2), count = Nparas, comment = prpara_comments)
        if Nparas eq 0 then stop
        prpara = json_parse(prpara)
        if ~prpara.haskey('WCFILE') then stop
        wcfile = prpara['WCFILE']
        if n_elements(prrefs) eq 0 then anchor = 'PRPARA'+strtrim(istep+1, 2) else anchor = prrefs[-1]
        stp = strtrim(istep+1, 2) + chars[ichar]
        red_fitsaddkeyword, header, 'PRREF'+stp $
                            , 'Align reference: '+wcfile $
                            , 'WB cube file name' $
                            , anchor = anchor
      endif
    endif
  
  ;; Rewrite the keyword in the header
  fxaddpar, header, 'PRSTEP'+strtrim(istep+1, 2), prsteps[istep], strtrim(prstep_comments[istep], 2)
  
  endfor                        ; istep

  if keyword_set(nowrite) then begin
    print
  endif else begin
    ;; Write the modified header to the file
    if ~array_equal(header, headfits(filename)) then red_fitscube_newheader, filename, header
  endelse

end

fname = '/scratch/mats/2016.09.19/CHROMIS-jan19/cubes_nb/nb_3950_2016-09-19T09:28:36_scans=0-3_corrected_im.fits'

red_fitscube_correct_prstep, fname, /verbose, /nowrite

end
