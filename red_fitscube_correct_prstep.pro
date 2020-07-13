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
  prmodes = fxpar(header, 'PRMODE*', count = Nsteps, comment = prmode_comments)
  
  ;; Loop through steps (starting from 1)
  for istep = 0, Nsteps-1 do begin

    oldprstep = prsteps[istep]

    ;; We can take PRPROC into account if we want to
    prproc = fxpar(header, 'PRPROC'+strtrim(istep+1, 2), count = Nproc)
    if Nproc eq 0 then prproc = ''
    
        
    case strtrim(prsteps[istep], 2)of

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
