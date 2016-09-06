; docformat = 'rst'

;+
; Add header info from meta data.
; 
; :Categories:
;
;    SST pipelines
; 
; 
; :Author:
; 
;    Mats Löfdahl, ISP
; 
; 
; :Returns:
; 
; 
; :Params:
; 
;    head : in, type=strarr
;
;       The header to which we want to add information. 
; 
; 
; :Keywords:
; 
;    metadata : in, type=struct
;
;	Additional metadata to include in the header
;	Form: metadata.[keyword_name] = [keyword_data]
;   
; 
; 
; :History:
;
;    2016-08-15 : MGL. Take relevant code from
;                 red_filterchromisheaders. Remove .momfbd extension
;                 from variable barefile. Recognize underscore as a
;                 delimiter when parsing filenames.
;
;    2016-08-23 : THI. Rename camtag to detector and channel to camera,
;                 so the names match those of the corresponding SolarNet
;                 keywords.
; 
;-
function red_meta2head, head, metadata=metaStruct

  ;; Was any meta data provided?
  if n_elements(metaStruct) eq 0 then return, head

  newhead = head

  ;; Extract metadata from the filename (for now it's all we have)
  if tag_exist(metaStruct, 'filename', /top_level) then begin

     ;; Strip directory...
     filename = file_basename(metaStruct.filename)
     ;; ...and extension if needed
     barefile = filename
     barefile = file_basename(barefile, '.fits')
     barefile = file_basename(barefile, '.momfbd')

     dummy = fxpar(newhead, 'FILENAME', count=count)
     if count eq 0 then begin
        ;; Preserve existing file name, this keyword should be the
        ;; name of the original file.
        metaStruct.filename = filename
        sxaddpar, newhead, 'FILENAME', filename $
                  , 'Name of original file.', before = 'COMMENT'
     endif

     ;; Detector ID (check that this is the correct keyword...)
     dummy = fxpar(newhead, red_keytab('detector'), count=count)
     if count eq 0 then begin

        ;; The camera tag consists of the string 'cam' followed by a
        ;; roman number.
        detector = ((stregex(barefile, '(_|\.|^)(cam[IVXL]+)(_|\.|$)', /extr, /subexp))[2,*])[0]

        if detector ne '' then sxaddpar, newhead, red_keytab('detector'), detector $
                                       , 'Inferred from filename.' $
                                       , before = 'COMMENT'
     endif                      ; CAMERA
     
     
     ;; FILTERn, WAVEBAND and WAVELNTH
     dummy = fxpar(newhead, red_keytab('pref'), count=count)
     if count eq 0 then begin

        ;; Get prefilter names from file name. For now, we have
        ;; the filter wheel position as a tag consisting of the
        ;; letter w followed by a single digit.
        wheelpos = ((stregex(barefile, '(_|\.|^)(w[0-9]{1})(_|\.|$)', /extr, /subexp))[2,*])[0]

        if wheelpos ne '' then begin
           
           camera = fxpar(head, red_keytab('camera'), count=count)
           if count eq 0 then begin
                                ;  TDB: can we extract filter-info without the camera-tag?
              sxaddpar, newhead, red_keytab('pref'), wheelpos, before = 'COMMENT'
           endif else begin

              if camera eq 'Chromis-N' then begin

                 ;; Chromis-N
                 case wheelpos of
                    'w1' : begin
                       filter1 = 'CaK-blue'
                       waveband = 'Ca II H & K'
                    end
                    'w2' : begin
                       filter1 = 'CaK-core'
                       waveband = 'Ca II H & K'
                    end
                    'w3' : begin
                       filter1 = 'CaH-core'
                       waveband = 'Ca II H & K'
                    end
                    'w4' : begin
                       filter1 = 'CaH-red'
                       waveband = 'Ca II H & K'
                    end
                    'w5' : begin
                       filter1 = 'CaH-cont'
                       waveband = 'Ca II H & K'
                    end
                    'w6' : begin
                       filter1 = 'Hb-core'
                       waveband = 'H-beta'
                    end
                 endcase
              endif else begin

                 ;; Chromis-W and Chromis-D
                 case wheelpos of
                    'w6' : begin
                       filter1 = 'Hb-cont'
                       waveband = 'H-beta'
                    end
                    else: begin
                       filter1 = 'CaHK-cont'
                       waveband = 'Ca II H & K'
                    end
                 endcase
              endelse
              
              comment = 'Inferred from filter wheel position in filename.'
              sxaddpar, newhead, red_keytab('pref'), filter1, comment, before = 'COMMENT'
              
           endelse

        endif                   ; Found a wheel position


     endif                      ; FILTER1, WAVEBAND

     ;; WAVELNTH
     prefilter = fxpar(newhead, red_keytab('pref'), count=count)
     if count gt 0 then begin
        if fxpar(newhead, 'WAVELNTH') eq 0.0 then begin

           case strtrim(prefilter,2) of
              'CaK-blue' : begin
                 wavelnth = 392.55e-9
                 fwhm = 0.37e-9
              end
              'CaK-core' : begin
                 wavelnth = 393.44e-9
                 fwhm = 0.37e-9
              end
              'CaH-core' : begin
                 wavelnth = 396.92e-9
                 fwhm = 0.37e-9
              end
              'CaH-red' : begin
                 wavelnth = 397.7e-9
                 fwhm = 0.37e-9
              end
              'CaH-cont' : begin
                 wavelnth = 400.01e-9
                 fwhm = 0.37e-9
              end
              'Hb-core' : begin
                 wavelnth = 486.20e-9
                 fwhm = 0.44e-9
              end
              'Hb-cont' : begin
                 wavelnth = 484.7e-9
                 fwhm = 0.6e-9
              end
              'CaHK-cont' : begin
                 wavelnth = 395.1e-9
                 fwhm = 1.3e-9
              end
              else: begin
		wavelnth = 'unknown'
		fwhm = 'unknown'
	      end
           endcase

           sxaddpar, newhead, 'WAVELNTH', strtrim(wavelnth, 2) $
                     , '[m] Prefilter peak wavelength', after = red_keytab('pref')
           ;; Add also the FWHM in a keyword?

        endif 
     endif                      ; WAVELNTH
        
     framenumber = fxpar(head, red_keytab('frame'), count = count)
     if count eq 0 then begin
        ;; The frame number is the last field iff it consists
        ;; entirely of digits. The third subexpression of the
        ;; regular expression matches only the end of the string
        ;; because that's where it is if it is present. We do not
        ;; know the length of the frame number field so if the
        ;; third subexpression were allowed to match a dot we would
        ;; get false matches with the scan and prefilter fields.
        ;;
        ;; Actually, this is the file number. Let's make the frame
        ;; number (of the first frame) 10000 times this.
        ;;
        filenumber = ((stregex(barefile, '(\.|_)([0-9]+)($)', /extr, /subexp))[2,*])[0]
        if filenumber ne '' then begin
           framenumber1 = 10000L*long(filenumber)
           ;;
           ;; But what keyword to use for this? Invent something
           ;; for now!
           sxaddpar, newhead, red_keytab('frame'), framenumber1 $
                     , '(number of first frame) Inferred from filename.' $
                     , before = 'COMMENT'
        endif
     endif 

     ;; The scan number is the only field that is exactly five
     ;; digits long:
     scannumber = ((stregex(barefile, '(_|\.|^)([0-9]{5})(_|\.|$)', /extr, /subexp))[2,*])[0]
     if scannumber ne '' then sxaddpar, newhead, red_keytab('scannumber'), long(scannumber) $
                                        , '(scan number) Inferred from filename.' $
                                        , before = 'COMMENT'

  endif                         ; filename

  ;; add to the metadata if the appropriate keywords aren't there,
  ;; don't overwrite metadata we received
  
;        ;; check these against solarnet standard
;        ;; waveband seems to be the right one for pref-filter 
;        ;; not sure about the others.
;        keys = tag_names(metaStruct)
;        nametags = ['camtag','scannum','waveband','framenum']
;        for i = 0, 3 do begin
;           junk = where(keys eq strupcase(nametags[i]),cnt)
;           
;           if cnt eq 0 then begin
;              keys = [keys,nametags[i]]
;              metaStruct = create_struct(metaStruct,nametags[i],states[i])
;           endif
;           
;        endfor

     
;        for ikey = 0, n_elements(keys)-1 do begin
;           if somecheck(keys[i]) then begin
;              sxaddpar,newhead,keys[i],metaStruct.(i),before='COMMENT'
;           endif
;        endfor                  ; ikey

  return, newhead

end
