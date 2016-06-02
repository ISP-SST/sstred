; docformat = 'rst'

;+
; Filter Guus' initial FITS headers to more or less SOLARNET standard
; keywords.
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, ISP
; 
; 
; :Returns:
; 
;    A corrected FITS header.
; 
; :Params:
; 
;    head : in, type=strarr
;
;       The header to correct
; 
; :Keywords:
;
;    metadata : in, type=struct
;
;	Additional metadata to include in the header
;	  Form: metadata.[keyword_name] = [keyword_data]
;
;    file : in, type=string
;
;	The file from which the header came.
;
; :History:
; 
;    2016-05-20 : MGL. First version.
;
;    2016-05-25 : JLF. Add keywords to the header.
;
;    2016-05-26 : MGL. Better interpretation/generation of DATE
;                 keywords. Get some info from the file name using
;                 regular expressions. Different filters in filter
;                 wheels for narrowband and wideband. Use tag_exist. 
;
;    2016-05-27 : MGL. Filter names in FILTERn keyword. Move the
;                 generation of DATE-END. Invented keywords for the
;                 scan number and for the number of the first frame,
;                 they should be replaced with proper keywords ASAP.
;                 Commented out some unfinished code.
;
;    2016-05-31 : MGL. Add filename to the header.
;
;    2016-05-31 : JLF. Added silent keyword (suppresses printing SOLARNET
;		  compliance message. 
;		  Started using red_keytab for keywords we aren't sure of.
;
;    2016-05-31 : MGL. Added filter wavelengths and FWHMs based on a
;                 table in an email from Aluxa. (If we knew the filter
;                 tilt angle, we could calculate the actual
;                 wavelength.)
;
;    2016-06-02 : MGL. New information can be added through the
;                 metadata also when the header is already SOLARNET
;                 compliant. 
; 
; 
;-
function red_filterchromisheaders, head, metadata=metaStruct, silent=silent

  if fxpar(head, 'SOLARNET') eq 0 then begin

     if ~keyword_set(silent) then begin
        print, 'red_filterchromisheaders : This header is not even partially SOLARNET compliant.'
        print, 'Make a new one.'
     endif

     Nlines = n_elements(head)

     ;; First make minimal header based on the old header

     simple = fxpar(head, 'SIMPLE', comment = simple_comment)
     if simple then begin

        bitpix = fxpar(head, 'BITPIX', comment = bitpix_comment)
        case bitpix of
           8: type = 1
           16: type = 2
           32: type = 3
           64: type = 14
           -32: type = 4
           -64: type = 5
           else: begin
              print, 'Invalid bitpix value: ', bitpix
              return, head
           end
        endcase

        naxisx = fxpar(head, 'NAXIS*')
        mkhdr, newhead, type, naxisx

        ;; Add SOLARNET keyword
        sxaddpar, newhead, 'SOLARNET', 0.5,  format = 'f3.1' $
                  , 'Fully SOLARNET-compliant=1.0, partially=0.5', before = 'COMMENT'
        
        ;; Add some more keywords:
        dummy = fxpar(newhead, 'TIMESYS',count=count)
        if count eq 0 then sxaddpar, newhead, 'TIMESYS', 'UTC', after = 'SOLARNET'
        dummy = fxpar(newhead, 'OBS_SHDU',count=count)
        if count eq 0 then $
           sxaddpar, newhead, 'OBS_SHDU', 1, 'This HDU contains observational data', after = 'SOLARNET'


        for i = 0, n_elements(head)-1 do begin
           
           name = strtrim((strsplit(head[i], '=', /extract))[0], 2)

           if name ne 'END' then begin

              value = fxpar(head, name, comment = comment)

              ;; Rewrite some keywords
              case name of
                 'DATE_OBS' : begin
                    ;; We'll use this as POINT-ID for now. It should
                    ;; really be something that we get from the turret or
                    ;; PIG systems so it is common to all instruments.
                    sxaddpar, newhead, 'POINT_ID', value, comment, before = 'COMMENT'
                 end
                 'EXPTIME' : begin
                    ;; Early versions of CHROMIS data used a non-SOLARNET
                    ;; keyword, rename it.
                    sxaddpar, newhead, 'XPOSURE',  value, comment, before = 'COMMENT' $
                              , format = 'f8.6'
                 end
                 'INTERVAL' : begin
                    ;; Early versions of CHROMIS data used a non-SOLARNET
                    ;; keyword, rename it.
                    sxaddpar, newhead, 'CADENCE',  value, comment, before = 'COMMENT' $
                              , format = 'f8.6'
                 end
                 else : begin
                    ;; Just transfer all other keywords to the new header.
                    sxaddpar, newhead, name,       value, comment, before = 'COMMENT'
                 end
              endcase
           endif                ; END

        endfor                  ; i

        ;; If there isn't already a DATE-END keyword, we'll have
        ;; construct one.
        dummy = fxpar(newhead, 'DATE-END', count=count)
        if count eq 0 then begin
           date_str = fxpar(newhead, 'DATE', count=count2)
           if count2 then sxaddpar, newhead, 'DATE-END', date_str, comment, after = 'DATE'
        endif
        dummy = fxpar(newhead, 'DATE-BEG', count=count)
        if count eq 0 then begin
           date_str = fxpar(newhead, 'DATE', count=count2)
           cadence = fxpar(newhead, 'CADENCE', count=cadcount)
           nframes = fxpar(newhead, 'NAXIS3', count=framecount)
           if count2 && cadcount && framecount then begin
              ;; Make DATE-BEG from DATE, counting backwards
              date = strsplit(date_str, 'T', /extract)
              time = red_time2double(date[1]) - cadence*nframes
              sxaddpar, newhead, 'DATE-BEG', date[0]+'T'+red_time2double(time, /dir), before = 'DATE-END'
           endif
        endif
        
     endif else stop            ; non-SIMPLE
  endif                         ; SOLARNET compliant

  ;; Add any additional metadata we were given.
  if n_elements(metaStruct) ne 0 then begin

     if n_elements(newhead) eq 0 then newhead = head

     if tag_exist(metaStruct, 'filename', /top_level) then begin


        ;; Extract metadata from the filename (for now it's all we
        ;; have)
        
        ;; Strip directory and extension if needed
        filename = file_basename(metaStruct.filename)
        barefile = file_basename(metaStruct.filename, '.fits')

        dummy = fxpar(newhead, 'FILENAME', count=count)
        if count eq 0 then begin
           ;; Preserve existing file name, this keyword should be the
           ;; name of the original file.
           metaStruct.filename = filename
           sxaddpar, newhead, 'FILENAME', filename $
                     , 'Name of original file.', before = 'COMMENT'
        endif

        ;; Camera tag (check that this is the correct keyword...)
        dummy = fxpar(newhead, red_keytab('camtag'), count=count)
        if count eq 0 then begin

           ;; The camera tag consists of the string 'cam' followed by a
           ;; roman number.
           camtag = ((stregex(barefile, '(\.|^)(cam[IVXL]+)(\.|$)', /extr, /subexp))[2,*])[0]

           if camtag ne '' then sxaddpar, newhead, red_keytab('camtag'), camtag $
                                          , 'Inferred from filename.' $
                                          , before = 'COMMENT'
        endif                   ; CAMERA
        
        
        ;; FILTERn, WAVEBAND and WAVELNTH
        dummy = fxpar(newhead, red_keytab('pref'), count=count)
        if count eq 0 then begin

           ;; Get prefilter names from file name. For now, we have
           ;; the filter wheel position as a tag consisting of the
           ;; letter w followed by a single digit.
           wheelpos = ((stregex(barefile, '(\.|^)(w[0-9]{1})(\.|$)', /extr, /subexp))[2,*])[0]

           if wheelpos ne '' then begin
              
              channel = fxpar(head, red_keytab('cam_channel'), count=count)
              if count eq 0 then begin
                                ;  TDB: can we extract filter-info without the channel-tag?
                 sxaddpar, newhead, red_keytab('pref'), wheelpos, before = 'COMMENT'
              endif else begin

                 if channel eq 'Chromis-N' then begin

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
                          wavelnth = 397.7e-9
                          fwhm = 0.37e-9
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

           endif                ; Found a wheel position


        endif                   ; FILTER1, WAVEBAND

        ;; WAVELNTH
        prefilter = fxpar(newhead, red_keytab('pref'), count=count)
        if count gt 0 then begin
           if fxpar(newhead, 'WAVELNTH') eq 0.0 then begin

              case fxpar(newhead, red_keytab('pref')) of
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
                    wavelnth = 484.7-e-9
                    fwhm = 0.6e-9
                 end
                 'CaHK-cont' : begin
                    wavelnth = 395.1e-9
                    fwhm = 1.3e-9
                 end
              endcase

              sxaddpar, newhead, 'WAVELNTH', strtrim(wavelnth, 2) $
                        , '[m] Prefilter peak wavelength', after = red_keytab('pref')
              ;; Add also the FWHM in a keyword?

           endif 
        endif                   ; WAVELNTH
        
        ;; The frame number is the last field iff it consists
        ;; entirely of digits. The third subexpression of the regular
        ;; expression matches only the end of the string because
        ;; that's where it is if it is present. We do not know the
        ;; length of the frame number field so if the third
        ;; subexpression were allowed to match a dot we would get
        ;; false matches with the scan and prefilter fields.
        ;;
        ;; Actually, this is the file number. Let's make the frame
        ;; number (of the first frame) 10000 times this.
        ;;
        filenumber = ((stregex(barefile, '(\.)([0-9]+)($)', /extr, /subexp))[2,*])[0]
        if filenumber ne '' then begin
           framenumber1 = 10000*long(filenumber)
           ;;
           ;; But what keyword to use for this? Invent something
           ;; for now!
           sxaddpar, newhead, red_keytab('frame'), framenumber1 $
                     , '(number of first frame) Inferred from filename.' $
                     , before = 'COMMENT'
        endif

        ;; The scan number is the only field that is exactly five
        ;; digits long:
        scannumber = ((stregex(barefile, '(\.|^)([0-9]{5})(\.|$)', /extr, /subexp))[2,*])[0]
        if scannumber ne '' then sxaddpar, newhead, red_keytab('scannumber'), long(scannumber) $
                                           , '(scan number) Inferred from filename.' $
                                           , before = 'COMMENT'
        ;; What keyword to use for this? Invented SCANNUM for now.


        ;; Don't trust positions if we don't have to.
        ;; states = strsplit(barefile,'.',/extract)
        ;; camtag = states[0]
        ;; fscan = states[1]
        ;; pref = states[2]
        ;; frame = states[3]

     endif                      ; filename

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

  endif else begin              ; metaStruct

     if ~keyword_set(silent) then begin
        print, 'red_filterchromisheaders : No changes, return old header.'
     endif

     return, head

  endelse


end

; /mnt => /storage (in case you aren't on polar)
fname = '/storage/sand15n/Incoming/2016.05.13/AR2542/15:57:10/Chromis-W/camXXVIII.00019.w2.0000076.fits'
head = red_readhead(fname)

print, head
print

head2 = red_filterchromisheaders(head)

print, head2

;head3 = red_filterchromisheaders(head2)
;
;print, head3
;
end
