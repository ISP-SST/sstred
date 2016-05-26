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
;                 regular expressions.
; 
; 
;-
function red_filterchromisheaders, head, metadata=metaStruct
  
  if fxpar(head, 'SOLARNET') gt 0 then begin
     print, 'red_filterchromisheaders : This header should already be (at least partially) SOLARNET compliant.'
     print, 'Return without changes.'
     return, head
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

     naxis  = fxpar(head, 'NAXIS', comment = naxis_comment)
     naxisx = lonarr(naxis)

     for i = 0, naxis-1 do naxisx[i] = fxpar(head, 'NAXIS'+strtrim(i+1, 2))
     
     mkhdr, newhead, type, naxisx

     for i = 0, n_elements(head)-1 do begin
        
        name = strtrim((strsplit(head[i], '=', /extract))[0], 2)

        if name ne 'END' then begin

           value = fxpar(head, name, comment = comment)

           case name of
              'DATE' : begin
                 ;; Rewrite the automatically generated DATE keyword
                 sxaddpar, newhead, 'DATE', value, comment, before = 'COMMENT'
                 ;; If there isn't already a DATE-END keyword,
                 ;; this is the best we have for it.
                 if fxpar(head, 'DATE-END') eq '' then $
                    sxaddpar, newhead, 'DATE-END', value, comment, before = 'COMMENT'
              end
              'DATE_OBS' : begin
                 ;; We'll use this as POINT-ID for now. It
                 ;; should really be something that we get from the
                 ;; turret or PIG systems so it is common to all
                 ;; instruments. 
                 sxaddpar, newhead, 'POINT_ID', value, comment, before = 'COMMENT'
              end
              'EXPTIME' : begin
                 ;; Early versions of CHROMIS data used a non-SOLARNET
                 ;; keyword, rename it.
                 sxaddpar, newhead, 'XPOSURE',  value, comment, before = 'COMMENT', format = 'f8.6'
              end
              'INTERVAL' : begin
                 ;; Early versions of CHROMIS data used a non-SOLARNET
                 ;; keyword, rename it.
                 sxaddpar, newhead, 'CADENCE',  value, comment, before = 'COMMENT', format = 'f8.6'
              end
              else : begin
                 ;; Just transfer all other keywords to the new header.
                 sxaddpar, newhead, name,       value, comment, before = 'COMMENT'
              end
           endcase
        endif                   ; END

     endfor                     ; i

     if fxpar(head, 'DATE-BEG') eq '' then begin
        ;; Make DATE-BEG from DATE, counting backwards
        date = strsplit(fxpar(head, 'DATE'), 'T', /extract)
        time = red_time2double(date[1]) - fxpar(head, 'INTERVAL')*fxpar(head, 'NAXIS3')
        sxaddpar, newhead, 'DATE-BEG', date[0]+'T'+red_time2double(time, /dir), after = 'DATE'
     endif
     
     ;; Add any additional metadata we were given.
     if n_elements(metaStruct) ne 0 then begin
      keys = tag_names(metaStruct)
     
      ;; check if the filename is present
      ifile = where(keys eq strupcase('filename'),cnt)
      ;; extract metadata from the filename (for now it's all we have)
      if cnt ne 0 then begin

         ;; Strip directory and extension if needed
         filename = file_basename(metaStruct.(ifile),)
         barefile = file_basename(metaStruct.(ifile), '.fits')
         metaStruct.(ifile) = filename

         ;; Camera tag (check that this is the correct keyword...)
         if fxpar(newhead, 'CAMERA') eq '' then begin

            ;; The camera tag consists of the string 'cam' followed by a
            ;; roman number.
            camtag = (stregex(barefile, '(\.|^)(cam[IVXL]+)(\.|$)', /extr, /subexp))[2,*]

            if camtag ne '' then sxaddpar, newhead, 'CAMERA', camtag $
                                           , 'Inferred from file name.' $
                                           , before = 'COMMENT'

         endif                  ; CAMERA
         
         
         ;; WAVEBAND and WAVELNTH
         if fxpar(newhead, 'WAVEBAND') eq '' then begin

            ;; Get prefilter names from file name. For now, we have
            ;; the filter wheel position as a tag consisting of the
            ;; letter w followed by a single digit.
            wheelpos = (stregex(barefile, '(\.|^)(w[0-9]{1})(\.|$)', /extr, /subexp))[2,*]
            if wheelpos ne '' then begin

               ;; The filter names here are for Chromis-N. For
               ;; Chromis-W and Chromis-D, w1-w5 means 'CaHK-cont' and
               ;; w6 means 'Hb-cont'.
               case wheelpos of
                  'w1' : begin
                     waveband = 'CaK-blue'
                     wavelnth = 0.
                  end
                  'w2' : begin
                     waveband = 'CaK-core'
                     wavelnth = 0.
                  end
                  'w3' : begin
                     waveband = 'CaH-core'
                     wavelnth = 0.
                  end
                  'w4' : begin
                     waveband = 'CaH-red'
                     wavelnth = 0.
                  end
                  'w5' : begin
                     waveband = 'CaH-cont'
                     wavelnth = 0.
                  end
                  'w6' : begin
                     waveband = 'Hb-core'
                     wavelnth = 0.
                  end
               endcase

               comment = 'Inferred from filter wheel position in file name.'
               sxaddpar, newhead, 'WAVEBAND', waveband, comment, before = 'COMMENT'
               sxaddpar, newhead, 'WAVELNTH', wavelnth, '[nm]', before = 'COMMENT'

            endif               ; Found a wheel position

         endif                  ; WAVEBAND & WAVELNTH
         
         ;; The frame number is the last field iff it consists
         ;; entirely of digits. The third subexpression of the regular
         ;; expression matches only the end of the string because
         ;; that's where it is if it is present. We do not know the
         ;; length of the frame number field so if the third
         ;; subexpression were allowed to match a dot we would get
         ;; false matches with the scan and prefilter fields.
         ;;
         ;; Frame number = (stregex(barefile, '(\.)([0-9]+)($)', /extr, /subexp))[2,*]
         ;;
         ;; Actually, this is the file number. Let's make the frame
         ;; number (of the first frame) 10000 times this.
         ;;
         ;; But what keyword to use for this?


         ;; The scan number is the only field that is exactly five
         ;; digits long:
         scannumber = (stregex(fname, '(\.|^)([0-9]{5})(\.|$)', /extr, /subexp))[2,*]
         if scannumber ne '' then sxaddpar, newhead, 'XXX', scannumber $
                                            , 'Inferred from file name.' $
                                            , before = 'COMMENT'
         ;; What keyword to use for this?


         ;; Don't trust positions if we don't have to.
         ;; states = strsplit(barefile,'.',/extract)
         ;; camtag = states[0]
         ;; fscan = states[1]
         ;; pref = states[2]
         ;; frame = states[3]


         ;; add to the metadata if the appropriate keywords aren't there,
         ;; don't overwrite metadata we received
	
         ;; check these against solarnet standard
         ;; waveband seems to be the right one for pref-filter 
         ;; not sure about the others.
         nametags = ['camtag','scannum','waveband','framenum']
         for i = 0, 3 do begin
	  junk = where(keys eq strupcase(nametags[i]),cnt)
	  
	  if cnt eq 0 then begin
	    keys = [keys,nametags[i]]
	    metaStruct = create_struct(metaStruct,nametags[i],states[i])
	  endif
	
	endfor
      endif 		; filename
     
      for i = 0, n_elements(keys)-1 do $
	sxaddpar,newhead,keys[i],metaStruct.(i),before='COMMENT'
     endif		; metaStruct
     
     ;; Add SOLARNET keyword
     sxaddpar, newhead, 'SOLARNET', 0.5,  format = 'f3.1' $
               , 'Fully SOLARNET-compliant=1.0, partially=0.5', before = 'COMMENT'

  endif else stop               ; Non-simple

  return, newhead

end

; /mnt => /storage (in case you aren't on polar)
fname = '/storage/sand15n/Incoming/2016.05.13/AR2542/15:57:10/Chromis-W/camXXVIII.00019.w2.0000076.fits'
head = red_readhead(fname)

print, head
print

head2 = red_filterchromisheaders(head)

print, head2

head3 = red_filterchromisheaders(head2)

print, head3

end
