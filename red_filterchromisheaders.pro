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
;    2016-06-03 : MGL. Modify regex for framenumber so both dots and
;                 underscores work as delimiters. Only generate new
;                 framenumber if it is not already in the header.
;
;    2016-08-15 : MGL. Move metadata processing from
;                 red_filterchromisheaders to a new function,
;                 red_meta2head. Remove blank lines from new header
;                 before returning.
; 
; 
;-
function red_filterchromisheaders, head, silent=silent

  ;; This filtering should only be done for files without SOLARNET
  ;; header keyword.
  if fxpar(head, 'SOLARNET') ne 0 then return, head

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
                 ;; Just transfer all other keywords to the new
                 ;; header.
                 sxaddpar, newhead, name, value, comment, before = 'COMMENT'
              end
           endcase
        endif                   ; END

     endfor                     ; i

     
     ;; If there isn't already a DATE-END keyword, we'll have to
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
     
  endif else stop               ; non-SIMPLE

  newhead = newhead[where(newhead ne blanks(80))] ; Remove blank lines

  return, newhead

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
