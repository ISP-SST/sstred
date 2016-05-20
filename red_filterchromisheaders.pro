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
; :History:
; 
;    2016-05-20 : MGL. First version.
; 
; 
;-
function red_filterchromisheaders, head
  
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
     if bitpix eq 16 then type = 2 else stop ; More cases here!

     naxis  = fxpar(head, 'NAXIS', comment = naxis_comment)
     naxisx = lonarr(naxis)

     for i = 0, naxis-1 do naxisx[i] = fxpar(head, 'NAXIS'+strtrim(i+1, 2))
     
     mkhdr, newhead, type, naxisx

     for i = 0, n_elements(head)-1 do begin
        
        name = strtrim((strsplit(head[i], '=', /extract))[0], 2)

        if name ne 'END' then begin

           value = fxpar(head, name, comment = comment)

           case name of
              'DATE'     : sxaddpar, newhead, 'DATE-BEG', value, comment, before = 'COMMENT'
              'DATE_OBS' : sxaddpar, newhead, 'POINT_ID', value, comment, before = 'COMMENT'
              'EXPTIME'  : sxaddpar, newhead, 'XPOSURE',  value, comment, before = 'COMMENT', format = 'f8.6'
              'INTERVAL' : sxaddpar, newhead, 'CADENCE',  value, comment, before = 'COMMENT', format = 'f8.6'
              else       : sxaddpar, newhead, name,       value, comment, before = 'COMMENT'
           endcase
        endif                   ; END

     endfor                     ; i

     ;; Make DATE-END
     dbeg = strsplit(fxpar(head, 'DATE'), 'T', /extract)
     date = dbeg[0]
     time = red_time2double(dbeg[1]) + fxpar(head, 'INTERVAL')*fxpar(head, 'NAXIS3')
     dend = date+'T'+red_time2double(time, /dir)
     sxaddpar, newhead, 'DATE-END', dend, after = 'DATE-BEG'

     ;; Add SOLARNET keyword
     sxaddpar, newhead, 'SOLARNET', 0.5,  format = 'f3.1' $
               , 'Fully SOLARNET-compliant=1.0, partially=0.5', before = 'COMMENT'

  endif else stop               ; Non-simple

  return, newhead

end


fname = '/mnt/sand15n/Incoming/2016.05.13/AR2542/15:57:10/Chromis-W/camXXVIII.00019.w2.0000076.fits'
head = red_readhead(fname)

print, head
print

head2 = red_filterchromisheaders(head)

print, head2

head3 = red_filterchromisheaders(head2)

print, head3

end
