; docformat = 'rst'

;+
; Get DATE-BEG, DATE-END, and DATE-AVG from FITS header, with default
; for DATE-AVG constructed from the other two.
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
; 
; 
; :Params:
; 
;    hdr : in, type=strarr
; 
;       A FITS header.
; 
; 
; :Keywords:
; 
;     date_beg : out, optional, type=string
;
;        DATE-BEG keyword from header.
;
;     date_end : out, optional, type=string
;
;        DATE-END keyword from header.
;
;     date_avg : out, optional, type=string
;
;        DATE-AVG keyword from header or reconstructed from the other
;        DATE-??? keywords.
;
;     count_beg : out, optional, type=integer
;
;        1 if DATE-BEG was found, 0 otherwise.
;
;     count_end : out, optional, type=integer
;
;        1 if DATE-END was found, 0 otherwise.
; 
;     count_avg : out, optional, type=integer
;
;        1 if DATE-AVG was found or reconstructed, 0 otherwise.
; 
;     comment_beg : out, optional, type=string
;
;        The DATE-BEG comment from the header.
;
;     comment_end : out, optional, type=string
;
;        The DATE-END comment from the header.
; 
;     comment_avg : out, optional, type=string
;
;        The DATE-AVG comment from the header, or a note that it was
;        constructed from the other DATE-??? keywords.
; 
;   
; :History:
; 
;    2017-06-13 : MGL. First version.
; 
; 
; 
; 
;-
pro red_fitspar_getdates, hdr $
                          , date_beg = date_beg $
                          , date_end = date_end $
                          , date_avg = date_avg $
                          , count_beg = count_beg $
                          , count_end = count_end $
                          , count_avg = count_avg $
                          , comment_beg = comment_beg $
                          , comment_end = comment_end $
                          , comment_avg = comment_avg

  ;; Read the existing data
  date_beg = fxpar(hdr, 'DATE-BEG', count = count_beg, comment = comment_beg)
  date_end = fxpar(hdr, 'DATE-END', count = count_end, comment = comment_end)
  date_avg = fxpar(hdr, 'DATE-AVG', count = count_avg, comment = comment_avg)  

  ;; Reconstruct missing avg if beg and end exist
  if count_avg eq 0 and count_beg eq 1 and count_end eq 1 then begin
    isodate = (strsplit(date_beg,'T',/extract))[0]
    time_beg = red_time2double((strsplit(date_beg,'T',/extract))[1])
    time_end = red_time2double((strsplit(date_end,'T',/extract))[1])
    time_avg = (time_beg+time_end)/2d
    date_avg = isodate + 'T' + red_timestring(time_avg)
    comment_avg = 'Average of DATE-BEG and DATE-END'
    count_avg = 1
  endif
  
end
