; docformat = 'rst'

;+
; Calculates the combined FITS header for files summed with
; red_sumfiles or rdx_sumfiles.
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl
; 
; 
; :Returns:
; 
;    The combined header.
; 
; :Params:
; 
;    files : in, type=strarr
;
;      The names of the summed files.
; 
;    sum : in, type=array
;
;      The array representing the summed files.
;
; :Keywords:
; 
;    Nsum : in, optional, type=integer
;      
;      The number of summed frames. Specify this if some were excluded
;      from the sum.
; 
; 
; :History:
; 
;    2017-03-13 : MGL. First version.
; 
;-
function red_sumheaders, files, sum $
                         , Nsum = Nsum
    
  Nfiles = n_elements(files)
  
  for ifile = 0, Nfiles-1 do begin

      head = red_readhead(files[ifile])

      tab_hdus = fxpar(head, 'TAB_HDUS')
      if tab_hdus ne '' then begin

        tab = readfits(files[ifile], theader, /exten, /silent)
        red_append, date_beg_array, ftget(theader,tab, 'DATE-BEG') 

      endif
      
  endfor                        ; ifile

  if n_elements(Nsum) eq 0 then begin
    Nsum = n_elements(date_beg_array)
  endif else begin
    if n_elements(date_beg_array) ne Nsum then begin
      ;; If we could get info about the rejected frame numbers from
      ;; both red_sumfiles and rdx_sumfiles we could use that here!
      print, 'red_sumheaders : n_elements(date_beg_array) ne Nsum'
      print, '                 The DATE_??? keywords may be slightly off.' 
    endif
  endelse

  ;; Assume all file headers are the same, except for frame numbers,
  ;; timestamps, etc. We use info from the latest read header and
  ;; modify it where necessary.
  
  ;; Update header with respect to summed array
  check_fits, sum, head, /UPDATE, /SILENT        

  ;; New date, at better position
  sxdelpar, head, 'DATE'
  sxaddpar, head, 'DATE', red_timestamp(/iso) $
            , ' Creation UTC date of FITS header ', after = 'SOLARNET'

  ;; Remove some irrelevant keywords
  sxdelpar, head, 'TAB_HDUS'    ; No tabulated headers in summed file
  sxdelpar, head, 'FRAMENUM'    ; No particular frame number
  sxdelpar, head, 'SCANNUM'     ; No particular scan number
  sxdelpar, head, 'FILENAME'    ; Add this later!
  sxdelpar, head, 'CADENCE'     ; Makes no sense to keep?

  ;; Exposure times etc.
  exptime = sxpar(head, 'XPOSURE', count=count, comment=exptime_comment)
  sxaddpar, head, 'XPOSURE', nsum*exptime, ' [s] Total exposure time'
  sxaddpar, head, 'TEXPOSUR', exptime $
            , ' [s] Single-exposure time', before = 'XPOSURE'
  if nsum gt 1 then $
     sxaddpar, head, 'NSUMEXP', nsum $
               , ' Number of summed exposures', before = 'XPOSURE'
  
  ;; DATE-??? keywords, base them on the tabulated timestamps
  date_obs = (strsplit(fxpar(head, 'DATE-OBS'), 'T',/extract))[0]
  times = red_time2double(strmid(DATE_BEG_ARRAY,11))
  time_beg = red_timestring(min(times))
  time_end = red_timestring(max(times)+exptime)
  time_ave = red_timestring(mean(times)+exptime/2)
  if n_elements(time_end) ne 0 then $
     sxaddpar, head, 'DATE-END', date_obs+'T'+time_end $
               , ' Date of end of observation', after = 'DATE'
  if n_elements(time_ave) ne 0 then $
     sxaddpar, head, 'DATE-AVE', date_obs+'T'+time_ave $
               , ' Average date of observation', after = 'DATE'
  if n_elements(time_beg) ne 0 then $
     sxaddpar, head, 'DATE-BEG', date_obs+'T'+time_beg $
               , ' Date of start of observation', after = 'DATE'

  ;; Add "global" metadata
  red_metadata_restore, head

  return, head

end
