; docformat = 'rst'

;+
; Return the header from a FITS format file, taking some
; pipeline-specific issues into account. 
; 
; :Categories:
;
;    SST pipeline
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
;    fname : in, type=string
;
;       The name of the data file.
; 
; :History:
; 
;    2017-03-10 : MGL. Moved reading of headers from ANA fz format
;                 files from red_readhead.pro.
; 
;    2017-03-13 : MGL. Deal with camera software OBS_PHDU bug.
; 
; 
; 
;-
function red_readhead_fits, fname, $
                            framenumber = framenumber, $
                            silent = silent, $
                            extension = extension

  compile_opt idl2

  ;; Are we here for an extension or the primary HDU?
  if n_elements(extension) eq 0 then begin
    ;; primary
    red_rdfits, fname, header = header

    if fxpar(header, 'SOLARNET') eq 0 then begin
      caminfo = red_camerainfo( red_detectorname(fname,head=header) )
      if strmatch(caminfo.model,'PointGrey*') then begin 
        ;; We could add a date check here as well.

        ;; Old PointGrey data header filtering to bring it
        ;; to solarnet compliance. 
        header = red_filterchromisheaders(header, silent=silent)
      endif
    endif

    ;; Quick and dirty table reading to get date-beg and
    ;; date-end. Should rewrite this to propagate the whole
    ;; table through the pipeline! /MGL
    date_beg = sxpar(header, 'DATE-BEG', count = Nbeg)
    if Nbeg eq 0 then begin
      tab_hdus = fxpar(header, 'TAB_HDUS')
      if tab_hdus ne '' then begin
        
        ;; At some point, implement removing bad frames from
        ;; the tabulated list. /MGL

        tab = readfits(fname, theader, /exten, /silent)
        date_beg_array = ftget(theader,tab, 'DATE-BEG') 
        
        fxaddpar, header, 'DATE-BEG', date_beg_array[0] $
                  , 'First in DATE-BEG table.', after = 'DATE'
        
        isodate = (strsplit(date_beg_array[0], 'T', /extract))[0]
        time_beg_array = red_time2double(red_strreplace(date_beg_array,isodate+'T',''))

        date_end = sxpar(header, 'DATE-END', count = Nend)
        if Nend eq 0 then begin
          ;; time_end = time_beg_array[-1] + sxpar(header, 'XPOSURE')
          time_end = time_beg_array[n_elements(time_beg_array)-1] + sxpar(header, 'XPOSURE')
          date_end = isodate + 'T' + red_time2double(time_end, /inv)
          fxaddpar, header, 'DATE-END', date_end $
                    , 'Last in DATE-BEG table + XPOSURE.' $
                    , after = 'DATE-BEG'
        endif

        date_avg = sxpar(header, 'DATE-AVG', count = Navg)
        if Navg eq 0 then begin
          time_avg = mean(time_beg_array) + sxpar(header, 'XPOSURE')/2.
          date_avg = isodate + 'T' + red_time2double(time_avg, /inv)
          sxaddpar, header, 'DATE-AVG', date_avg $
                    , 'Average of DATE-BEG table + XPOSURE/2.' $
                    , after = 'DATE-BEG'
        endif

      endif
    endif

    ;; Hack to get the prefilter from the file name in data
    ;; from 2016.08.30.
    pref = fxpar( header, red_keytab('pref'), count=count )
    if count eq 0 then begin
      state = fxpar(header, 'STATE', count=count )
      if count eq 0 then begin
        ;; Try to read from file name
        fname_split = strsplit(file_basename(fname,'.fits'),'_',/extr)
        if n_elements(fname_split) gt 4 then begin
          ;; Shorter and it might be a dark
          prefilter = (fname_split)[-1]
          ;; Translate to previously used filter names
          case prefilter of
            'hbeta-core' : prefilter = '4862'
            'hbeta-cont' : prefilter = '4846'
            'cah-core'   : prefilter = '3969'
            else:
          endcase
          fxaddpar, header, red_keytab('pref'), prefilter, 'Extracted from file name'
        endif
      endif else begin          ; STATE keyword exists but not FILTER1 (e.g. 2016.08.31 data)
        state_split = strsplit( state, '_',  /extr )
        if n_elements(state_split) gt 1 then begin
          state1 = state_split[0] ;  for 2016.08.31:  state = 'wheel00002_hrz32600'
          camera = fxpar(header, red_keytab('camera'), count=count)
          if count gt 0 then begin
            if camera eq 'Chromis-N' then begin

              ;; Chromis-N
              case state1 of
                'wheel00001' : prefilter = '3925' ; Ca II K blue wing
                'wheel00002' : prefilter = '3934' ; Ca II K core
                'wheel00003' : prefilter = '3969' ; Ca II H core
                'wheel00004' : prefilter = '3978' ; Ca II H red wing
                'wheel00005' : prefilter = '3999' ; Ca II H continuum
                'wheel00006' : prefilter = '4862' ; H-beta core
                else :
              endcase
            endif else begin
              ;; Chromis-W and Chromis-D
              case state1 of
                'wheel00006' : prefilter = '4846' ; H-beta continuum
                else: prefilter = '3950'          ; Ca II HK wideband
              endcase
            endelse
            if n_elements(prefilter) gt 0 then begin
              ;; Not defined for darks
              fxaddpar, header, red_keytab('pref'), prefilter, 'Extracted from state keyword'
            endif
          endif
        endif else begin
          case state of
            'hbeta-core' : prefilter = '4862'
            'hbeta-cont' : prefilter = '4846'
            'cah-core'   : prefilter = '3969'
            'wheel00005' : prefilter = '3950' ; WB: Ca II HK continuum
            'wheel00006' : prefilter = '4846' ; WB: H-beta continuum
            else:
          endcase
          if n_elements(prefilter) gt 0 then begin
            ;; Not defined for darks
            fxaddpar, header, red_keytab('pref'), prefilter, 'Extracted from state keyword'
          endif
        endelse
      endelse
    endif


    ;; Correct OBS_PHDU bug in camera software
    OBS_PHDU = sxpar( header, 'OBS_PHDU', comment = pcomment, count=count )
    if count gt 0 then begin
      if n_elements(pcomment) eq 0 || strtrim(pcomment,2) eq '' then $
         pcomment = ' Observational SOLARNET Header and Data Unit'
      fxaddpar, header, 'OBS_SHDU', OBS_PHDU, pcomment, after = 'OBS_PHDU'
      sxdelpar, header, 'OBS_PHDU'
    endif

    if n_elements(framenumber) ne 0 then begin
      ;; We may want to change or remove some header keywords
      ;; here, like FRAME1, CADENCE, and DATE-END.
    endif
  endif else begin
    ;; EXTENSION header
    elun = fxposit(fname,extension,/readonly,/no_fpack)
    fxhread,elun,header
    free_lun,elun
  endelse

  return, header

end

