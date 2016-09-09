; docformat = 'rst'

;+
; Read headers from data files in a camera independent way.
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;     J. Lewis Fox, ISP, 2016-05-18
; 
; 
; :Returns:
; 
;     A FITS compatible header (from a data file of any instrument)
; 
; :Params:
; 
;    fname : in, type=string
;
;       The data file.
;
; :Keywords:
;    
;    filetype : in, type=string
;
;	The type of file to read. Allowed values are: 
;		fits
;		fz
;               momfbd
;	If not set auto-detection will be attempted.
;
;    framenumber : in: optional, type=integer
;
;       Specify to get header corresponding to a single frame in a
;       multi-frame image file.
;
;    structheader : in, type=flag
;
;	If set output the header as a structure.
;
;    status : out, type=signed int
;
;	Output the return status, 0 for success, -1 for failure.
;
; :History:
; 
;   2016-05-18 : JLF. Created.
;
;   2016-05-31 : JLF. Added keyword silent to suppress informational
;                messages. 
;
;   2016-08-12 : MGL. Make header for momfbd-format file. New keyword
;                framenumber. Call red_filterchromisheaders only for
;                FITS files without SOLARNET header keyword. Call new
;                function red_meta2head.
;
;   2016-08-19 : MGL. Only test for .momfbd extension if filetype is
;                not specified.      
;
;   2016-08-23 : MGL. Deal with ANA headers lacking data size
;                information.      
;
;   2016-08-31 : MGL. Use red_detectorname instead of red_getcamtag.
;                Add rudimentary handling of tabulated DATE-BEG.       
;
;   2016-09-01 : MGL. Change FRAME1 to FRAMENUM if needed.       
;
;   2016-09-08 : MGL. Allow for headers without prefilter info (like
;                darks). 
;
;
;-
function red_readhead, fname, $
                       filetype = filetype, $
                       framenumber = framenumber, $
                       structheader = structheader, $
                       status = status, $
		       silent = silent
		       
    if( file_test(fname) eq 0 ) then begin
        message, 'File does not exist: '+fname,/info
        status = -1
        return, 0B
    endif

    ;; Remove this when rdx_filetype can recognize .momfbd files:
    if( n_elements(filetype) eq 0 ) then begin
    
        if file_basename(fname,'.momfbd') ne file_basename(fname) then begin
           filetype = 'momfbd'
        endif else begin
           filetype = rdx_filetype(fname)
        endelse
       
        if( filetype eq '' ) then begin
            message, 'Cannot detect filetype. Pass it manually as', /info
            message, "head = red_readhead('"+fname+"',filetype=fits')", /info
            status = -1
            return, 0B
        endif
        
    endif                         ; filetype

        
    case strupcase(filetype) of

        'ANA' : begin
            ;; Data stored in ANA files.
            anaheader = fzhead(fname)
            if n_elements(anaheader) ne 0 then begin
               if strmatch( anaheader, "SIMPLE*" ) gt 0 then begin
                   ;; it's actually a fits-header, split into strarr with length 80
                   len = strlen(anaheader)
                   for i=0,len-1, 80 do begin
                       card = strmid(anaheader,i,80)
                       red_append, tmpheader, card
                       if strmid( card, 0, 2 ) eq 'END' then break
                   endfor
                   header = tmpheader
               endif else begin
                  ;; Does the ana header have size information?
                  posw = strpos(anaheader, ' W=')
                  posh = strpos(anaheader, ' H=')
                  if posw eq -1 or posh eq -1 then begin
                     ;;Get it from the data instead.
                     fzread, data, fname
                     header = red_anahdr2fits(anaheader, img = data)
                  endif else begin
                     ;; Convert ana header to fits header
                     header = red_anahdr2fits(anaheader)
                  endelse
               endelse
            endif
         end

        'FITS' : begin
            ;; Data stored in fits files
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

                  tab = readfits(fname, theader, /exten, /silent)
                  date_beg = ftget(theader,tab, 'DATE-BEG') 
                  
                  fxaddpar, header, 'DATE-BEG', date_beg[0] $
                            , 'First in table.', after = 'DATE'
                  
                  date_end_split = strsplit(date_beg[-1], 'T', /extract)
                  date_end_split[1] = red_time2double(red_time2double(date_end_split[1]) $
                                                      + sxpar(header, 'NAXIS3')*sxpar(header, 'CADENCE') $
                                                      + sxpar(header, 'XPOSURE') $
                                                      , /inv)
                  date_end = strjoin(date_end_split, 'T')
                  
                  fxaddpar, header, 'DATE-END', date_end $
                            , 'Last in table + NAXIS3*CADENCE + XPOSURE.' $
                            , after = 'DATE-BEG'

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
                         'hbeta-core' : prefilter = 'Hb-core'
                         'hbeta-cont' : prefilter = 'Hb-cont'
                         'cah-core'   : prefilter = 'CaH-core'
                         else:
                      endcase
                      fxaddpar, header, red_keytab('pref'), prefilter, 'Extracted from file name'
                   endif
               endif else begin         ; STATE keyword exists but not FILTER1 (e.g. 2016.08.31 data)
                   state_split = strsplit( state, '_',  /extr )
                   if n_elements(state_split) gt 0 then begin
                       state1 = state_split[0]                      ;  for 2016.08.31:  state = 'wheel00002_hrz32600'
                       camera = fxpar(header, red_keytab('camera'), count=count)
                       if count gt 0 then begin
                          if camera eq 'Chromis-N' then begin

                             ;; Chromis-N
                             case state1 of
                                'wheel00001' : prefilter = 'CaK-blue'
                                'wheel00002' : prefilter = 'CaK-core'
                                'wheel00003' : prefilter = 'CaH-core'
                                'wheel00004' : prefilter = 'CaH-red'
                                'wheel00005' : prefilter = 'CaH-cont'
                                'wheel00006' : prefilter = 'Hb-core'
                                else :
                             endcase
                          endif else begin
                             ;; Chromis-W and Chromis-D
                             case state1 of
                                'wheel00006' : prefilter = 'Hb-cont'
                                else: prefilter = 'CaHK-cont'
                             endcase
                          endelse
                          if n_elements(prefilter) gt 0 then begin
                             ;; Not defined for darks
                             fxaddpar, header, red_keytab('pref'), prefilter, 'Extracted from state keyword'
                          endif
                       endif
                   endif else begin
                      case state of
                         'hbeta-core' : prefilter = 'Hb-core'
                         'hbeta-cont' : prefilter = 'Hb-cont'
                         'cah-core'   : prefilter = 'CaH-core'
                         else:
                      endcase
                   endelse
               endelse
            endif

            if n_elements(framenumber) ne 0 then begin
               ;; We may want to change or remove some header keywords
               ;; here, like FRAME1, CADENCE, and DATE-END.
            endif

         end

        'MOMFBD' : begin
           mr = momfbd_read(fname, /names) ; Use /names to avoid reading the data parts
           mkhdr, header, 4, [mr.clip[0,1,1]-mr.clip[0,1,0]+1 $
                              , mr.clip[0,0,1]-mr.clip[0,0,0]+1]
           header = header[where(header ne replicate(' ',80))] ; Remove blank lines
           date_ave = mr.date + 'T' + mr.time
           sxaddpar, header, 'DATE-AVE', date_ave, ' ', before='COMMENT'
        end

    endcase

    ;; Keyword FRAME1 changed to FRAMNUM. Rewrite old headers
    ;; to match the new standard.
    frnm = sxpar(header, 'FRAMENUM', count = Npar)
    if Npar eq 0 then begin
       fr1 = sxpar(header, 'FRAME1', count = Npar)
       if Npar eq 1 then sxaddpar, header, 'FRAMENUM', fr1
    endif


    header = header[where(header ne replicate(' ',80))] ; Remove blank lines

    header = red_meta2head(header, meta={filename:fname})
    
    if keyword_set(structheader) then begin
       header = red_paramstostruct(header)
    endif

    status = 0

    return, header

end
