; docformat = 'rst'

;+
; Correct various mistakes in the fitscube header, that have been made
; during the development of the pipeline. Set the SOLARNET compliance
; status to 1.
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
; :Params:
; 
;    hdr : in, type=string array
; 
;      Fitscube header. 
; 
; 
; :Keywords:
;
;    keywords : in, optional, type=struct
;
;      Keywords and values to be added to the header. Keywords are
;      upcased and checked against a list of SOLARNET-approved
;      keywords.
;
;    no_checksum : in, optional, type=boolean
;
;      Do not calculate DATASUM and CHECKSUM.
;
;    release_date : in, optional, type=string
;
;      The value of the RELEASE keyword, a date after which the data
;      are not proprietary.
;
;    release_comment : in, optional, type=string
;
;      The value of the RELEASEC keyword, a comment with details about
;      the proprietary state of the data, and whom to contact.
;
;
; :History:
;
;   2022-04-05 : OA. Imported the code from 'red__fitscube_finalize'.
; 
;-
pro red::fitscube_header_finalize, hdr $
                            , keywords = keywords $
                            , no_checksum = no_checksum $ ; just for prpara
                            , coordinates = coordinates $
                            , release_date = release_date $
                            , release_comment = release_comment $
                            , feature = feature $
                            , observer = observer $
                            , point_id = point_id $
                            , status = status
  
  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)

  ;; We could look for PRPROCn eq 'red::fitscube_finalize' in the
  ;; header and if we find it show the PRREFn date and ask if we want
  ;; to finalize again.  
  
  ;; Make prpara
  red_make_prpara, prpara, keywords     
  red_make_prpara, prpara, no_checksum 
  red_make_prpara, prpara, release_date      
  red_make_prpara, prpara, release_comment        

  ;; This will be used for a new DATE header.
  new_DATE = red_timestamp(/utc,/iso)

  ;; Delete header keywords that should not be there
  red_fitsdelkeyword, hdr, 'STATE'

  ;; Add keywords after this one:
  red_fitsaddkeyword, anchor = 'EXTNAME', hdr, 'TIMESYS', "UTC"
  anchor = 'TIMESYS'

  ;; Update the DATE keyword
  red_fitsaddkeyword, anchor = anchor, hdr, 'DATE', new_DATE

  for itag = 0, n_tags(keywords)-1 do begin
    ;; Add FITS header keywords as requested
  if max(strmatch(approved_keywords, strupcase((tag_names(keywords))[itag]))) then begin
      red_fitsaddkeyword, anchor = anchor, hdr $
                          , strupcase((tag_names(keywords))[itag]), keywords.(itag)
    endif else begin
      print, inam + ' : The tag '+(tag_names(keywords))[itag]+' is not in the list of approved keywords: ' 
      print, approved_keywords
      status = 0
      return
    endelse
  endfor                        ; itag
  
  ;; Proprietary data?
  if n_elements(release_date) eq 0 then begin
    s = ''
    read, inam + ' : No release date given. Mark as "Not proprietary"? [Yn]', s
    if s eq '' then s = 'Y'
    if strupcase(strmid(s, 0, 1)) eq 'Y' then begin
      release_date = ''
      release_comment = 'Data is not proprietary.'
    endif else begin
      status = 0
      return
    endelse
  endif else begin
    ;; Check that release_date is a proper date
    if ~(strmatch(release_date,'[12][90][0-9][0-9]-[012][0-9]-[0123][0-9]') $
         || release_date eq '') then begin
      print, inam + ' : Release date is not a proper ISO date: '+release_date
      status = 0
      return
    endif
    ;; Proprietary data
    if n_elements(release_comment) eq 0 then begin
      releasec = 'Proprietary data, see RELEASE keyword for release date.'
    endif
    red_fitsaddkeyword, anchor = anchor, hdr, 'RELEASE',  release_date
    red_fitsaddkeyword, anchor = anchor, hdr, 'RELEASEC', release_comment
  endelse

  prev_feature = fxpar(hdr,'FEATURE')
  is_prev_feature = 0B
  if prev_feature then $
     if strtrim(prev_feature,2) ne '' then is_prev_feature = 1B 
  ans=''
  download_ok = red_geturl('https://dubshen.astro.su.se/sst_tags/features.txt',content=approved_features)
  case 1 of
    n_elements(feature) gt 0 and is_prev_feature : begin
      print, 'There is FEATURE keyword in the header: ', prev_feature       
      read,'Would you like to overwrite it? (Y/N): ',ans
      if strupcase(ans) eq 'Y' then begin
        if download_ok then begin
          ftrs = strsplit(feature,',',/extract)
          for j=0, n_elements(ftrs)-1 do begin
            if ~max(strmatch(approved_features, strtrim(ftrs[j],2))) then begin
              print, ftrs[j],' is not in the list of approved FEATURE values: '
              for j=0, n_elements(approved_features)-1 do $
                print,'[',strtrim(string(j),2),'] ', approved_features[j]
              print,'Please provide correct value and try again.'
              status = 0
              return
            endif
          endfor
        endif else begin
          print, "Failed to download list of FEATURE values. Can't verify supplied value."
          status = 0
          return
        endelse
        red_fitsaddkeyword, anchor = anchor, hdr, 'FEATURE', feature
     endif else undefine, feature     
    end
    n_elements(feature) gt 0 and ~is_prev_feature : begin
      if download_ok then begin
        ftrs = strsplit(feature,',',/extract)
        for j=0, n_elements(ftrs)-1 do begin
          if ~max(strmatch(approved_features, strtrim(ftrs[j],2))) then begin
            print
            print, ftrs[j],' is not in the list of approved FEATURE values:'            
            for j=0, n_elements(approved_features)-1 do $
              print,'[',strtrim(string(j),2),'] ', approved_features[j]
            print,'Please provide correct value and try again.'
            status = 0
            return          
          endif
        endfor
      endif else begin
        print, "Failed to download list of FEATURE values. Can't verify supplied value."
        status = 0
        return
      endelse
      red_fitsaddkeyword, anchor = anchor, hdr, 'FEATURE', feature
    end
    n_elements(feature) eq 0 and ~is_prev_feature : begin        
      if download_ok then begin
        print, 'You have to set FEATURE keyword. Please choose from the list: '
        for j=0, n_elements(approved_features)-1 do $
          print,'[',strtrim(string(j),2),'] ', approved_features[j]
        read,'Like "1,4-6": ',ans
        if ans eq '' then begin
          red_fitsaddkeyword, anchor = anchor, hdr, 'FEATURE', 'Missing'  
        endif else begin
          ind = red_expandrange(ans)
          feature = ''
          for j=0, n_elements(ind)-1 do $
             feature += approved_features[ind[j]] + ', '
          feature = strmid(feature,0,strlen(feature)-2)
          red_fitsaddkeyword, anchor = anchor, hdr, 'FEATURE', feature
        endelse
      endif else begin
        print, "Failed to download list of FEATURE values.  Can't print it."
        status = 0
        return
      endelse  
    end
    else : 
  endcase

  prev_point_id = fxpar(hdr,'POINT_ID')
  is_prev_point_id = 0B
  if prev_point_id then $
     if strtrim(prev_point_id,2) then is_prev_point_id = 1B
  case 1 of 
    n_elements(point_id) gt 0 and is_prev_point_id : begin
      print, 'There is POINT_ID keyword in the header: ', prev_point_id        
      read,'Would you like to overwrite it? (Y/N): ',ans
      if strupcase(ans) eq 'Y' then $
        red_fitsaddkeyword, anchor = anchor, hdr, 'POINT_ID', point_id $
      else $
        undefine, point_id
    end
    n_elements(point_id) gt 0 and ~is_prev_point_id : begin        
      red_fitsaddkeyword, anchor = anchor, hdr, 'POINT_ID', point_id
    end
    n_elements(point_id) eq 0 and ~is_prev_point_id : begin
      point_id = fxpar(hdr,'DATE-OBS')
      red_fitsaddkeyword, anchor = anchor, hdr, 'POINT_ID', point_id
    end
    else : 
  endcase

  prev_observer = fxpar(hdr,'OBSERVER')
  is_prev_observer = 0B
  if prev_observer then $
     if strtrim(prev_observer,2) ne '' then is_prev_observer = 1B
  case 1 of 
    n_elements(observer) gt 0 and is_prev_observer : begin
      print, 'There is OBSERVER keyword in the header: ', prev_observer        
      read,'Would you like to overwrite it? (Y/N): ',ans
      if strupcase(ans) eq 'Y' then $
        red_fitsaddkeyword, anchor = anchor, hdr, 'OBSERVER', observer $
      else $
        undefine, observer
    end
    n_elements(observer) gt 0 and ~is_prev_observer : begin        
      red_fitsaddkeyword, anchor = anchor, hdr, 'OBSERVER', observer
    end
    n_elements(observer) eq 0 and ~is_prev_observer : begin
      read, "You have to enter observers' names:", ans
      observer = strtrim(ans,2)
      red_fitsaddkeyword, anchor = anchor, hdr, 'OBSERVER', ans
    end
    else : 
  endcase

  ;; Add FITS keywords with info available in the file but not as
  ;; headers. 

  ;; WCS info: DATE-BEG and DATE-END take care of the temporal
  ;; coordinates, WAVEMIN and WAVEMAX of the spectral coordinates. But
  ;; we need to specify the pointing as an additional WCS coordinate.
  if ~keyword_set(coordinates) then begin
    filename = fxpar(hdr,'FILENAME')
    if file_dirname(filename) eq '.' then begin
      point_id = fxpar(hdr,'POINT_ID')
      if strmatch(point_id, '*grouped*') or strmatch(point_id, '*mosaic*') then $
        filename = 'cubes_concatenated/' + filename $
      else begin
        if strmid(filename,0,2) eq 'nb' then $
          filename = 'cubes_nb/' + filename $
        else $
          filename = 'cubes_wb/' + filename
      endelse
    endif
    red_fitscube_getwcs, filename, coordinates = coordinates
  endif
  red_fitsaddkeyword, anchor = anchor, hdr $
                      , 'WCSNAMEA', 'AVERAGED APPROXIMATE HPLN-TAN/HPLT-TAN CENTER POINT'
  red_fitsaddkeyword, anchor = anchor, hdr $
                      , 'CRPIX1A', fxpar(hdr,'NAXIS1')/2. + 0.5, 'Center pixel of image array'
  red_fitsaddkeyword, anchor = anchor, hdr $
                      , 'CRPIX2A', fxpar(hdr,'NAXIS2')/2. + 0.5, 'Center pixel of image array' 
  red_fitsaddkeyword, anchor = anchor, hdr $
                      , 'CRVAL1A', mean(coordinates.hpln), '[arcsec] Coordinates of center of image array'
  red_fitsaddkeyword, anchor = anchor, hdr $
                      , 'CRVAL2A', mean(coordinates.hplt), '[arcsec] Coordinates of center of image array' 
  red_fitsaddkeyword, anchor = anchor, hdr $
                      , 'CDELT1A', 0.0, 'Zero FOV extent'
  red_fitsaddkeyword, anchor = anchor, hdr $
                      , 'CDELT2A', 0.0, 'Zero FOV extent'

  
  ;; After this, this FITS file should be SOLARNET compliant
  red_fitsaddkeyword, anchor = 'DATE', hdr, 'SOLARNET', 1, 'SOLARNET compliant file'

  ;; Add header info about this step
  prstep = 'HEADER-CORRECTION'
  self -> headerinfo_addstep, hdr $
                              , prstep = prstep $
                              , prpara = prpara $
                              , prproc = inam

  status = 1
end
