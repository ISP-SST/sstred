; docformat = 'rst'

;+
; Make a fits header out of an ANA header
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
;     A FITS compatible header from ANA header data.
; 
; :Params:
; 
;    anahdr : in, type=string
;
;       The ANA header (from fz file).
;
; :Keywords:
;
;    img : in, type=intarr
;
;	The ANA image that goes with the header. If known it
;       simplifies the process of setting up the header.
;
;    naxisx : in, optional, type=array
;
;       The dimensions of the array.
;
;
;    datatype : in, optional, type=integer
;
;       The data type of the array.
;
; :History:
; 
;   2016-05-18 : JLF. Created.
;
;   2016-05-20 : MGL. Make more SOLARNET compliant. Get also camera
;                and instrument from header.
;
;   2016-05-25 : MGL. Check if Ts and Te were found.
;
;   2016-05-31 : JLF. Start using red_keytab to keep track of SOLARNET 
; 		 keywords.
;
;   2016-08-05 : JLF. Make DATE-BEG and DATE-END ISO-8601 compliant.
;
;   2016-08-23 : MGL. Deal with ANA headers lacking data size
;                information.
;
;   2016-08-23 : THI. Rename camtag to detector and channel to camera,
;                so the names match those of the corresponding SolarNet
;                keywords.
;
;   2016-09-15 : MGL. New keywords naxisx and datatype. Get what
;                little extra info there is from momfbd program
;                output. 
;
;    2017-03-13 : MGL. Recognize TRIPPEL data.
; 
;    2017-06-01 : MGL. Use red_fitsaddpar and red_mkhdr.
;
;    2017-07-03 : THI. Also try to get camera from the ANA header.
;
;    2017-09-07 : MGL. Changed red_fitsaddpar --> red_fitsaddkeyword. 
;
;-
function red_anahdr2fits, anahdr $
                          , img = img $
                          , naxisx = naxisx $
                          , datatype = datatype
  
  if n_elements(datatype) eq 0 then datatype = 2 ; default?

  if n_elements(img) ne 0 then begin
    red_mkhdr, hdr, img 
  endif else if n_elements(naxisx) ne 0 then begin
    red_mkhdr, hdr, datatype, naxisx
  endif else begin
    posw = strpos(anahdr, ' W=')
    if posw eq -1 then NAXIS1 = 0 else NAXIS1 = strtrim(long(strmid(anahdr, posw+3)), 2)
    posh = strpos(anahdr, ' H=')
    if posh eq -1 then NAXIS2 = 0 else NAXIS2 = strtrim(long(strmid(anahdr, posh+3)), 2)
    red_mkhdr, hdr, 2, [naxis1,naxis2]
  endelse

  anchor = 'DATE'
  
  ;; Header date (as in the FZ output from Michiel's momfbd program)
  dpos = strpos(anahdr, 'DATE=')
  if dpos ne -1 then begin
    date = (red_strreplace(strmid(anahdr, dpos+5, 19), ' ', 'T'))[0]
    red_fitsaddkeyword, anchor = anchor, hdr, 'DATE', date, ''
  endif
  
  ;; Time info (as in CRISP raw data)
  tspos = strpos(anahdr, 'Ts=')
  tepos = strpos(anahdr, 'Te=')
  if tspos ne -1 and tepos ne -1 then begin
    Ts = strmid(anahdr, tspos+3, 26)
    Te = strmid(anahdr, tepos+3, 26)
    red_fitsaddkeyword, anchor = anchor, hdr $
                        ,'DATE-BEG' $
                        , (red_strreplace((red_strreplace(Ts, ' ', 'T')),'.','-',n=2))[0],' '
    red_fitsaddkeyword, anchor = anchor, hdr $
                        ,'DATE-END' $
                        , (red_strreplace((red_strreplace(Te, ' ', 'T')),'.','-',n=2))[0],' '
    exptime = red_time2double(strmid(Te, strpos(Te, ' '))) $
              - red_time2double(strmid(Ts, strpos(Ts, ' ')))
    red_fitsaddkeyword, anchor = anchor, hdr $
                        , 'XPOSURE', exptime, '[s]'
  end

  ;; Camera (as in CRISP raw data)
  campos = strpos(anahdr, '"Camera')
  if campos ne -1 then begin
    detector = 'cam' + (strsplit(strmid(anahdr, campos+8), ' ', /extract))[0]
    red_fitsaddkeyword, anchor = anchor, hdr, red_keytab('detector'), detector, 'Camera identifier'
  end

  ;; Instrument (as in CRISP raw data)
  case 1 of
    strmatch(anahdr,'*"Spectrograph*"*') : instrument = 'TRIPPEL'
    strmatch(anahdr,'*CRISP-*') :          instrument = 'CRISP'
    else:
  endcase
  if n_elements(instrument) gt 0 then red_fitsaddkeyword, anchor = anchor, hdr $
     , 'INSTRUME', instrument, ' Name of instrument'
;  if ipos ne -1 then begin
;     instrument = 'Crisp-'+(strsplit(strmid(anahdr, ipos+6), ']', /extract))[0]
;     red_fitsaddkeyword, hdr, 'INSTRUME', instrument, 'Name of instrument'
;  end
  
  if strmatch(anahdr,'*CRISP-W*',/FOLD) then red_fitsaddkeyword, anchor = anchor, hdr, 'CAMERA', 'Crisp-W' $
  else if strmatch(anahdr,'*CRISP-R*',/FOLD) then red_fitsaddkeyword, anchor = anchor, hdr, 'CAMERA', 'Crisp-R' $
  else if strmatch(anahdr,'*CRISP-T*',/FOLD) then red_fitsaddkeyword, anchor = anchor, hdr, 'CAMERA', 'Crisp-T' $
  else if strmatch(anahdr,'*CRISP-D*',/FOLD) then red_fitsaddkeyword, anchor = anchor, hdr, 'CAMERA', 'Crisp-D'

  ;; Observations date (as in the FZ output from Michiel's momfbd program)
  dpos = strpos(anahdr, 'DATE_OBS')
  if dpos ne -1 then begin
    date_obs = strmid(anahdr, dpos+9, 10)
    ;; Would like to add a time but TIME_OBS is usually empty
    tpos = strpos(anahdr, 'TIME_OBS')
    if tpos ne -1 then begin
      time_obs = strmid(anahdr, tpos+9,dpos-(tpos+9))
      if strlen(time_obs) gt 1 then date_obs += 'T' + time_obs
    endif 
    red_fitsaddkeyword, anchor = anchor, hdr, 'DATE-AVG', date_obs, '', after = 'DATE'
  endif

  ;; Should extract more info from anahdr: states of prefilter, liquid
  ;; crystals, LRE, and HRE. But first find out what keywords to use
  ;; for them in the FITS header.
  
  ;; Add SOLARNET keyword
  red_fitsaddkeyword, hdr, 'SOLARNET', 0.5,  format = 'f3.1' $
                      , 'Fully SOLARNET-compliant=1.0, partially=0.5', before = 'DATE'

  return,hdr

end
