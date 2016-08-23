; docformat = 'rst'

;+
; Make a fits header out of and ANA header
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
;	The ANA image that goes with the header.
;	If known it simplifies the process of setting up the header.
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
;-
function red_anahdr2fits, anahdr $
		    , img = img

  if n_elements(img) ne 0 then $
     mkhdr,hdr,img $
  else begin
     posw = strpos(anahdr, ' W=')
     if posw eq -1 then NAXIS1 = 0 else NAXIS1 = strtrim(long(strmid(anahdr, posw+3)), 2)
     posh = strpos(anahdr, ' H=')
     if posh eq -1 then NAXIS2 = 0 else NAXIS2 = strtrim(long(strmid(anahdr, posh+3)), 2)
     mkhdr,hdr,2,[naxis1,naxis2]
  endelse

  ;; Time info
  tspos = strpos(anahdr, 'Ts=')
  tepos = strpos(anahdr, 'Te=')
  if tspos ne -1 and tepos ne -1 then begin
     Ts = strmid(anahdr, tspos+3, 26)
     Te = strmid(anahdr, tepos+3, 26)
     exptime = red_time2double(strmid(Te, strpos(Te, ' '))) $
               - red_time2double(strmid(Ts, strpos(Ts, ' ')))
     sxaddpar,hdr,'XPOSURE',exptime,' [s]', before='COMMENT'
     sxaddpar,hdr,'DATE-BEG',(red_strreplace((red_strreplace(Ts, ' ', 'T')),$
      '.','-',n=2))[0],' ', before='COMMENT'
     sxaddpar,hdr,'DATE-END',(red_strreplace((red_strreplace(Te, ' ', 'T')),$
      '.','-',n=2))[0],' ', after='DATE-BEG'
  end

  ;; Camera
  campos = strpos(anahdr, '"Camera')
  if campos ne -1 then begin
     cam = 'cam' + (strsplit(strmid(anahdr, campos+8), ' ', /extract))[0]
     sxaddpar, hdr, red_keytab('cam'), cam, 'Name of camera'
  end

  ;; Instrument
  ipos = strpos(anahdr, 'CRISP-')
  if ipos ne -1 then begin
     instrument = 'Crisp-'+(strsplit(strmid(anahdr, ipos+6), ']', /extract))[0]
     sxaddpar, hdr, 'INSTRUME', instrument, 'Name of instrument'
  end
  
  ;; Should extract more info from anahdr: states of prefilter, liquid
  ;; crystals, LRE, and HRE. But first find out what keywords to use
  ;; for them in the FITS header.
  
  ;; Add SOLARNET keyword
  sxaddpar, hdr, 'SOLARNET', 0.5,  format = 'f3.1' $
            , 'Fully SOLARNET-compliant=1.0, partially=0.5', before = 'COMMENT'

  return,hdr

end
