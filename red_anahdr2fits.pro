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
;-
function red_anahdr2fits, anahdr $
		    , img = img

  if n_elements(img) ne 0 then $
     mkhdr,hdr,img $
  else begin
     NAXIS1 = strtrim(long(strmid(anahdr, strpos(anahdr, ' W=')+3)), 2)
     NAXIS2 = strtrim(long(strmid(anahdr, strpos(anahdr, ' H=')+3)), 2)
     mkhdr,hdr,2,[naxis1,naxis2]
  endelse

  ;; Time info
  Ts = strmid(anahdr, strpos(anahdr, 'Ts=')+3, 26)
  Te = strmid(anahdr, strpos(anahdr, 'Te=')+3, 26)
  exptime = red_time2double(strmid(Te, strpos(Te, ' '))) $
            - red_time2double(strmid(Ts, strpos(Ts, ' ')))
  sxaddpar,hdr,'XPOSURE',exptime,' [s]', before='COMMENT'
  sxaddpar,hdr,'DATE-BEG',(red_strreplace(Ts, ' ', 'T'))[0],' ', before='COMMENT'
  sxaddpar,hdr,'DATE-END',(red_strreplace(Te, ' ', 'T'))[0],' ', after='DATE-BEG'

  ;; Camera
  campos = strpos(h, '"Camera')
  if campos ne -1 do begin
     cam = 'cam' + (strsplit(strmid(h, campos+8), ' ', /extract))[0]
     sxaddpar, hdr, 'CAMERA', cam, 'Name of camera'
  end

  ;; Instrument
  ipos = strpos(h, 'CRISP-')
  if ipos ne -1 do begin
     instrument = 'Crisp-'+(strsplit(strmid(h, ipos+6), ']', /extract))[0]
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
