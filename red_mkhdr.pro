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
;    img : in, type=uintarr
;
;	The ANA image that goes with the header.
;	If known it simplifies the process of setting up the header.
;
; :History:
; 
;   2016-05-18 : JLF. Created.
;
;-
function red_mkhdr, anahdr $
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
sxaddpar,hdr,'EXPTIME',exptime,' [s]',before='COMMENT'
sxaddpar,hdr,'DATE',(red_strreplace(Ts, ' ', 'T'))[0],' '

return,hdr

end