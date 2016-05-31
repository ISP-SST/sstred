; docformat = 'rst'

;+
; One-stop shop for standard SOLARNET FITS keyword names.
;
; Send it our own private (unchanging) name for some meta data
; Get back the (current) appropriate SOLARNET keyword name.
; 
; Extensible, so if you change the magic word in your
; code you add it here and can still get the correct FITS keyword.
;
; If you want to add a keyword to a header and you are not 100% certain
; of the proper SOLARNET name for it, pick your own name and put it in the
; the translator below. Then use it where-ever you want the FITS keyword name
; as red_keytab('mykeyname')
;
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;     J. Lewis Fox, ISP
; 
; 
; :Returns:
; 
;    keynames : out, type=strarr
; 
; :Params:
; 
;    strings : in, type=strarr
;   
;      A list of strings to translate to SOLARNET keywords
;   
; :Keywords:
; 
;
; 
; :History:
; 
;   2016-05-31 : JLF. First version.
; 
function red_keytab,metaname

compile_opt idl2

; did you pass me strings?
if size(metaname,/type) ne 7 then message,'Input must be type string (array)'

nstr = n_elements(metaname)

fitsname = metaname

for i = 0, nstr-1 do begin
  
  test = strlowcase(metaname[i])
  
  case 1 of 
    
    ; camera name
    test eq 'cam' || test eq 'camtag': fitsname[i] = 'CAMERA'
    
    ; frame number
    test eq 'frame' || test eq 'framenumber': fitsname[i] = 'FRAME1'

    ; prefilter
    test eq 'pref' || test eq 'prefilter': fitsname[i] = 'FILTER1'
    
    ; camdir
    test eq 'camdir' || test eq 'cam_channel': fitsname[i] = 'INSTRUME'
    
    ; scan number
    test eq 'scan' || test eq 'scannumber': fitsname[i] = 'SCANNUM'
    
  endcase

endfor

return,fitsname

end

print,red_keytab('cam')
print,red_keytab('Camtag')
print,red_keytab(['frame','PREF','scannumber'])
print,red_keytab('cam_channel')
print,red_keytab(cam)

end