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

;; The translation table.
;; Add your own translations here.
keyword_translator_table = [['cam',		'CAMERA'  ],$
			    ['camtag',		'CAMERA'  ],$
			    ['frame',		'FRAME1'  ],$
			    ['framenumber',	'FRAME1'  ],$
			    ['pref',		'FILTER1' ],$
			    ['prefilter',	'FILTER1' ],$
			    ['camdir',		'INSTRUME'],$
			    ['cam_channel',	'INSTRUME'],$
			    ['scan',		'SCANNUM' ],$
			    ['scannumber',	'SCANNUM' ]]

for i = 0, nstr-1 do begin
  
  test = strlowcase(metaname[i])
  
  idx = where(keyword_translator_table[0,*] eq test,cnt)
  
  if cnt eq 1 then fitsname[i] = keyword_translator_table[1,idx] $
  else message,'Translation error, no unambiguous keyword translation found.'
  
endfor

return,fitsname

end

print,red_keytab('cam')
print,red_keytab('Camtag')
print,red_keytab(['frame','PREF','scannumber'])
print,red_keytab('cam_channel')
print,red_keytab('camp')

end