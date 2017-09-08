; docformat = 'rst'

;+
; Return the keywords in the VAR_KEYS keyword from a FITS header.
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
; :Returns:
; 
;    The VAR_KEYS string value.
; 
; :Params:
; 
;    hdr : in, type=strarr
;  
;       The FITS header.
; 
; 
; :Keywords:
; 
;    extensions : out, optional, type=strarr   
;   
;       The extension names in which the returned eywords can be
;       found. 
; 
; 
; :History:
; 
;   2017-09-05 : MGL. First version.
; 
; 
; 
;-
function red_fits_var_keys, hdr $
                            , extensions = extensions

  var_keys = fxpar(hdr, 'VAR_KEYS', count = Nmatch, comment = var_keys_comment)

  ;; Split into keywords (some of which have extension names in
  ;; front).
  keywords = strsplit(var_keys, ',', /extract)
  Nkeywords = n_elements(keywords)

  ;; We need the corresponding extension names
  extensions = strarr(Nkeywords)
  for ikey = 0, Nkeywords-1 do begin
    ;; Identify extension name for this keyword, either before a
    ;; semi-colon or the same extension as the pervious keyword.
    splt = strsplit(keywords[ikey], ';', /extract)
    case n_elements(splt) of
      1 :                       ; Leave the keyword alone
      2 : begin                 ; Remove the extension name from the keyword string 
        latest_extension = splt[0]
        keywords[ikey] = splt[1]
      end
      else : print, 'Should not happen!'
    endcase
    extensions[ikey] = latest_extension
  endfor                        ; ikey

  return, keywords

  
end
