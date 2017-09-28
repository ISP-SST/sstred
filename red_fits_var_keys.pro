; docformat = 'rst'

;+
; Return the keywords in the VAR_KEYS keyword from a FITS header.
; Optionally add a keyword, both to output and to the keyword in the
; header. 
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
;       The extension names in which the returned keywords can be
;       found. 
; 
;    new_keyword : in, optional, type=string
;
;       A keyword to be added.
;
;    new_extension : out, optional, type=string
;
;       The extension of the new keyword, constructed from the
;       new_keyword.
;
;    var_keys : out, optional, type=string
;
;       The VAR_KEYS keyword value (with optional new_keyword added).  
;
; :History:
; 
;   2017-09-05 : MGL. First version.
; 
;   2017-09-08 : MGL. New keywords new_keyword and new_extension. 
; 
; 
;-
function red_fits_var_keys, hdr $
                            , extensions = extensions $
                            , new_keyword = new_keyword $
                            , new_extension = new_extension $
                            , var_keys = var_keys

  ;;  Format : var_keys = 'VAR-EXT-1;KEYWD_1,KEYWD_2[He_I],VAR-EXT-2;KEYWD_3'

  ;; Get the existing var_keys keyword
  var_keys = fxpar(hdr, 'VAR_KEYS', count = Nmatch, comment = var_keys_comment)

  if Nmatch eq 0 then begin
    ;; This header does not have a VAR_KEYS keyword.
    if n_elements(new_keyword) ne 0 then begin
      ;; But we want to add one!
      new_extension = 'VAR-EXT-'+new_keyword ; Variable keyword extension for keyword <keyword_name>.
      keywords   = [ new_keyword   ]
      extensions = [ new_extension ]
      var_keys = new_extension + ';' + new_keyword
      fxaddpar, hdr, 'VAR_KEYS', var_keys, 'SOLARNET variable-keywords', after = 'DATE'
      return, keywords
    endif else begin
      ;; Just return 0
      undefine, extensions
      var_keys = ''
      return, 0
    end
  endif
  
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

  if n_elements(new_keyword) then begin

    ;; Check if it is already there
    match = new_keyword eq keywords
    if max(match) eq 1 then begin
      ;; Already there
      new_extension = extensions[where(match)]
    endif else begin
      ;; Add it
      new_extension = 'VAR-EXT-'+new_keyword ; Variable keyword extension for keyword <keyword_name>.
      match = new_extension eq extensions
      if max(match) eq 1 then begin
        ;; Extension exists already, splice in the new keyword
        pos = (where(match))[0]
        case pos of
          0 : begin
            ;; First
            extensions = [ new_extension, extensions ]
            keywords   = [ new_keyword,   keywords   ]
          end
          else : begin
            extensions = [ extensions[0:pos-1], new_extension, extensions[pos:*] ]
            keywords   = [ keywords[0:pos-1],   new_keyword,   keywords[pos:*]   ]
          end
        endcase
        ;; Rebuild the var_keys string
        var_keys = extensions[0] + ';' + keywords[0]
        for ikey = 1, Nkeywords-1 do begin
          if extensions[ikey] eq extensions[ikey-1] then begin
            var_keys += ',' + keywords[ikey]
          endif else begin
            var_keys += ',' + extensions[ikey] + ';' + keywords[ikey]
          endelse
        endfor                  ; ikey
      endif else begin
        ;; Does not exist, add at the end.
        var_keys += ',' +  new_extension + ';' + new_keyword
        red_append, keywords,   new_keyword
        red_append, extensions, new_extension
      endelse
    endelse

    ;; Put the new var_keys string into the header
    fxaddpar, hdr, 'VAR_KEYS', var_keys, var_keys_comment

  endif
  
  return, keywords

  
end
