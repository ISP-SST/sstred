; docformat = 'rst'

;+
; Add a parameter to a FITS header, wrapper around fxaddpar.
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
; 
; 
; :Params:
; 
;     header : in, type=strarr
; 
;        The FITS header. See documentation for fxaddpar.
; 
;     name : in, type=string
;     
;        The name of the added parameter. See documentation for fxaddpar.  
;     
;     value : in
;     
;        Value for the added parameter. See documentation for fxaddpar. 
;     
;     comment : in, optional, type=string
;     
;        A comment. See documentation for fxaddpar. 
; 
; 
; 
; 
; :Keywords:
; 
;     anchor : in, out, optional, type=string
;   
;        The name of the parameter relative to which this parameter
;        will be positioned. Decide on before/after based on logical
;        value of those two keywords. (Don't use COMMENT as anchor!) 
; 
;     force : in, optional, type=boolean
; 
;        Set this to force positioning for existing parameters.
;        Implied by anchor keyword.
; 
; :History:
; 
;    2017-05-30 : MGL. First version.
; 
;    2017-07-17 : MGL. Make robust to name and value being unit length
;                 arrays. 
; 
; 
; 
; 
;-
pro red_fitsaddpar, header, name, value, comment $
                    , anchor = anchor $
                    , after = after $
                    , before = before $
                    , force = force $
                    , _ref_extra = extra
  
  if n_elements(anchor) ne 0 then begin
    ;; If an anchor is given, this implies /force.
    force = 1
    case 1 of
      ;; If keyword anchor is given and non-nil, interpret
      ;; before/after as logical and place the new keyword
      ;; before/after the previous keyword. (As in fxaddpar, keyword
      ;; after takes precedence.)
      keyword_set(after)  : aft = anchor
      keyword_set(before) : bef = anchor
      else                : aft = anchor
    endcase
  endif else begin
    ;; Without keyword anchor, use the given before/after, if any.
    if keyword_set(after)  then aft = after
    if keyword_set(before) then bef = before
  endelse

  if keyword_set(force) then begin
    ;; Remove existing occurrence of parameter so positioning will
    ;; work.
    if n_elements(comment) eq 0 then begin
      ;; If the comment keyword is not given, use existing comment if any.
      oldvalue = fxpar(header, name[0], comment = oldcomment, count = oldcount)
      if oldcount gt 0 then comment = oldcomment[0]
    endif
    sxdelpar, header, name[0]
  endif
  
  if n_elements(comment) eq 0 then begin
    fxaddpar, header, name[0], value[0] $
              , after = aft, before = bef $
              , _strict_extra = extra
  endif else begin
    ;; Strip heading and trailing spaces from comment, fxaddpar will
    ;; add one space at the beginning.
    fxaddpar, header, name[0], value[0], strtrim(comment, 2) $
              , after = aft, before = bef $
              , _strict_extra = extra
  endelse
  
  ;; Set this to use in next call
  anchor = name[0]

  ;; Print what we did
  if 0 then begin
    case arg_present(anchor) of
      0 : anchstring = ''
      1 : anchstring = ' ; New anchor=' + anchor
    endcase
    case 1 of
      n_elements(aft) : print, 'Added ' + name[0] + ' after '  + aft + anchstring
      n_elements(bef) : print, 'Added ' + name[0] + ' before ' + bef + anchstring
      else:  print, 'Added ' + name[0] + anchstring
    endcase
  endif
  
end

;; Ideas:
;;
;; * Allow BEFORE and AFTER to be arrays of parameter names, to be
;;   tested for existence in order. Use first match.
