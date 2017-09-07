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
;     name : in, type=strarr
;     
;        The name of the added parameter/keyword. See documentation
;        for fxaddpar, but here it can be an array.
;     
;     value : in, type=array
;     
;        Value for the added parameter. See documentation for
;        fxaddpar, but here it can be an array. The array notation for
;        multiple keywords means only keywords with the same value
;        type (string, float, integer, ...) can be added in a single
;        call. 
;     
;     comment : in, optional, type=strarr
;     
;        A comment. See documentation for fxaddpar, but here it can be
;        an array.  
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
;    2017-09-05 : MGL. Parameters name, value, comment can be arrays.
;                 Override some of the positioning of fxaddpar by
;                 protecting blank and COMMENT keywords. 
; 
;    2017-09-06 : MGL. Line-wrap blank and COMMENT keywords if needed.
;                 Remove trailing blank lines.
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

  ;; Protect input and set defaults
  if n_elements(name)    eq 0 then names    = '' else names    = name
  if n_elements(value)   eq 0 then values   = '' else values   = value
  if n_elements(comment) eq 0 then comments = '' else comments = comment

  ;; Find the number of keywords to add
  Nkeys = max([n_elements(names), n_elements(values), n_elements(comments)])
  ;; All of names, values, and comments must have the same number of
  ;; elements of be scalar/single-element.
  if max(n_elements(names)    eq [1, Nkeys]) ne 1 then stop
  if max(n_elements(values)   eq [1, Nkeys]) ne 1 then stop
  if max(n_elements(comments) eq [1, Nkeys]) ne 1 then stop
  ;; Make them all arrays of the same length.
  if n_elements(names)    ne Nkeys then names    = replicate(names,    Nkeys)
  if n_elements(values)   ne Nkeys then values   = replicate(values,   Nkeys)
  if n_elements(comments) ne Nkeys then comments = replicate(comments, Nkeys)
  
  ;; Remove trailing blank lines
  Nlines = where(strmatch(header, 'END *'), Nmatch)
  if Nmatch eq 0 then stop
  header = header[0:Nlines]
  
  ;; Protect existing blank, COMMENT and CONTINUE keywords by
  ;; replacing them with keywords of the form +i+, -i-, or *i*, where
  ;; i is a number and the surrounding characters depend on which
  ;; keyword it is that is protected. This will override
  ;; fxaddpar's placement rules.
  namefields = strmid(header,0,8)
  pnames = ['COMMENT ', '        ']
  pchars = ['+',        '-'       ]
  Nchars = n_elements(pchars)
  iprotect = 0
;  pindx = where(namefields eq pnames[0] or $
;                namefields eq pnames[1], Nprotect)
;  for iprotect = 0, Nprotect-1 do begin
;    ichar = where(namefields[pindx[iprotect]] eq pnames)
;    header[pindx[iprotect]] = string(pchars[ichar]+strtrim(iprotect, 2)+pchars[ichar], '(A-8)')+'= ' $
;                              + "'"+strmid(header[pindx[iprotect]], 8)+"'"
;  endfor                        ; iprotect
;stop
  for ikey = 0, Nkeys-1 do begin
    
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

    if keyword_set(force) and max(strtrim(names[ikey],2) eq strtrim(pnames, 2)) eq 0 then begin
      ;; Remove existing occurrence of parameter so positioning will
      ;; work.
      if comments[ikey] eq '' then begin
        ;; If the given comments keyword is empty, use existing comment if any.
        oldvalue = fxpar(header, names[ikey], comment = oldcomment, count = oldcount)
        if oldcount gt 0 then comments[ikey] = oldcomment[0]
      endif
      sxdelpar, header, names[ikey]
    endif

    ;; Special handling of blank and COMMENT keywords.
    match = strtrim(names[ikey],2) eq strtrim(pnames, 2)
    if max(match) eq 0 then begin
      ;; No special handling needed.
      names_ikey = names[ikey]
      ;; Strip heading and trailing spaces from comment, fxaddpar will
      ;; add one space at the beginning.
      fxaddpar, header, names_ikey, values[ikey], strtrim(comments[ikey], 2) $
                , after = aft, before = bef $
                , _strict_extra = extra
      iprotect++
    endif else begin
      ;; Protect keyword name from fxaddpar's positioning.
      ichar = where(match)
      if strlen(values[ikey]) le 72 then begin
        names_ikey = pchars[ichar]+strtrim(iprotect, 2)+pchars[ichar]
        fxaddpar, header, names_ikey, values[ikey] $
                  , after = aft, before = bef $
                  , _strict_extra = extra
        iprotect++
      endif else begin
        ;; Line-wrap the value string
        words = strsplit(values[ikey], /extract)
        rep = 0
        repeat begin
          if rep gt 0 then begin
            aft = names_ikey
            undefine, bef
          endif
          ii = (where(total(strlen(words)+1, /cumulative) gt 68, N68))[0]
          if N68 eq 0 then begin
            wrapped_line = strjoin(words, ' ')
          endif else begin
            wrapped_line = strjoin(words[0:ii-1], ' ')
            words = words[ii:*]
          endelse
          names_ikey = pchars[ichar]+strtrim(iprotect, 2)+pchars[ichar]
          fxaddpar, header, names_ikey, wrapped_line $
                    , after = aft, before = bef $
                    , _strict_extra = extra
          if N68 eq 0 then break ; Exit repeat loop
          anchor = names_ikey
          iprotect++
          rep++
        endrep until 0
      endelse
    endelse

    ;; Set anchor to use for next keyword.
    anchor = names_ikey

    ;; Print what we did
    if 0 then begin
      case arg_present(anchor) of
        0 : anchstring = ''
        1 : anchstring = ' ; New anchor=' + anchor
      endcase
      case 1 of
        n_elements(aft) : print, 'Added ' + names[ikey] + ' after '  + aft + anchstring
        n_elements(bef) : print, 'Added ' + names[ikey] + ' before ' + bef + anchstring
        else:  print, 'Added ' + names[ikey] + anchstring
      endcase
    endif

  endfor                        ; ikey

  ;; Now undo the protection
  firstchars = strmid(header,0,1)
  pindx = where(firstchars eq pchars[0] or $
                firstchars eq pchars[1], Nprotect)
  
  for iprotect = 0, Nprotect-1 do begin
    ichar = where(firstchars[pindx[iprotect]] eq pchars)
    header[pindx[iprotect]] = pnames[ichar] + strmid(header[pindx[iprotect]], 11)
    ;; Strmid from pos 11 removed the first single quote. We need to
    ;; also remove the final quote and any slash that might come
    ;; after.
    pos=strpos(header[pindx[iprotect]],"'",/reverse_search)
    if pos ne -1 then header[pindx[iprotect]] = strmid(header[pindx[iprotect]],0,pos)
    ;; Quotes in the COMMENT line may have been doubled
    header[pindx[iprotect]] = red_strreplace(header[pindx[iprotect]], "''", "'", n = 40)
    ;; Make sure the length is 80 characters
    if strlen(header[pindx[iprotect]]) ne 80 then $
       header[pindx[iprotect]] = strmid(header[pindx[iprotect]]+strjoin(replicate(' ', 80)), 0, 80)
  endfor                        ; iprotect

  ;; Remove trailing empty lines
  header = header[0:where(strmatch(header, 'END *'), Nmatch)]
 
end

;; Ideas:
;;
;; * Allow BEFORE and AFTER to be arrays of parameter names, to be
;;   tested for existence in order. Use first match.
