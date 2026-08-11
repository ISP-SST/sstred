; docformat = 'rst'

;+
; Add a parameter to a FITS header, wrapper around fxaddpar.
;
; The routine implements record-valued keywords, see the Calabretta et
; al. draft from 2004, "Representations of distortions in FITS world
; coordinate systems". Contrary to "normal" keywords, record-valued
; keywords can occur added multiple times (with different field
; specifiers) in a header. But with red_fitsaddkeyword you can only
; add them in the same call, any old occurences are removed. It is the
; user's responsibility to make sure any combination of keyword and
; field-specifier is only added once in this call.
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
;        for fxaddpar, but here it can be an array. If name has two
;        parts, separated by a space character, the first part is
;        interpreted as the name of a record-valued keyword with the
;        second part as the field-specifier.
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
; :Keywords:
; 
;     anchor : in, out, optional, type=string
;   
;        The name of the parameter relative to which this parameter
;        will be positioned. Decide on before/after based on logical
;        value of those two keywords. (Don't use COMMENT or HISTORY as
;        anchor!)
; 
;     force : in, optional, type=boolean
; 
;        Set this to force positioning for existing keywords. Implied
;        by anchor keyword.
; 
;     nodelete : in, optional, type=boolean
; 
;        The default behavior is to delete existing instances of the
;        same HIERARCH keyword. With /nodelete this is not done.
; 
; 
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
;    2017-09-06 : MGL. Renamed from red_fitsaddpar.
; 
;    2017-10-26 : MGL. Implemented record-valued keywords.
; 
;    2017-11-13 : MGL. Implemented HISTORY keywords.
; 
;    2018-03-27 : MGL. Record-valued keywords, values with zero
;                 decimals are printed as integers.
;
;    2019-09-30 : MGL. New keyword nodelete.
;
;-
pro red_fitsaddkeyword, header, name, value, comment $
                        , anchor = anchor $
                        , after = after $
                        , before = before $
                        , force = force $
                        , nodelete = nodelete $
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
  Nlines = where(strmatch(header, 'END *') or header eq 'END', Nmatch)
  if Nmatch eq 0 then stop
  header = header[0:Nlines]

  ;; Protect existing blank, COMMENT and CONTINUE keywords by
  ;; replacing them with keywords of the form +i+, -i-, or *i*, where
  ;; i is a number and the surrounding characters depend on which
  ;; keyword it is that is protected. This will override
  ;; fxaddpar's placement rules.
  namefields = strmid(header,0,8)
  pnames = ['COMMENT ', '        ', 'HISTORY ']
  pchars = ['+',        '-'       , '@'       ]
  red_append, pchars, '*'       ; One extra protection character for record-valued keywords
;  Nchars = n_elements(pchars)
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
        ;; before/after as boolean and place the new keyword
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
      if strlen(names[ikey]) le 8 then red_fitsdelkeyword, header, names[ikey]
    endif

    ;; Special handling of blank and COMMENT keywords.
    match = (strtrim(names[ikey],2) eq strtrim(pnames, 2))
    ;; Also record-valued keywords have to be handled. They consist of
    ;; two parts separated by a space character. The first part is a
    ;; regular FITS header keyword, the second part (the field
    ;; specifier) is either a string with capital letters A-Z and the
    ;; characters minus and underscore (the field identifier) OR a
    ;; field identifier followed by a period and a field index (a
    ;; number).
    record_regex = '^([-_A-Z0-9]*) ([-_.A-Z0-9]*)'
    case 1 of
      max(match) gt 0: begin
        ;; Protect COMMENT and blank keywords from fxaddpar's
        ;; positioning.
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
      end
      stregex(names[ikey], record_regex, /boolean) : begin
        ;; This is a record-valued keyword

        ichar = n_elements(pchars)-1 ; The last element in pchars is used for this

        ;; Split the record-valued keyword
;        rec_parts     = stregex(names[ikey], record_regex, /subexpr, /extract)
        rec_parts     = strsplit(names[ikey], ' ', /extract)
        rec_name      = rec_parts[0]
        rec_fieldspec = rec_parts[1]

        ;; Construct the value string
        if values[ikey] eq long(values[ikey]) then begin
          rec_value = rec_fieldspec + ': ' + strtrim(long(values[ikey]), 2)          
        endif else begin
          rec_value = rec_fieldspec + ': ' + strtrim(values[ikey], 2)
        endelse
        
        ;; Remove any existing occurrences of rec_name in the header
        if ~keyword_set(nodelete) then red_fitsdelkeyword, header, rec_name
        
        ;; Add the record_valued keyword in a protected form, so as to
        ;; not have it removed if it is added more than once in this
        ;; call.
        names_ikey = pchars[ichar]+strtrim(iprotect, 2)+pchars[ichar]
        
        ;; Store the real keyword name in a hash so we can put it into
        ;; the header during the cleaning phase below.
        if n_elements(rec_hash) eq 0 then rec_hash = hash()
        rec_hash[names_ikey] = rec_name
        
        ;; Strip heading and trailing spaces from comment, fxaddpar will
        ;; add one space at the beginning.
        fxaddpar, header, names_ikey, rec_value, strtrim(comments[ikey], 2) $
                  , after = aft, before = bef $
                  , _strict_extra = extra

        iprotect++
      end
      else : begin
        ;; No special handling of the name needed.
        names_ikey = names[ikey]
        
        ;; If placement is specified, then delete any old occurence of
        ;; the keyword.
        if n_elements(bef) eq 0 and n_elements(aft) eq 0 then $
           red_fitsdelkeyword, header, names_ikey
        
        ;; Strip heading and trailing spaces from comment, fxaddpar will
        ;; add one space at the beginning.
        fxaddpar, header, names_ikey, values[ikey], strtrim(comments[ikey], 2) $
                  , after = aft, before = bef $
                  , _strict_extra = extra
        iprotect++
      end
    endcase

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
  for ich = 0, n_elements(pchars)-2 do begin
    wh = where(firstchars eq pchars[ich], Nwh)
    if Nwh gt 0 then red_append, pindx, wh
  endfor                        ; ich
  Nprotect = n_elements(pindx)

  ;; First COMMENT, HISTORY, and blank keywords
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
  
  ;; Replace any protected record-valued keyword names with the actual
  ;; names.
  pindx = where(firstchars eq pchars[n_elements(pchars)-1], Nprotect)
  for iprotect = 0, Nprotect-1 do begin
    pname = strtrim(strmid(header[pindx[iprotect]], 0, 8), 2)
    header[pindx[iprotect]] = string(rec_hash[pname]+'         ',format='(a8)') $
                              + strmid(header[pindx[iprotect]], 8)
  endfor                        ; iprotect

  
  ;; Remove trailing empty lines
  header = header[0:where(strmatch(header, 'END *'), Nmatch)]
 
end


; Test implementation of record-valued keywords and /nodelete
mkhdr, hdr, 0
hprint, hdr
print


names = ['TEST1', 'TEST2 field.1', 'TEST2 field.2']
values = [14, 3.1, 4]
comments = ['Normal keyword', 'Record valued', 'Record valued again']
red_fitsaddkeyword, hdr, names, values, comments
hprint, hdr
print



names = ['TEST1', 'TEST2 field.1', 'TEST2 field.2']
values = [18, 6.4, 42]
comments = ['Normal keyword', 'Record valued', 'Record valued again']
red_fitsaddkeyword, hdr, names, values, comments, /nodelete
hprint, hdr
print


end

;; Ideas:
;;
;; * Allow BEFORE and AFTER to be arrays of parameter names, to be
;;   tested for existence in order. Use first match.

mkhdr, hdr, 0
print, hdr


names = ['HISTORY', '', 'COMMENT', 'HISTORY']
values = ['Some history that is relevant to the file. It could be long enough to span several lines. Or at least I hope it can.', '', "I'd like to add a comment", 'Another piece of history that is relevant to the file. It could be long enough to span several lines. Or at least I hope it can.']
red_fitsaddkeyword, hdr, names, values, comments
print, hdr
stop



;; Test implementation of record-valued keywords
mkhdr, hdr, 0
hprint, hdr
print

names = ['TEST1', 'TEST2 field.1', 'TEST2 field.2']
values = [14, 3.1, 4]
comments = ['Normal keyword', 'Record valued', 'Record valued again']
red_fitsaddkeyword, hdr, names, values, comments
print, hdr
print

names = ['TEST1', 'TEST2 field.1', 'TEST2 field.7']
values =  [17, 30.1, 40.5]
comments = ['Normal keyword', 'Record valued', 'Record valued again'] + ': new value'
red_fitsaddkeyword, hdr, names, values, comments

hprint, hdr
print

end
