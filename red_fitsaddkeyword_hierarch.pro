; docformat = 'rst'

;+
; Add a HIERARCH keyword to a FITS header.
;
; The routine implements HIERARCH keywords, see Wicenec, A., Grosbol,
; P., & Pence, W. 2009, "The ESO HIERARCH Keyword Conventions",
; https://fits.gsfc.nasa.gov/registry/hierarch/hierarch.pdf. Contrary
; to "normal" keywords, HIERARCH keywords can occur multiple times
; (with different field specifiers) in a header. With
; red_fitsaddkeyword you can only add them in the same call, any old
; occurences are removed. It is the user's responsibility to make sure
; any combination of keyword and field-specifier is only added once in
; this call.
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
;     fields : in, type="array of lists"
;     
;        Fields for the added parameter, the lists should consist of
;        two or three elements: the field name, the value, and an
;        optional comment.
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
;   2018-04-06 : MGL. First version.
; 
;   2019-09-30 : MGL. New keyword nodelete.
;
;-
pro red_fitsaddkeyword_hierarch, hdr, name, fields $
                                 , anchor = anchor $
                                 , after = after $
                                 , before = before $
                                 , force = force $
                                 , nodelete = nodelete $
                                 , _ref_extra = extra

  
  ;; Find the number of hierarch records to add
  Nkeys = n_elements(fields)
 
  ;; Remove trailing blank lines in the input header 
  Nlines = where(strmatch(hdr, 'END *'), Nmatch)
  if Nmatch eq 0 then stop
  hdr = hdr[0:Nlines]

  ;; Remove any existing occurrences of rec_name in the header
  if ~keyword_set(nodelete) then red_fitsdelkeyword, hdr, name, /hierarch_only
  
  pchar = '+'

  rec_hash = hash()

  iprotect = 0

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

;    if keyword_set(force) and max(strtrim(names[ikey],2) eq strtrim(pnames, 2)) eq 0 then begin
;      ;; Remove existing occurrence of parameter so positioning will
;      ;; work.
;      if comments[ikey] eq '' then begin
;        ;; If the given comments keyword is empty, use existing comment if any.
;        oldvalue = fxpar(hdr, names[ikey], comment = oldcomment, count = oldcount)
;        if oldcount gt 0 then comments[ikey] = oldcomment[0]
;      endif
;      sxdelpar, hdr, names[ikey]
;    endif

     
     
    ;; Construct the fields string
    hierarch_field = ' ' + name + ' ' + (fields[ikey])[0] + ' = '
    case size((fields[ikey])[1], /tname) of ; Add the fields
      'STRING' : hierarch_field += "'" + (fields[ikey])[1] +  "'" 
      else : hierarch_field += strtrim((fields[ikey])[1], 2)
    endcase
    if n_elements(fields[ikey]) gt 2 then begin ; Comment?
      if strlen(hierarch_field) lt 22 then begin
        hierarch_field += strjoin(replicate(' ', 22-strlen(hierarch_field)))
      endif
      hierarch_field +=  ' / ' + strtrim((fields[ikey])[2], 2)
    endif
          
    ;; Add the record_valued keyword in a protected form, so as
    ;; to not have it removed if it is added more than once in
    ;; this call.
    names_ikey = pchar + strtrim(iprotect, 2) + pchar

    ;; Add it to the header
    fxaddpar, hdr, names_ikey, hierarch_field $
              , after = aft, before = bef $
              , _strict_extra = extra
    

    iprotect++

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
  firstchars = strmid(hdr,0,1)
  pindx = where(firstchars eq pchar, Nprotect)
  
  for iprotect = 0, Nprotect-1 do begin
    hdr[pindx[iprotect]] = 'HIERARCH' + strmid(hdr[pindx[iprotect]], 11)
    ;; Strmid from pos 11 removed the first single quote. We need to
    ;; also remove the final quote and any slash that might come
    ;; after.
    pos=strpos(hdr[pindx[iprotect]],"'",/reverse_search)
    if pos ne -1 then hdr[pindx[iprotect]] = strmid(hdr[pindx[iprotect]],0,pos)
    ;; Quotes in the COMMENT line may have been doubled
    hdr[pindx[iprotect]] = red_strreplace(hdr[pindx[iprotect]], "''", "'", n = 40)
    ;; Make sure the length is 80 characters
    if strlen(hdr[pindx[iprotect]]) ne 80 then $
       hdr[pindx[iprotect]] = strmid(hdr[pindx[iprotect]]+strjoin(replicate(' ', 80)), 0, 80)
  endfor                        ; iprotect

  ;; Remove trailing empty lines
  hdr = hdr[0:where(strmatch(hdr, 'END *'), Nmatch)]
 
end



mkhdr, hdr, 0
hprint, hdr
print

undefine, hierarch_dw3
red_append, hierarch_dw3, list('EXTVER',      1,      'Extension version number')
red_append, hierarch_dw3, list('NAXES',       3.3,    '3 axes in the lookup table')
red_append, hierarch_dw3, list('AXIS ONE TWO THREE FOUR',    'aaa',  'Spatial X')

anchor = 'DATE'
red_fitsaddkeyword_hierarch, anchor = anchor, hdr, 'DW3', hierarch_dw3

undefine, hierarch_dw3
red_append, hierarch_dw3, list('EXTVER',      2,      'Extension version number')
red_append, hierarch_dw3, list('NAXES',       3.3,    '3 axes in the lookup table')
red_append, hierarch_dw3, list('AXIS TWO THREE FOUR ONE',    'bbb',  'Spatial X')

;anchor = 'DATE'
red_fitsaddkeyword_hierarch, anchor = anchor, hdr, 'DW3', hierarch_dw3, /nodelete



hprint, hdr

end
