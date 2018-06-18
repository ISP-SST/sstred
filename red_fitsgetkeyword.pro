; docformat = 'rst'

;+
; Read a keyword from a FITS header.
;
; Like fxpar but also checks for SOLARNET variable-keywords. It also
; reads record-valued keywords, see the Calabretta et al. draft
; from 2004, "Representations of distortions in FITS world coordinate
; systems", and HIERARCH keywords, see  Wicenec, A., Grosbol,
; P., & Pence, W. 2009, "The ESO HIERARCH Keyword Conventions",
; https://fits.gsfc.nasa.gov/registry/hierarch/hierarch.pdf.
;
; In case the same keyword was written with more than one of these
; mechanisms, the returned values will be taken from HIERARCH with the
; highest priority, record-valued second, and normal keyword last
; (possibly SOLARNET-variable).
;
; * HIERARCH keywords are returned as the same kind of array of lists
;   used as input by red_fitsaddkeyword_hierarch.
;
; * Record-valued keywords are returned as a strarr unless the
;   field_specifiers keyword is present, in which case the returned
;   value is an array of values and the field specifiers are returned
;   in the keyword.
;
; * SOLARNET variable keywords are returned as just the
;   average/min/max value in the "normal" keyword with the same name.
;   The variable info is returned in keywords is specified. Which
;   keywords depends on the type of variable keyword.
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
; :Returns:
;
;    The keyword value.
; 
; :Params:
; 
;     filename_or_header : in, type="string or strarr"
; 
;        If strarr: The FITS header. See documentation for fxpar.
;        If string: The file name. Must be the file name to get
;        variable keywords, unless the file in in the current
;        directory. 
; 
;     name : in, type=string
;     
;        The name of the keyword. See documentation for fxpar.  
;     
; 
; :Keywords:
; 
;     comment : out, optional, type=string
;     
;        A comment. See documentation for fxaddpar. 
; 
;     dir : in, optional, type=string, default="./"
;
;        The directory in which the file can be found. Only needed
;        when a variable-keyword needs to be read and the filename is
;        not given as the first parameter.
;     
;     field_specifiers : out, optional, type=strarr
;     
;        If this keyword is present and the "name" is a record-valued
;        header keyword, the returned values are stripped of the
;        field-specifiers, which are instead returned here.
;     
;     variable_values : out, optional, type=struct
;
;        The variable-keyword values in the form of a struct. The
;        struct depends on the association, but it should always have
;        the tag names "association", "name", and "values". Additional
;        tags may aid the interpretation of the values.
; 
; :History:
; 
;    2017-08-30 : MGL. First version.
; 
;    2017-10-26 : MGL. Renamed from red_fitskeyword. Implemented
;                 record-valued keywords.
; 
;    2017-12-07 : MGL. New keyword comment.
; 
;    2018-04-06 : MGL. Read also HIERARCH keywords.
; 
;    2018-06-15 : MGL. Return variable keywords with the full fitscube
;                 dimensions.
; 
;-
function red_fitsgetkeyword, filename_or_header, name $
                             , count = count $
                             , dir = dir $
                             , comment = comment $
                             , coordinate_values = coordinate_values $
                             , coordinate_names = coordinate_names $
                             , ignore_variable = ignore_variable $
                             , field_specifiers = field_specifiers $
                             , variable_values = variable_values $
                             , _ref_extra = extra

  if n_elements(dir) eq 0 then dir = './'
  
  ;; Get the header
  case n_elements(filename_or_header) of
    0: stop
    1: begin
      filename = filename_or_header
      if ~file_test(filename) then begin
        message, 'File does not exist: ' + filename
        stop
      endif
      hdr = headfits(filename)
    end
    else : begin
      hdr = filename_or_header
      maybe_filename = fxpar(hdr, 'FILENAME', count = cnt)
      if cnt gt 0 then filename = dir+maybe_filename
    end
  endcase

  ;; Check for HIERARCH with the wanted name.
  hindx = where(strmatch(hdr, 'HIERARCH ' + name + ' *'), Nmatch)
  if Nmatch gt 0 then begin
    for imatch = 0, Nmatch-1 do begin
      line = red_strreplace(hdr[hindx[imatch]], 'HIERARCH ' + name + ' ', '')
      split_line = strsplit(line,'=',/extract)
      field_spec = strtrim(split_line[0], 2)
      split_line = strsplit(split_line[1], '/', /extract)
      tmp = execute('field_value='+split_line[0]) ; Get the value with the correct type
      if n_elements(split_line) ge 2 then begin
        field_comment = strjoin(split_line[1:*], ' ')
        field_list = list(field_spec, field_value, field_comment)
      endif else begin
        field_list = list(field_spec, field_value)
      endelse
      red_append, hierarch_fields, field_list
    endfor                      ; imatch

    ;; No need to look further if we found HIERARCH keywords!
    return, hierarch_fields

  endif
  
  ;; Check for multiple occurences. 
  value = fxpar(hdr, name, comment = comment, _strict_extra = extra, count = count)

  ;; Missing keyword
  if count eq 0 then return, 0

  if count gt 1 then begin
    ;; If this keyword occurs more than once, it might be a
    ;; record-valued keyword. Return all values as a strarr, also
    ;; return all comments in a strarr. fxpar cannot do this so we
    ;; will have to work around this.
    xhdr = hdr                  ; Protect the header
    namefields = strtrim(strmid(hdr,0,8),2)
    pos = where(strmatch(namefields, strtrim(name, 2)))
    value = strarr(count)
    comment = strarr(count)
    for i = 0, count-1 do begin
      ;; Replace all occurences of the keyword with a temporary,
      ;; numbered name.
      xhdr[pos[i]] = string('TMP'+strtrim(i+1, 2)+'         ',format='(a8)') $
                     + strmid(hdr[pos[i]], 8)
    endfor                      ; i
    value = fxpar(xhdr, 'TMP*', comment = comment, _strict_extra = extra)
    if arg_present(field_specifiers) then begin
      ;; Strip the values of field specifiers and return them in this keyword.
      field_specifiers = strarr(count)
      for i = 0, count-1 do begin
        rec_parts = strsplit(value[i], ':', /extract)
        field_specifiers[i] = strtrim(rec_parts[0], 2)
        value[i]            = strtrim(rec_parts[1], 2)
      endfor                    ; i
    endif
    return, value
  endif
  
  ;; Deal with the SOLARNET variable-keywords mechanism.
  if arg_present(variable_values) then begin

    ;; Is this a variable-keyword?
    var_keys = red_fits_var_keys(hdr, extensions = extensions)

    matches = strmatch(var_keys, name)
    if max(matches) eq 0 then return, value ; This is not a variable-keyword
    if total(matches) gt 1 then stop

    if ~file_test(filename) then begin
      message, 'File does not exist: ' + filename
      stop
    endif
    
    ;; Open the binary extension and read the extension header
    fxbopen, tlun, filename, extensions[where(matches)], bdr
    naxis2 = fxpar(bdr,'NAXIS2')
    if naxis2 ne 1 then stop    ; Not a SOLARNET var-key extension
    tfields = fxpar(bdr,'TFIELDS')
    

    ;; Identify type of variable-keyword
    ;; "static" if there are iCNAn WCS keywords starting with STATIC
    ;; "coordinate" if there are iCTYPn keywords
    ;; "pixel-to-pixel" otherwise
    ctype1 = fxpar(bdr,'1CTYP1',count=Nctype1)
    cname1 = fxpar(bdr,'1CNAM1',count=Ncname1)
    case 1 of
      Nctype1 gt 0 : $
         association = 'coordinate'
      Ncname1 gt 0 && strmid(cname1, 0, 5) eq 'STATIC' : $
         association = 'static'
      else : $
         association = 'pixel-to-pixel'
    endcase
    
    case association of
      'pixel-to-pixel' : begin
        ;; This should do it for single-valued pixel-to-pixel associated
        ;; variable-keywords:
        fxbread, tlun, values, name

        tdim1 = fxpar(bdr, 'TDIM1')
        values = reform(values,(long(strsplit(strmid(tdim1,1,strlen(tdim1)-2),',',/extr)))[1:*])
        
        ;; The variable_values keyword should be a struct
        variable_values = { association:association   $
                            , name:name               $
                            , values:values           $
                          }
        
        fxbclose, tlun
      end
      'coordinate' : begin
        fxbread, tlun, values, name
        Ncoords = tfields-1

        ;; The variable_values keyword should be a struct
        variable_values = { association:association   $
                            , name:name               $
                            , values:values           $
                            , ncoords:Ncoords         $
                          }

        ;; Need to return values for all keywords used, so it is easy
        ;; to copy data from one file to another.

        coordinate_names = strarr(Ncoords)
        coordinate_units = strarr(Ncoords)
                          
        for icoord = 0, Ncoords-1 do begin
          ;; Read coordinate #icoord, assume it is in the same
          ;; extension (can be checked in strtrim(icoord, 2)+'PS1_0'
          colname = strtrim(fxpar(bdr, strtrim(icoord+1, 2)+'PS1_1'), 2)
          coordinate_names[icoord] = (strsplit(colname, '-', /extract))[0]
          coordinate_units[icoord] = strtrim(fxpar(bdr, strtrim(icoord+1, 2)+'CUNI1'), 2)
          fxbread, tlun, coordinate_values, colname
          
          
          ;; Put the coordinate into the struct
          variable_values = create_struct(variable_values $
                                          , 'coordinate'+strtrim(icoord, 2), coordinate_values)
        endfor

        ;; Put the coordinate names into the struct
        variable_values = create_struct(variable_values $
                                        , 'coordinate_names', coordinate_names $
                                        , 'coordinate_units', coordinate_units $
                                       )

        fxbclose, tlun

      end
      else : message, 'Association "'+association+'" not implemented yet.'
    endcase
    
  endif

  return, value                 ; Return the value from the header.
  
end


;; Test implementation of record-valued keywords
mkhdr, hdr, 0
names = ['TEST1', 'TEST2 field.1', 'TEST2 field.2']
values = [14, 3.1, 4.5]
comments = ['Normal keyword', 'Record valued', 'Record valued again']
red_fitsaddkeyword, hdr, names, values, comments
print, hdr
print

a = red_fitsgetkeyword(hdr, 'TEST2')
print, a

b = red_fitsgetkeyword(hdr, 'TEST2', field_specifiers = field_specifiers, count = cnt)
print, b
print, field_specifiers

c = red_fitsgetkeyword(hdr, 'TEST1')

end

