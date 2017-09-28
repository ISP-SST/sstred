; docformat = 'rst'

;+
; Read a keyword from a FITS header. Like fxpar but also checks for
; SOLARNET variable-keywords.
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
;        If string: The file name.
; 
; 
; :Keywords:
; 
;     name : in, type=string
;     
;        The name of the keyword. See documentation for fxpar.  
;     
;     value : in
;     
;        Value for the added parameter. See documentation for fxaddpar. 
;     
;     comment : in, optional, type=string
;     
;        A comment. See documentation for fxaddpar. 
; 
;     dir : in, optional, type=string, default="./"
;
;        The directory in which the file can be found. Only needed
;        when a variable-keyword needs to be read and the filename is
;        not given as the first parameter.
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
; 
; 
; 
;-
function red_fitskeyword, filename_or_header, name $
                          , dir = dir $
                          , coordinate_values = coordinate_values $
                          , coordinate_names = coordinate_names $
                          , ignore_variable = ignore_variable $
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
      filename = dir+fxpar(hdr, 'FILENAME')
    end
  endcase

  ;; Read the keyword
  value = fxpar(hdr, name, value, _strict_extra = extra)

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
