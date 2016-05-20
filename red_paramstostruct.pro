; docformat = 'rst'

;+
; 
;   Converts a parameter string array (eg, header for FITS files) into
;   an IDL structure where each keyword is a tag of a structure and
;   each of these structures contain a value and a comment. The values
;   should by datatyped (float, int, string, etc) correctly. Logical T
;   is transformed to 1byte and F to 0byte. Does not work on PDS
;   labels because it requires unique keynames and OBJECT and
;   END_OBJECT can be repeated.
; 
; :Categories:
;
;    Datafile handling
; 
; 
; :Author:
; 
;   Ed Shaya / U. of Maryland [April 27, 2012] (as paramtostruct)
; 
; 
; :Returns:
; 
;    A structure in which each keyword is a tag holding a structure.
;    Each of these structure hold a value and a comment.
;
; :Params:
; 
;   input : in, type=string
;
;      A string consisting of 'key = value'. Value can have spaces.
;      Comments (which describe the keyvalue) are on the same line and
;      after the value and start with the commentChar character. (eg,
;      ' keyword = 17.3 / A keyword with floating point value')
;
;
; :Keywords:
; 
;   commentChar : in, optional, type=string, default='/'
;
;      Character that begins the comment on each line. 
;   
;   
;   inverse : in, optional, type=boolean
; 
;      Set this to convert from the struct form to a string array
;      suitable for making a FITS header.
; 
; :History:
;
;     2016-05-18 : MGL. Change documentation header format. Add
;                  keyword inverse, this will make a string array out
;                  of a struct. Change input variable name to input.
;                  Check if commentChar is given with n_elements, not
;                  keyword_set.
; 
; 
;-
function red_paramstostruct, input, commentChar=commentChar, inverse = inverse

  IF n_elements(commentChar) eq 0 THEN commentChar = '/'

  if keyword_set(inverse) then begin

     ;; Convert a structure into a string array

     tags = tag_names(input)
     Nparams = n_tags(input)

     vallengths = lonarr(Nparams)
     for i = 0, Nparams-1 do begin
        if tags[i] ne 'COMMENT' and tags[i] ne 'HISTORY' then vallengths[i] = strlen(input.(i).value)
     endfor
     vallength = max(vallengths) >30
     commentpos = 10 + vallength + 1

     header = strarr(Nparams)

     for i = 0, Nparams-1 do begin

        line = blanks(80)

        strput, line, tags[i], 0
        case 1 of
           tags[i] eq 'COMMENT' or tags[i] eq 'HISTORY' : strput, line, input.(i).value, 8
           tags[i] eq 'SIMPLE' : begin
              strput, line, '=', 8
              if input.(i).value eq 'T' or input.(i).value eq '1B' then begin
                 strput, line, string('T', format = '(a'+strtrim(vallength, 2)+')'), 10
              endif else begin
                 strput, line, string('F', format = '(a'+strtrim(vallength, 2)+')'), 10
              endelse
           end
           else : begin
              strput, line, '=', 8
              strput, line, string(input.(i).value, format = '(a'+strtrim(vallength, 2)+')'), 10
              if n_tags(input.(i)) gt 1 then begin
                 strput, line, commentChar, commentpos
                 strput, line, strtrim(input.(i).comment, 2), commentpos+2
              endif             ; E comment?
           end
        endcase

        header[i] = line

     endfor                     ; i

     return, header

  endif else begin

     Nparams = N_ELEMENTS(input)
     executestring = ''
     FOR i = 0, Nparams-1 DO BEGIN
        ;; This function should provide 3-element array
        ;; [keyword,value,comment] Remove Carriage return at the end of
        ;; value
        par = input[i]
        pos = STREGEX(par,string(10b)+'$') 
        IF  (pos NE -1) THEN  par = strmid(par,0,pos)
        pos = STREGEX(par,string(13b)+'$') 
        IF  (pos NE -1) THEN  par = strmid(par,0,pos)
                                ; Remove blanks
        par = strtrim(par,2) 
                                ; If it is now an empty string go to next line
        if (par eq '') then continue

        keyval = red_parseparameter(par,commentChar=commentChar)

        ;; Reached the END?
        IF (keyval[0] EQ 'END') THEN BREAK

        ;; Replace invalid IDL characters with underscore
        keyval[0] = IDL_VALIDNAME(keyval[0],/convert_all)

        ;; Here we use IDL's on-the-fly datatyping to handle the datatype
        IF (size(keyval,/dimension) GE 2) THEN BEGIN
           IF (keyval[1] EQ 'T') THEN keyval[1] = '1B'
           IF (keyval[1] EQ 'F') THEN keyval[1] = '0B'
           IF (keyval[1] EQ '') THEN keyval[1] = "''"
;	Result = EXECUTE('keyval1 = '+keyval[1])
;	if (Result eq 0) then print,' Problem with keyword ',keyval[0]
        ENDIF

        ;; Create structures for each keyword with value and comment if
        ;; there.
        CASE size(keyval,/dimension)  OF
           3: Result = EXECUTE(keyval[0]+'_valstruct = {value : keyval[1], comment : keyval[2]}')
           2: Result = EXECUTE(keyval[0]+'_valstruct = {value : keyval[1]}')
           1: Result = EXECUTE(keyval[0]+'_valstruct = {value : ''}')
        ENDCASE

        executestring = executestring+','+keyval[0]+' : '+keyval[0]+'_valstruct'
     ENDFOR
     ;; Remove extra comma at the beginning
     executestring = STRMID(executestring,1)

     ;; Finally create the header structure
     Result = EXECUTE('struct = {'+executestring+'}')

     RETURN, struct
  endelse

END
