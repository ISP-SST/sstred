; docformat = 'rst'

;+
; Convert an array of numbers into a string of comma and dash
; separated (sub-)ranges, e.g., [1,2,3,4,5,8,9,10,14,17,20,21,22] ->
; '[1-5,8-10,14,17,20-22]'
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, ISP
; 
; 
; :Returns:
; 
; 
; :Params:
; 
;    arr : in, type=array
;
;       The array of numbers.
; 
; :Keywords:
; 
;    ld : in, optional, type=string, default='['
;   
;       The delimiter on the left hand side of the resulting string. 
; 
;    rd : in, optional, type=string, default=']'
;   
;       The delimiter on the right hand side of the resulting string. 
; 
; 
; :History:
; 
;     2003-06-19 : MGL. First version (in ANA)
;
;     2009-06-11 : MGL. Ported to IDL.
;
;     2013-10-03 : MGL. Renamed to the red_ namespace.
; 
;     2017-06-29 : MGL. Use rdx_ints2str, will create colon syntax
;                  with increments.
;
;-
function red_collapserange, arr, ld = ld, rd = rd
  
  if n_elements(ld) eq 0 then ld = '['
  if n_elements(rd) eq 0 then rd = ']'

  return, ld + rdx_ints2str(arr) + rd

  ;-----------------------
  
  ;; Simple cases, one or two elements
  if n_elements(arr) eq 1 then return, ld + strtrim(arr[0], 2) + rd
  if n_elements(arr) eq 2 then return, ld + strjoin(string(arr, format = '(i0)'), ',') + rd
  
  strng = ld+strcompress(string(arr[0]), /rem)
  
  for i = 1, n_elements(arr)-2 do begin
    if arr[i] eq arr[i-1]+1 then begin
      if arr[i] ne arr[i+1]-1 then begin
                                ;strng += '!'
        if ~strmatch(strng, '*[,-]') then strng += ','	
        strng += strcompress(string(arr(i)), /rem) + ','
      end else begin
        if ~strmatch(strng, '*-') then strng += '-'
      endelse
    end else begin
      if ~strmatch(strng, '*,') then strng += ','	
      strng += strcompress(string(arr[i]), /rem)
    endelse
  endfor
  
  ;;if (last(strng) NE ',') AND (last(strng) NE '-') THEN strng += ','	
  if ~strmatch(strng, '*[,-]') then strng += ',' 

  ;;return, strng+strcompress(string(last(arr)), /rem)+rd
  return, strng + strcompress(string(arr[n_elements(arr)-1]), /rem) + rd
  
end                             ;
