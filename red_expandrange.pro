; docformat = 'rst'

;+
; Expand a comma separated list of subranges given as a string into an
; integer array.
;
; Subranges can be of three kinds:
;    1. Single numbers, expands to the number..
;    2. Dash separated numbers, '10-15' --> [10,11,12,13,14,15].
;    3. Colon separated numbers with increments, '10:2:20' --> [10,12,14,16,18,20].
; 
; :Categories:
;
;    SST observations
; 
; 
; :Author:
; 
;    Mats Löfdahl, ISP
; 
; 
; :Returns:
; 
;    An integer array.
; 
; :Params:
; 
;    st_in : in, type=str
; 
;      This is the comma and dash separated list.
; 
; :History:
; 
;     2004-03-25 : MGL. First version (in ANA).
;
;     2009-06-11 : MGL. Ported to IDL.
;
;     2013-10-03 : MGL. Renamed to the red_ namespace.
; 
;     2014-01-22 : MGL. Adapt to string functions moved to the str_
;                  namespace. 
; 
;     2017-06-29 : MGL. Use rdx_str2ints, will accept colon syntax
;                  with increments.
; 
;-
function red_expandrange, st_in
  
  st=st_in  ; Non-destructive
  
  ; Remove non-numbers from beginning and end of string
  st = byte(st)
  while st(0) lt byte('0') or st(0) gt byte('9') do st = st(1:*)
  st = reverse(st) ; Reverse
  while st(0) lt byte('0') or st(0) gt byte('9') do st = st(1:*)
  st = string(reverse(st)) ; Reverse string

  return, rdx_str2ints(st)

  ;-----------------------
  
  
  n_commas = red_strcount(st,',')
  n_dashes = red_strcount(st,'-')
  n_sep = n_commas+n_dashes
  
  st=st+',-'
  
  range = [long(st)]
  for i=0,n_sep-1 do begin
    if strpos(st,',') lt strpos(st,'-') then begin
      ; Comma: Just append the number
      st=red_strskp(st,',')
      range = [range,long(st)]
    end else begin
      ; Dash, append range
      st=red_strskp(st,'-')
      first_number = range((size(range,/dimensions))[0]-1)
      last_number  = long(st)
      new_range = IndGen(last_number-first_number) + long(first_number+1)
      range = [range,new_range]
    end
  end
  
  return, range

end
