; docformat = 'rst'

;+
; Expand a comma and dash separated list of numbers given as a string
; into an integer array, e.g. '1,3-5,7-9,11' --> [1,3,4,5,7,8,9,11].
; 
; :Categories:
;
;    SST observations
; 
; 
; :author:
; 
;    Mats Löfdahl, 2004-03-25 (originally in ANA)
; 
; 
; :returns:
; 
;    An integer array.
; 
; :Params:
; 
;    st_in : in, type=str
; 
;      This is the comma and dash separated list.
; 
; :history:
; 
;     2013-10-03 : MGL. Renamed to the red_ namespace.
; 
; 
; 
;-
function red_expandrange,st_in
  
  st=st_in  ; Non-destructive
  
  ; Remove non-numbers from beginning and end of string
  st = byte(st)
  while st(0) lt byte('0') or st(0) gt byte('9') do st = st(1:*)
  st = reverse(st) ; Reverse
  while st(0) lt byte('0') or st(0) gt byte('9') do st = st(1:*)
  st = string(reverse(st)) ; Reverse string

  n_commas = strcount(st,',')
  n_dashes = strcount(st,'-')
  n_sep = n_commas+n_dashes
  
  st=st+',-'
  
  range = [long(st)]
  for i=0,n_sep-1 do begin
    if strpos(st,',') lt strpos(st,'-') then begin
      ; Comma: Just append the number
      st=strskp(st,',')
      range = [range,long(st)]
    end else begin
      ; Dash, append range
      st=strskp(st,'-')
      first_number = range((size(range,/dimensions))[0]-1)
      last_number  = long(st)
      new_range = IndGen(last_number-first_number) + long(first_number+1)
      range = [range,new_range]
    end
  end
  
  return, range

end;
