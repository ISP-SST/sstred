; docformat = 'rst'

;+
; Convert from a string representing a roman numeral to an integer or
; from integer to roman numeral.
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author: 
;
;    Mats Löfdahl, ISP, 2016-05-01
; 
; :Returns:
;
;    A roman numeral or its integer representation.
; 
; 
; :Params:
; 
;    x : in, type="string or integer"
;   
;      If x is a string, convert from roman number to integer, else
;      convert from integer to roman number. If a string, the roman
;      numeral is allowed to be prepended by the substring 'cam'.
; 
; 
; :History:
; 
;   2016-05-01 : MGL. First version.
; 
; 
;-
function red_romannumber, x

  forward_function red_romannumber ; Enable recursion (doesn't seem to be necessary?)
  
  if size(x, /tname) eq 'STRING' then begin

     ;; Translate from roman numeral to long integer

    roman = strupcase(red_strreplace(x, 'cam', ''))
     
    if roman eq '' then return, 0                                       
    if strmid(roman,0,1) eq "M"  then return, 1000 + red_romannumber(strmid(roman, 1))
    if strmid(roman,0,2) eq "CM" then return,  900 + red_romannumber(strmid(roman, 2))
    if strmid(roman,0,1) eq "D"  then return,  500 + red_romannumber(strmid(roman, 1))
    if strmid(roman,0,2) eq "CD" then return,  400 + red_romannumber(strmid(roman, 2))
    if strmid(roman,0,1) eq "C"  then return,  100 + red_romannumber(strmid(roman, 1))
    if strmid(roman,0,2) eq "XC" then return,   90 + red_romannumber(strmid(roman, 2))
    if strmid(roman,0,1) eq "L"  then return,   50 + red_romannumber(strmid(roman, 1))
    if strmid(roman,0,2) eq "XL" then return,   40 + red_romannumber(strmid(roman, 2))
    if strmid(roman,0,1) eq "X"  then return,   10 + red_romannumber(strmid(roman, 1))
    if strmid(roman,0,2) eq "IX" then return,    9 + red_romannumber(strmid(roman, 2))
    if strmid(roman,0,1) eq "V"  then return,    5 + red_romannumber(strmid(roman, 1))
    if strmid(roman,0,2) eq "IV" then return,    4 + red_romannumber(strmid(roman, 2))
    if strmid(roman,0,1) eq "I"  then return,    1 + red_romannumber(strmid(roman, 1))
    
  endif else begin

     arabics = [1000,  900, 500,  400, 100,   90,  50,   40,  10,    9,   5,    4,   1 ]
     romans  = ["M",  "CM", "D", "CD", "C", "XC", "L", "XL", "X", "IX", "V", "IV", "I" ]

     roman  = ''
     arabic = x
     FOR i = 0, 12 do begin
        WHILE arabic ge arabics[i] do begin
           roman  = roman  + romans[i]
           arabic = arabic - arabics[i]
        endwhile
     endfor
     return, roman

  endelse

end
