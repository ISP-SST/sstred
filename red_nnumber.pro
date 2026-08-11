FUNCTION red_Nnumber, value, ndigits
;+
; NAME:
;       NNUMBER
; PURPOSE:
;       Convert an integer to a string with leading 0s
; CATEGORY:
;       Type conversion
; CALLING SEQUENCE:
;       RESULT = NNUMBER ( INT, DIGITS )
; INPUTS:
;       INT     : Value for string, converted to integer.
;       DIGITS  : Number of max. digits for the string
; OUTPUTS:
;       String of length DIGITS with leading zeroes
; MODIFICATION HISTORY:
;       15-Dec-1992  P.Suetterlin, KIS
;-

sv = strtrim(fix(value), 2)
num = n_elements(sv)
sval = strarr(num)

FOR k = 0, num-1 DO BEGIN
    svl = strlen(sv(k))
    sva = '0'
    FOR i = 2, ndigits DO sva = sva+'0'
    strput, sva, sv(k), ndigits-svl
    sval(k) = sva
ENDFOR
IF num EQ 1 THEN $
  return, sval(0) $
ELSE $
  return,sval
END

