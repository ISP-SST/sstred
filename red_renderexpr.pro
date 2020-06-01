; docformat = 'rst'

;+
; Render an expression of the type used with mpfitexpr as a string
; with the parameters inserted.
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
; :Returns:
; 
; 
; :Params:
; 
;   expr, in, type=string
;
;      The expression as a string, with X denoting the independent
;      variable and P[i] the ith element of the parameter array P.
;      Example: 'P[0] + X*P[1] + X*X*P[2]'.
; 
;   x, in, type=array
; 
;      The independent variable X.
; 
;   p, in, type=
; 
;      The parameter array.
; 
; 
; :History:
; 
;   2019-07-15 : MGL. first version.
; 
;-
function red_renderexpr, expr, p, format = format

  Nparam = n_elements(p)
  strexpr = strlowcase(expr)
  
  for i = 0, Nparam-1 do begin
    if n_elements(format) eq 0 then $
       pstring = strtrim(p[i], 2) else $
          pstring = string(p[i], format = format) 
    strexpr = red_strreplace(strexpr, 'p['+strtrim(i, 2)+']', pstring)
  endfor                        ; i
  
  return, strexpr
  
end

