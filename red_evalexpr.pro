; docformat = 'rst'

;+
; Evaluate an expression of the type used with mpfitexpr.
;
; Based on mpevalexpr, a utility function defined in mpfitxpr.pro,
; that sometimes does not compile automatically
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
function red_evalexpr, expr, x, p

  cmd = 'f = ' + expr
  err = execute(cmd)
  return, f
  
end

