; docformat = 'rst'

;+
; Round x to # of digits or # of decimals.
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
;   x : in, type=number
; 
; 
; 
; 
; :Keywords:
; 
;   n_decimals : in, optional, type=integer
;   
;      The number of decimals to preserve. If negative, return an
;      integer with the -n_decimals least significant digits set to
;      zero.
; 
;   n_digits : in, optional, type=integer
;   
;      The number of digits to preserve. If no decimals, return an
;      integer.
; 
; 
; :History:
; 
;   2026-06-04 : MGL. First version.
; 
;-
function red_round, x, n_digits = n_digits, n_decimals = n_decimals

  case !true of

    n_elements(n_decimals) gt 0 : begin
      if n_decimals gt 0 then begin
        return, round(x*10.^n_decimals)/10.^n_decimals
      endif else begin
        ;; If no decimals, return an integer
        return, round(round(x*10.^n_decimals)/10.^n_decimals)
      endelse
    end
    
    n_elements(n_digits) gt 0 : begin
      return, red_round(x, n_decimals = n_digits - floor(alog10(x)) - 1)
    end
    
    else : return, round(x)

  endcase
  
end


n = 2
x = !dpi

print
print, x, 1e3*x, 1e-3*x
print
print, 'Normal round'
print, red_round(x), red_round(1e3*x), red_round(1e-3*x)
print
print, strtrim(n, 2)+ ' decimals:'
print, red_round(x, n_decimals = n), red_round(1e3*x, n_decimals = n), red_round(1e-3*x, n_decimals = n)

print
print, strtrim(n, 2)+ ' digits:'
print, red_round(x, n_digits = n), red_round(1e3*x, n_digits = n), red_round(1e-3*x, n_digits = n)

print
print, strtrim(-n, 2)+ ' decimals:'
print, red_round(x, n_decimals = -n), red_round(1e3*x, n_decimals = -n), red_round(1e-3*x, n_decimals = -n)


end

