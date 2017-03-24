; docformat = 'rst'

;+
; Build a list of parameter settings useful as prpara in the
; headerinfo_addstep method.
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
;    prpara : in, out, type=strarr
; 
;       The list of parameter settings.
;   
;    paraname : in, type=string
;
;       The name of the added parameter.
; 
;    parameter : in, type=various
;
;       The parameter value.
; 
; :History:
; 
;    2017-03-24 : MGL. First version.
; 
; 
;-
pro red_make_prpara, prpara, paraname, parameter
  
  if n_elements(paraname) eq 0 then return

  case n_elements(parameter) of
    0 : red_append, prpara, paraname
    1 : red_append, prpara, paraname + '=' + strtrim(parameter, 2)
    else : begin
      if typename(parameter) eq 'STRING' then begin
        red_append, prpara, paraname + '=' + '[' + strjoin(parameter, ',') + ']'
      endif else begin
        red_append, prpara, paraname + '=' + red_collapserange(parameter)
      endelse
    end
  endcase

end



