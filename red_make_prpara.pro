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
;    prpara : in, out, type=dictionary
; 
;       The list of parameter settings.
;   
;    parameter : in, type=various
;
;       The parameter value.
; 
; :Keywords:
;
;    paraname : in, type=string
;
;       The name of the added parameter.
; 
; :History:
; 
;    2017-03-24 : MGL. First version.
; 
;    2017-12-20 : MGL. Re-implementation: generate dictionary instead
;                 of string array. Get parameter name from
;                 scope_varname on the parameter.
; 
; 
;-
pro red_make_prpara, prpara, parameter, paraname = paraname

  if n_elements(paraname) eq 0 then paraname = (strupcase(scope_varname(parameter, level = -1)))[0]
  
  if n_elements(prpara) eq 0 then prpara = dictionary() ; Initialize with an empty dictionary.
  
  if n_elements(paraname) eq 0 then return  ; May have run to just define the dictionary?
  if strtrim(paraname, 2) eq '' then return ; No parameter name
  
  if n_elements(parameter) eq 0 then return ; Set to empty string if no value!

  prpara[strtrim(paraname, 2)] = parameter

end



