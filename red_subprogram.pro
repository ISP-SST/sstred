; docformat = 'rst'

;+
; Return the name of the current level/subprogram and (optionally) the
; name of the level that called it.
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
;    The name of the subprogram/level.
; 
; :Params:
; 
; 
; 
; 
; 
; 
; :Keywords:
; 
;   calling_program : out, optional, type=string
;
;      The name of the calling level/subprogram.
;
;   lowcase : in, optional, type=string
;   
;      Return lowcase strings.
; 
; 
; :History:
; 
;    2017-08-30 : MGL. First verison.
; 
; 
; 
;-
function red_subprogram, calling_program = calling_program, lowcase = lowcase

  levels = (scope_traceback(/structure)).routine
  if keyword_set(lowcase) then levels = strlowcase(levels)

  if n_elements(levels) ge 3 then calling_program = levels[-3] else calling_program = ''
  
  return, levels[-2]
  
end
