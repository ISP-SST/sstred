; docformat = 'rst'

;+
; Select subset from a list.
; 
; :Categories:
;
;    SST observations
; 
; 
; :Author:
; 
;    Mats Löfdahl, 2013-10-03
; 
; 
; :Returns:
; 
;    The list of selected elements.
; 
; :Params:
; 
;    list : in, type=strarr
; 
;      This is the list to chose from,
; 
; :Keywords:
; 
;    count : out, optional, type=integer
; 
;      The number of elements in the output.
; 
;    help : in, optional, type=boolean
;
;      Set this to get instructions.
;
; :History:
; 
;     2013-10-03 : MGL. Use red_expandrange instead of expandrange.
; 
;     2013-12-05 : MGL. A single dash will select nothing. New keyword
;                  indx for the selected indices.
;
;     2016-09-26 : MGL. If none selected, undefine indx.
; 
; 
; 
;-
function red_select_subset, list, count = count, help = help, qstring = qstring, default = default, indx = indx

  Nelements = n_elements(list)
  
  if Nelements eq 0 then begin
     count = 0
     return, ''
  endif

  if n_elements(qstring) eq 0 then qstring = 'Select some'
  if n_elements(default) eq 0 then default = '*'

  numbers = string(indgen(Nelements), format = '(i4)')
  
  print, numbers+' : '+list, format='(a0)'
  
  ;; Select a subset of the time-stamped directories?
  selection = ''
  if keyword_set(help) then begin
     print, 'Make a selection from the list by entering either of:'
     print, '  1. A number in the range [0, '+string(Nelements-1, format = '(i0)')+'],'
     print, '  2. A comma and dash separated list of numbers in the range,'
     print, '  3. An asterisk to select the entire list.'
     print, '  4. An dash to select nothing.'
  endif 
  read, qstring+' [default='+default+']: ', selection

  if n_elements(selection) eq 0 then selection = default
  if selection eq '' then selection = default

  case selection of 
    
    '-' : begin
      count = 0
      undefine, indx 
      return, ''
    end
    '*' : indx = lindgen(Nelements)
    else : indx = red_expandrange(selection)
  endcase 
  
  count = n_elements(indx)
  
  return, list(indx)

end
