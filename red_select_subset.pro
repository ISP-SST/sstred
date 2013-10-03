; docformat = 'rst'

;+
; Select subset from a list.
; 
; :Categories:
;
;    SST observations
; 
; 
; :author:
; 
;    Mats Löfdahl, 2013-10-03
; 
; 
; :returns:
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
; :history:
; 
;     2013-10-03 : MGL. First version.
; 
; 
; 
;-
function red_select_subset, list, count = count, help = help

  Nelements = n_elements(list)
  
  numbers = string(indgen(Nelements), format = '(i4)')
  
  print, numbers+' : '+list, format='(a0)'
  
  if Nelements gt 1 then begin
     ;; Select a subset of the time-stamped directories?
     selection = ''
     if keyword_set(help) then begin
        print, 'Make a selection from the list by entering either of:'
        print, '  1. A number in the range [0, '+string(Nelements-1, format = '(i0)')+'],'
        print, '  2. A comma and dash separated list of numbers in the range,'
        print, '  3. An asterisk to select the entire list.'
     endif 
     read, 'Select some [default=*]: ', selection
     if selection eq '' or selection eq '*' then selection = '0-'+strtrim(string(Nelements-1), 2)
     sindx = expandrange(selection)
  endif else begin
     sindx = indgen(Ndirs)
  endelse

  count = n_elements(sindx)
  
  return, list(sindx)

end
