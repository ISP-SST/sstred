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
;    maxcount : in, optional, type=integer
; 
;      The maximum number of elements in the output.
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
;     2016-10-05 : MGL. Allow non-string defaults.
;
;     2018-09-17 : MGL. New keyword maxcount.
; 
;-
function red_select_subset, list $
                            , count = count $
                            , default = default $
                            , help = help $
                            , indx = indx $
                            , maxcount = maxcount $
                            , qstring = qstring 
  
  Nelements = n_elements(list)
  
  if Nelements eq 0 then begin
    count = 0
    return, ''
  endif
  
  if n_elements(qstring) eq 0 then qstring = 'Select some'
  if n_elements(default) eq 0 then dflt = '*' else dflt = strtrim(default, 2)
  if n_elements(maxcount) eq 0 then maxcount = Nelements
  
  numbers = string(indgen(Nelements), format = '(i4)')

  repeat begin
  
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
    read, qstring+' [default='+dflt+']: ', selection

    if n_elements(selection) eq 0 then selection = dflt
    if selection eq '' then selection = dflt

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

    if count gt maxcount then begin
      print
      print, '-----> Note that the number of selected items may not exceed '+strtrim(maxcount, 2)
      print
    end
      
  endrep until count le maxcount
  
  return, list(indx)

end


items = ['a', 'b', 'c', 'd']

selection = red_select_subset(items, count = count, default = 1, maxcount = 2)

print, selection

end
