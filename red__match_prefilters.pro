; docformat = 'rst'

;+
; Check if two lists of wideband and narrowband prefilters belong to
; matching scans.
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
;    An array of booleans. 
;
; :Params:
; 
;    pf1 : type=(string or strarr)
;
;       A list of prefilters, wb or nb.
;
;    pf2 : type=(string or strarr)
;
;       A list of prefilters, wb or nb.
; 
; 
; :History:
; 
;    2020-05-28 : MGL. Allow different lengths if one of them is 1.
; 
;-
function red::match_prefilters, pf1, pf2

  N1 = n_elements(pf1)
  N2 = n_elements(pf2)

  if N1 eq 0 or N2 eq 0 then return, 0

  if N1 eq N2 then begin
    Npref = N1
    list1 = pf1
    list2 = pf2
  endif else begin
    case 1 of
      N1 : begin
        Npref = N2
        list1 = replicate(pf1[0], Npref)
        list2 = pf2
      end
      N2 : begin
        Npref = N1
        list1 = pf1
        list2 = replicate(pf2[0], Npref)
      end
      else: begin
        print, 'red::match_prefilters : prefilter lists must be of the same length.'
        return, 0
      end
    endcase
  endelse

    
  return, pf1 eq pf2
  
end
