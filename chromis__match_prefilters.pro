; docformat = 'rst'

;+
; Check if two lists of wideband and narrowband prefilters belong to
; matching scans.
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;    Tomas Hillberg
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
; 
; 
; :History:
; 
;    2016-09-30 : MGL. Allow different lengths if one of them is 1.
; 
;    2016-11-25 : MGL. Allow NB filters to match themselves.
; 
;-
function chromis::match_prefilters, pf1, pf2
;                         NB      WB
  prefilter_table = [['3999', '3950'], $
                     ['3969', '3950'], $
                     ['3978', '3950'], $
                     ['3934', '3950'], $
                     ['3925', '3950'], $
                     ['4862', '4846']]

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
           print, 'chromis::match_prefilters: prefilter lists must be of the same length.'
           return, 0
        end
     endcase
  endelse

  ret = bytarr(Npref)
  for ip=0, Npref-1 do begin
     ret[ip] = max(where( (prefilter_table[0,*] eq list1[ip] and prefilter_table[1,*] eq list2[ip]) $
                          or (prefilter_table[1,*] eq list1[ip] and prefilter_table[0,*] eq list2[ip]) $
                          or (prefilter_table[1,*] eq list1[ip] and prefilter_table[1,*] eq list2[ip]) $
                          or (prefilter_table[0,*] eq list1[ip] and prefilter_table[0,*] eq list2[ip]) ) $
                  ) ge 0
  endfor

  if n_elements(ret) eq 1 then return, ret[0]
  
  return, ret

end

a = chromisred()

prefilter_table = [['3999', '3950'], $
                   ['3969', '3950'], $
                   ['3978', '3950'], $
                   ['3934', '3950'], $
                   ['3925', '3950'], $
                   ['4862', '4846']]

pf1 = reform(prefilter_table[0, *])
pf2 = reform(prefilter_table[1, *])
help, pf1, pf2
print, pf1
print, pf2
print, a -> match_prefilters(pf1, pf2)

print

pf1 = reform(replicate('3978', n_elements(pf2)))
print, pf1
print, pf2
print, a -> match_prefilters(pf1, pf2)


pf1 = '3978'
print, pf1
print, pf2
print, a -> match_prefilters(pf2, pf1)


pf1 = ['3978', '4846']
print, pf1
print, pf2
print, a -> match_prefilters(pf2, pf1) ; Supposed to fail!

print, a -> match_prefilters('3999', '3999')

end
