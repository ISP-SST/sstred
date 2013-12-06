;	   File:       collapserange_f.ana
;	   Created:    <2003-06-19 15:20:16 mats>
;	   Author:     mats@astro.su.se
;	   Time-stamp: <2009-06-11 15:47:51 mats>

; Ported from ANA  (not done yet)

FUNCTION red_collapserange, arr, ld = ld, rd = rd
  
  if n_elements(ld) eq 0 then ld = '['
  if n_elements(rd) eq 0 then rd = ']'
  
  ;; Simple cases, one or two elements
  if n_elements(arr) eq 1 then return, ld + strtrim(arr[0], 2) + rd
  if n_elements(arr) eq 2 then return, ld + strjoin(string(arr, format = '(i0)'), ',') + rd
  
  strng = ld+strcompress(string(arr[0]), /rem)
  
  FOR i = 1, n_elements(arr)-2 DO BEGIN
     IF arr[i] EQ arr[i-1]+1 THEN BEGIN
        IF arr[i] NE arr[i+1]-1 THEN BEGIN
                                ;strng += '!'
           IF ~strmatch(strng, '*[,-]') THEN strng += ','	
           strng += strcompress(string(arr(i)), /rem) + ','
        END ELSE BEGIN
           if ~strmatch(strng, '*-') then strng += '-'
        ENDELSE
     END ELSE BEGIN
        IF ~strmatch(strng, '*,') then strng += ','	
        strng += strcompress(string(arr[i]), /rem)
     ENDELSE
  ENDFOR
  
  ;;IF (last(strng) NE ',') AND (last(strng) NE '-') THEN strng += ','	
  if ~strmatch(strng, '*[,-]') then strng += ',' 

  ;;return, strng+strcompress(string(last(arr)), /rem)+rd
  return, strng + strcompress(string(arr[n_elements(arr)-1]), /rem) + rd
  
end                             ;
