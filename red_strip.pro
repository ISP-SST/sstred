; docformat = 'rst'

;+
; Remove elements from an array.
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl & Tomas Hillberg
; 
; 
; :Params:
; 
; 
;   array : in, type=array
;   
;      The 'original' array.
; 
;   data : in, type=array
;   
;      The elements to remove from the original array. If data has
;      elements that are not in array, they are silently ignored.
; 
; 
; :History:
; 
;     2016-10-13 : MGL & THI. First version
;
;-
pro red_strip, array, data

   match2, array, data, suba, subb
   
   indx = where(suba eq -1, nindx)
   if nindx gt 0 then array = array[indx]     

end
