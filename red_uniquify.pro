; docformat = 'rst'

;+
; Remove duplicates from a 1D array.
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
;   A sorted list where all elements are unique.
; 
; :Params:
; 
;    arr : in, type=array
; 
;      The input 1D array.
; 
; :Params:
; 
;    indx : out, optional, type=array 
; 
;      The indices of the selected elements.
; 
; :History:
; 
;   2024-12-18 : MGL. First version.
; 
;-
function red_uniquify, arr, indx = indx

  indx = uniq(arr, sort(arr))
  
  return, arr[indx]

end

;; Usage:

a = [2,1,3,2,5,1,0,4,0]
print, a
print, red_uniquify(a)
print
b = ['cat',      'gnu' ,      'dog',      'ape',  'dog',      'dog',      'cat'     ]
c = ['domestic', 'wild',      'domestic', 'wild', 'domestic', 'domestic', 'domestic']
print, b
print, c
print, red_uniquify(b, indx = indx)
print, c[indx]


end
