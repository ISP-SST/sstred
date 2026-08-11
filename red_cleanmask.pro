; docformat = 'rst'

;+
; Clean a mask from large corner areas. Intended for masks used with
; red_fillpix.
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, ISP, 2012-10-07
; 
; 
; :Returns:
; 
; 
; :Params:
; 
;    dirtymask : in, type="2D array"
;   
;      A 2D mask interpreted as "true" or "on" where non-zero and
;      "false" or "off" where zero. So it could be a truly boolean
;      mask or just a gain table with zeroed bad pixels.
;   
; 
; :Keywords:
; 
; 
; 
; :History:
; 
;   2012-10-07 : Written by Mats Löfdahl (MGL).
;
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2013-06-10 : MGL. Rename the argument and add some documentation. 
;
;-
function red_cleanmask, dirtymask, circular_aperture=circ

  nkern = 5
  dims = size(dirtymask, /dim)

  mask = dirtymask ne 0

  oc = bytarr(dims[0]+2, dims[1]+2)
  oc[1, 1] = morph_open(morph_close(mask,replicate(1,nkern,nkern)),replicate(1,nkern,nkern))
  labels = (label_region(oc, /all_neighbors))[1:dims[0], 1:dims[1]]

  for i = 0, 100 do begin 
    aa = where(labels eq i, cnt) 
    if cnt eq 0 then break 
;    print, i
  endfor                        ; i
  mg = fltarr(i)
  for ii = 0, i-1 do mg[ii] = mean(dirtymask(where(labels eq ii)))

  mmm = min(mg, minlabel)

  if keyword_set(circ) then begin
    labels[0:nkern/2-1,*] = labels[nkern/2, nkern/2]
    labels[*,0:nkern/2-1] = labels[nkern/2, nkern/2]
    labels[dims[0]-nkern/2-1:dims[0]-1,*] = labels[nkern/2, nkern/2]
    labels[*,dims[1]-nkern/2-1:dims[1]-1] = labels[nkern/2, nkern/2]
    
    circ= (labels ne labels[nkern/2, nkern/2])

    idx = where(labels eq labels[0,0], count)
    labels += 1
    labels[idx]= 0
    
    return, mask eq 0 or labels eq 0
  endif
  

  
  
  return, mask or (labels eq minlabel)
end
