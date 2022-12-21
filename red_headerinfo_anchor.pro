pro red_headerinfo_anchor, hdr, anchor, next_stepnumber = next_stepnumber

  ;; Find existing steps in the header
  indx = where(strmatch(hdr,'PRSTEP*'), Nexisting)

  if Nexisting eq 0 then begin
    ;; No existing steps. Default anchor, if not input:
    if n_elements(anchor) eq 0 then anchor = 'OBS_HDU'
    next_stepnumber = 1
    return
  endif 

  ;; Get the step number for the last step
  keywrd = strmid(hdr, 0, 8)
  prsteps_existing = strtrim(keywrd[indx], 2)
  stp = max(strmid(prsteps_existing,6)) ; The step number as a string
  next_stepnumber = long(stp)+1
  
  ;; List of defined PR* keywords.
  prkeys = red_headerinfo_prkeys(count = Nkeys) 
   
  ;; Look for last PR* keyword with the last number, to use as anchor
  pos = 0
  iend = n_elements(keywrd)
  for ikey = 0, n_elements(prkeys)-1 do begin
    foreach letter, ['', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J'] do begin
      name = prkeys[ikey]+stp+letter
      thispos = fxparpos(keywrd, iend, before = name)
      if thispos eq iend then break
      if thispos gt pos then begin
        pos = thispos
        anchor = name
      endif
    end                         ; letter
  endfor                        ; ikey
  
end

h = headfits('/scratch/mats/2016.09.19/CRISP-aftersummer/cubes_wb/wb_6302_2016-09-19T09:30:20_09:30:20=1,2_mixed_corrected_im.fits')

;hprint, h

red_headerinfo_anchor, h, anchor, next_stepnumber = next_stepnumber

print, anchor
print, next_stepnumber

end
