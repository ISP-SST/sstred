; docformat = 'rst'

;+
; Copy FITS info about a processing step from one FITS header to
; another. 
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
; 
; :Params:
; 
;    hdr : in, type=strarr
;
;      The FITS header in which to write the step info.
; 
;    oldhdr : in, type=strarr
;
;      The FITS header from which to read the step info.
; 
; 
; :Keywords:
; 
;    all : in, optional, type=boolean
;   
;       Copy info for all available steps.
; 
;    anchor  : in, optional, type=string
;   
;       Position the step info after the specified keyword.
; 
;    stepnum  : in, optional, type=integer
;   
;       Copy info for the step with this (these) number(s). 
; 
;    last  : in, optional, type=boolean
;   
;       Copy info for the last available step.
; 
;    prstep  : in, optional, type=string
;   
;       Copy info for steps matching this string.
; 
; 
; :History:
; 
;    2018-10-31 : MGL. First version.
; 
;    2020-07-01 : MGL. Reimplement.
; 
;-
pro red::headerinfo_copystep, hdr, oldhdr $
                              , all = all $
                              , anchor = anchor $
                              , stepnum = stepnum $
                              , last = last $
                              , prstep = prstep 

  ;; Existing steps
  prsteps_existing = strtrim(fxpar(hdr,'PRSTEP*', count = Nexisting), 2)
  
  ;; Available steps
  prsteps_available = strtrim(fxpar(oldhdr,'PRSTEP*', count = Navailable), 2)

  ;; Any steps to copy?
  if Navailable eq 0 then return

  prkeys = ['PRSTEP', 'PRPROC', 'PRMODE', 'PRPARA', 'PRLIB', 'PRVER', 'PRREF', 'PRBRA']
  
  ;; Combine steps indicated by the keywords
  if keyword_set(all) then red_append, stepnums, indgen(Navailable)+1
  if keyword_set(last) then red_append, stepnums, Navailable
  if n_elements(stepnum) gt 0 then red_append, stepnums, stepnum
  for istep = 0, n_elements(prstep)-1 do begin
    ii = where(prsteps_available eq prstep[istep], Nmatch)
    if Nmatch gt 0 then red_append, stepnums, ii+1
  endfor                        ; istep

  ;; Some step could be indicated by multiple keywords
  stepnums = stepnums(uniq(stepnums, sort(stepnums))) 
   
  inew = Nexisting
  if Nexisting gt 0 then begin

    ;; Look for last PR* keyword to use as anchor
    pos = 0
    keywrd = strmid(oldhdr, 0, 8)
    iend = n_elements(keywrd)
    for ikey = 0, n_elements(prkeys)-1 do begin
      for istep = 0, n_elements(stepnums)-1 do begin
        foreach letter, ['', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J'] do begin
          name = prkeys[ikey]+strtrim(stepnums[istep], 2)+letter
          thispos = fxparpos(keywrd, iend, before = name)
          
          if thispos eq iend then break
          if thispos gt pos then begin
            pos = thispos
            anchor = name
          endif
        end
      endfor                    ; istep
    endfor                      ; ikey
    
  endif else begin
    anchor = 'OBS_HDU'          ;'SOLARNET'
  endelse

  for istep = 0, n_elements(stepnums)-1 do begin
    inew++                      ; Step number to add to new header
    for ikey = 0, n_elements(prkeys)-1 do begin
      foreach letter, ['', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J'] do begin

        name = prkeys[ikey]+strtrim(stepnums[istep], 2)+letter

        print, name
        
        value = red_fitsgetkeyword(oldhdr, name, comment = comment, count = count)
        if count eq 0 then break
        red_fitsaddkeyword, anchor = anchor, hdr, prkeys[ikey]+strtrim(inew, 2), value, comment
        
      end                       ; letter 
    endfor                      ; ikey
  endfor                        ; istep

end


a = crispred(/dev)

hdr = headfits('/scratch/mats/2016.09.19/CRISP-aftersummer/cubes_wb/wb_6302_2016-09-19T09:30:20_scans=2-8_corrected_im.fits')

mkhdr, newhdr, fltarr(10, 10)

hgrep, hdr, 'PR'

print

a -> headerinfo_copystep, newhdr, hdr, stepnum = 1
hprint, newhdr

a -> headerinfo_copystep, newhdr, hdr, stepnum = 2
hprint, newhdr

end
