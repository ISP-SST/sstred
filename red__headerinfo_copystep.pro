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
; 
; 
; 
;-
pro red::headerinfo_copystep, hdr, oldhdr $
                              , all = all $
                              , anchor = anchor $
                              , stepnum = stepnum $
                              , last = last $
                              , prstep = prstep 

  ;; Existing steps
  prsteps_existing = fxpar(hdr,'PRSTEP*')
  Nexisting = n_elements(prsteps_existing)
  
  ;; Available steps
  prsteps_available = fxpar(oldhdr,'PRSTEP*')
  Navailable = n_elements(prsteps_available)

  ;; Any steps to copy?
  if n_elements(prsteps_available) eq 0 then return
  
  ;; Combine steps indicated by the keywords
  if keyword_set(all) then red_append, stepnums, indgen(Navailable)+1
  if keyword_set(last) then red_append, stepnums, Navailable
  if n_elements(stepnum) gt 0 then red_append, stepnums, stepnum
  if n_elements(prstep) gt 0 then begin
    for i = 0, n_elements(prstep)-1 do begin
      ii = where(prsteps eq prstep, Nmatch)
      if Nmatch gt 0 then red_append, stepnums, ii
    endfor   
  endif 
  ;; Some step could be indicated by multiple keywords
  stepnums = stepnums(uniq(stepnums, sort(stepnums))) 
   
  inew = Nexisting
  if Nexisting gt 0 then begin
    anchor = 'PRBRA'+strtrim(Nexisting, 2) ; This should be the last existing
  endif else begin
    anchor = 'OBS_HDU'          ;'SOLARNET'
  endelse
  
  for i = 0, n_elements(stepnums)-1 do begin
    
    iold = stepnums[i]          ; Step number to get from old header
    inew++                      ; Step number to add to new header

    ;; We want to preserve the order these keywords were written by
    ;; headerinfo_addstep so we process them in that order:

    value = red_fitsgetkeyword(oldhdr, 'PRSTEP'+strtrim(iold, 2), comment = comment, count = count)
    if count eq 1 then red_fitsaddkeyword, anchor = anchor, hdr, 'PRSTEP'+strtrim(inew, 2), value, comment
    
    value = red_fitsgetkeyword(oldhdr, 'PRPROC'+strtrim(iold, 2), comment = comment, count = count)
    if count eq 1 then red_fitsaddkeyword, anchor = anchor, hdr, 'PRPROC'+strtrim(inew, 2), value, comment

    value = red_fitsgetkeyword(oldhdr, 'PRMODE'+strtrim(iold, 2), comment = comment, count = count)
    if count eq 1 then red_fitsaddkeyword, anchor = anchor, hdr, 'PRMODE'+strtrim(inew, 2), value, comment

    value = red_fitsgetkeyword(oldhdr, 'PRPARA'+strtrim(iold, 2), comment = comment, count = count)
    if count eq 1 then red_fitsaddkeyword, anchor = anchor, hdr, 'PRPARA'+strtrim(inew, 2), value, comment

    value = red_fitsgetkeyword(oldhdr, 'PRLIB'+strtrim(iold, 2), comment = comment, count = count)
    if count eq 1 then red_fitsaddkeyword, anchor = anchor, hdr, 'PRLIB'+strtrim(inew, 2), value, comment

    value = red_fitsgetkeyword(oldhdr, 'PRVER'+strtrim(iold, 2), comment = comment, count = count)
    if count eq 1 then red_fitsaddkeyword, anchor = anchor, hdr, 'PRVER'+strtrim(inew, 2), value, comment

    foreach letter, ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J'] do begin

      value = red_fitsgetkeyword(oldhdr, 'PRLIB'+strtrim(iold, 2)+letter, comment = comment, count = count)
      if count eq 0 then break
      if count eq 1 then red_fitsaddkeyword, anchor = anchor, hdr, 'PRLIB'+strtrim(inew, 2)+letter, value, comment

      value = red_fitsgetkeyword(oldhdr, 'PRVER'+strtrim(iold, 2)+letter, comment = comment, count = count)
      if count eq 1 then red_fitsaddkeyword, anchor = anchor, hdr, 'PRVER'+strtrim(inew, 2)+letter, value, comment

    end 

    value = red_fitsgetkeyword(oldhdr, 'PRBRA'+strtrim(iold, 2), comment = comment, count = count)
    if count eq 1 then red_fitsaddkeyword, anchor = anchor, hdr, 'PRBRA'+strtrim(inew, 2), value, comment


  endfor                        ; i

end
