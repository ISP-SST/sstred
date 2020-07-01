; docformat = 'rst'

;+
; Delete processing step keywords (PR*na) from a FITS header. 
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
;    2020-06-30 : MGL. First version.
; 
;-
pro red_headerinfo_deletestep, hdr $
                               , all = all $
                               , anchor = anchor $
                               , stepnum = stepnum $
                               , last = last $
                               , prstep = prstep 

  ;; Existing steps
  prsteps_existing = fxpar(hdr,'PRSTEP*')
  Nexisting = n_elements(prsteps_existing)

  if Nexisting eq 0 then return

  if keyword_set(last) then begin
    print, 'Not implemented'
    stop
  endif
  if keyword_set(stepnum) then begin
    print, 'Not implemented'
    stop
  endif

  prkeys = ['PRSTEP', 'PRPROC', 'PRMODE', 'PRPARA', 'PRLIB', 'PRVER', 'PRREF', 'PRBRA']

  if keyword_set(all) then stepnums = indgen(Nexisting)

  if keyword_set(last) then red_append, stepnums, Nexisting-1

  if n_elements(stepnum) gt 0 then red_append, stepnums, stepnum


  for istep = 0, n_elements(prstep)-1 do begin
    ii = where(prsteps_existing eq prstep, Nmatch)
    if Nmatch gt 0 then red_append, stepnums, ii
  endfor   
  
  ;; Some step could be indicated by multiple keywords
  stepnums = stepnums(uniq(stepnums, sort(stepnums))) 

  print, stepnums
  
  for istep = 1, n_elements(stepnums) do begin
    for ikey = 0, n_elements(prkeys)-1 do begin
      foreach letter, ['', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J'] do begin

        name = prkeys[ikey]+strtrim(istep, 2)+letter

        print, name
        
        value = red_fitsgetkeyword(hdr, name, count = count)
        if count eq 0 then break
        red_fitsdelkeyword, hdr, name
        
      end                       ; letter
    endfor                      ; ikey
  endfor                        ; istep
  
end


hdr = headfits('/scratch/mats/2016.09.19/CRISP-aftersummer/cubes_wb/wb_6302_2016-09-19T09:30:20_scans=2-8_corrected_im.fits')

hgrep, hdr, 'PR'
red_headerinfo_deletestep, hdr, /all
print
hgrep, hdr, 'PR'

end

