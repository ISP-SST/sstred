; docformat = 'rst'

;+
; Read step info from a fits header.
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
;   The PRSTEP info as a hash.
; 
; :Params:
; 
;   header : in, type=strarr
; 
;     A FITS file header. If only one element or a scalar string,
;     interpreted as the name of a file from which to read the header.
; 
; :Keywords:
; 
;   count : out, optional, type=integer
;    
;      The number of returned steps.
;
;   prkey : in, optional, type=string
;
;      Return only keywords matching this string.
; 
;   prstep : in, optional, type=string
;
;      Return only steps matching this step name.
; 
; :History:
; 
;   2021-05-27 : MGL. First version.
; 
;-
function red_headerinfo_getstep, header $
                                 , count = count $
                                 , prkey = prkey $
                                 , prstep = prstep

  if n_elements(header) eq 0 then begin
    count = 0
    print, 'No header provided.'
    return, ''
  endif

  hdr = header                  ; Protect input string
  
  if n_elements(hdr) eq 1 then begin
    ;; Is it a file name?
    if file_test(hdr[0]) then begin
      hdr = headfits(hdr[0])
    endif else begin
      ;; Nope
      print, 'Not a file name.'
      print, hdr[0]
      count = 0
      return, ''
    endelse
  endif

  ;; Read the step names
  prsteps = fxpar(hdr, 'PRSTEP*', count = count)
  if count eq 0 then return, '' ; No steps in the header

  ;; Did we call with keyword prstep?
  if n_elements(prstep) eq 1 then begin
    sindx = where(prsteps eq prstep, count)
    if count eq 0 then begin
      print, 'No such step name.'
      print, prstep
      return, ''                ; No matching step names    
    endif
  endif

;  prkeys = ['PRSTEP', 'PRPROC', 'PRMODE', 'PRPARA', 'PRLIB', 'PRVER', 'PRREF', 'PRBRA', 'PRHSH', 'PRENV']
;  Nkeys = n_elements(prkeys)

  prkeys = red_headerinfo_prkeys(count = Nkeys) ; List of defined PR* keywords.
    
  h8 = strmid(hdr, 0, 8)        ; Just the header keywords, 8 first characters.

  
  ;; Find the keywords
  for i = 0, count-1 do begin   ; Loop over selected steps
    print, i
    undefine, keys
    for ikey = 0, Nkeys-1 do begin ; Loop over possible keys
      indx = where(strmatch(h8, prkeys[ikey]+strtrim(sindx[i]+1, 2)+'*'), Nmatch)
      if Nmatch gt 0 then begin
        for imatch = 0, Nmatch-1 do begin
          this_key = strtrim(h8[indx[imatch]], 2)
          if n_elements(prkey) then begin
            if ~strmatch(this_key, prkey+'*') then continue
          endif
          red_append, keys, this_key
        endfor                  ; imatch
      endif
    endfor                      ; key

  endfor                        ; i

  count = n_elements(keys)
  if count eq 0 then return, ''

  ;; Get the keyword info for the selected keys
  output = hash()
  for ikey = 0, count-1 do begin
    value = fxpar(hdr, keys[ikey], comment = comment)
    output[strtrim(keys[ikey], 2)] = {value:value, comment:comment}
  endfor                        ; ikey

  return, output
  stop
  
end

fname = '/scratch/mats/prefilter_bug/2019-04-19/CRISP/cubes_nb/nb_6173_2019-04-19T11:29:24_scans=0,1_stokes_corrected_im.fits'
header = headfits(fname)

print, 888

stepinfo = red_headerinfo_getstep(header, prstep = 'CALIBRATION-INTENSITY-TEMPORAL', count = count, prkey = 'PRPROC')

end
