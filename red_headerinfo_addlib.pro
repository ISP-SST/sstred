; docformat = 'rst'

;+
; Add library info to a FITS header.
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, ISP
; 
; :Params:
; 
;    head : in, out, type=strarr
; 
;      A FITS header.
; 
;    prlib : in, type=string
; 
;      The library name.
;
;    prver : in, type=string
; 
;      The library version number.
;
; :Keywords:
; 
;    prevkey : in, out, optional, type=string
;   
;      The keyword to place the next keyword after.
; 
; 
; :History:
; 
;    2017-03-16 : MGL. First version.
; 
;-
pro red_headerinfo_addlib, head, prlib, prver, prevkey = prevkey

  ;; What is the first nonexisting step number?
  stepnumber = 0
  repeat begin 
    stepnumber++ 
    stp = strtrim(stepnumber, 2) 
    tmp = fxpar(head, 'PRSTEP'+stp, count = count)
  endrep until count eq 0 
  ;; Go back to the first that does exist.
  stepnumber-- 
  stp = strtrim(stepnumber, 2)

  if stepnumber eq 0 then begin
    print, 'You should add PRSTEP (and possibly PRPROC) before adding libraries.'
    stop
  endif
  
  prproc = fxpar(head, 'PRPROC'+stp, count = prcount)

  if n_elements(prevkey) eq 0 then begin
    if prcount gt 0 then prevkey = 'PRPROC'+stp else prevkey = 'PRSTEP'+stp
  endif

  ;; Now check if this should be the "main" library or an "additional"
  ;; library with A,B,... added.
  
  key = 'PRLIB'+stp
  tmp = fxpar(head, key, count = count)
  if count eq 0 then begin
    letter = ''
  endif else begin
    ;; We should add a letter. But which letter?
    letternumber = 64B           ; One before 'A'
    repeat begin
      letternumber++
      letter = string(letternumber)
      key = 'PRLIB'+stp+letter
      tmp = fxpar(head, key, count = count)
    endrep until count eq 0
    prevkey = key
  endelse


  ;; Add the library
  key = 'PRLIB'+stp+letter
  if letter eq '' then begin
    if n_elements(prlibcomment) eq 0 then prlibcomment = ' Software library'
    if prcount gt 0 then prlibcomment += ' containing '+prproc
  endif else prlibcomment = ' Additional software library'
  fxaddpar, head, key, prlib, prlibcomment, after = prevkey
  prevkey = key

  ;; Add the version
  if n_elements(prver) gt 0 then begin
    key = 'PRVER'+stp+letter
    fxaddpar, head, key, after = prevkey $
              , prver, ' Library version/MJD of last update' 
    prekey = key
  endif

end



mkhdr, head, 0
fxaddpar, head, 'PRSTEP1', 'First step!'
fxaddpar, head, 'PRPROC1', 'firstproc'

red_headerinfo_addlib, head, 'monkey','39';, prevkey = 'DATE'
red_headerinfo_addlib, head, 'kitten','0.345';, prevkey = prevkey
red_headerinfo_addlib, head, 'hedgehog','1.4.3';, prevkey = prevkey

fxaddpar, head, 'PRSTEP2', 'Second step!'
fxaddpar, head, 'PRPROC2', 'secondproc'

red_headerinfo_addlib, head, 'monkey','39';, prevkey = 'DATE'
red_headerinfo_addlib, head, 'kitten','0.345';, prevkey = prevkey
red_headerinfo_addlib, head, 'hedgehog','1.4.3';, prevkey = prevkey

printarr,head[where(head ne strjoin(replicate(' ',80)))]

end
