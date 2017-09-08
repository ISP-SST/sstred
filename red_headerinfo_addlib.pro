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
;    prbranch : in, optional, type=string
; 
;      If given, this will be added as a new PRBRAiia FITS header parameter.
; 
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
;    2017-05-31 : MGL. New keyword prbranch, optionally add the
;                 version control branch.
; 
;    2017-06-01 : MGL. Use red_fitsaddpar.
;
;    2017-09-07 : MGL. Changed red_fitsaddpar --> red_fitsaddkeyword. 
; 
;-
pro red_headerinfo_addlib, head, prlib, prver, prevkey = prevkey, prbranch = prbranch

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
    ;; This is the "main" library, no letter added
    letter = ''
  endif else begin
    ;; We should add a letter. But which letter? And after which
    ;; keyword should we add the next?
    altkey = 'PRVER'+stp
    tmp = fxpar(head, altkey, count = count2)
    if count2 gt 0 then key = altkey
    letternumber = 64B           ; One before 'A'
    repeat begin
      letternumber++
      letter = string(letternumber)
      prevkey = key
      key = 'PRLIB'+stp+letter
      tmp = fxpar(head, key, count = count)
      altkey = 'PRVER'+stp+letter
      tmp = fxpar(head, altkey, count = count2)
      if count2 gt 0 then key = altkey
    endrep until count eq 0
  endelse


  ;; Add the library
  key = 'PRLIB'+stp+letter
  if letter eq '' then begin
    if n_elements(prlibcomment) eq 0 then prlibcomment = 'Software library'
    if prcount gt 0 then prlibcomment += ' containing '+prproc
  endif else prlibcomment = 'Additional software library'
  red_fitsaddkeyword, anchor = anchor, head, key, prlib, prlibcomment, after = prevkey

  ;; Add the version
  if n_elements(prver) gt 0 then begin
    key = 'PRVER'+stp+letter
    red_fitsaddkeyword, anchor = anchor, head, key $
              , prver, 'Library version/MJD of last update' 
  endif

  ;; Optionally add branch
  if n_elements(prbranch) gt 0 then begin
    key = 'PRBRA'+stp+letter
    red_fitsaddkeyword, anchor = anchor, head, key $
              , prbranch, 'Version control branch' 
  endif
  
end



mkhdr, head, 0
red_fitsaddkeyword, head, 'PRSTEP1', 'First step!'
red_fitsaddkeyword, head, 'PRPROC1', 'firstproc'

red_headerinfo_addlib, head, 'monkey','39';, prevkey = 'DATE'
red_headerinfo_addlib, head, 'kitten','0.345';, prevkey = prevkey
red_headerinfo_addlib, head, 'hedgehog','1.4.3';, prevkey = prevkey

red_fitsaddkeyword, head, 'PRSTEP2', 'Second step!'
red_fitsaddkeyword, head, 'PRPROC2', 'secondproc'

red_headerinfo_addlib, head, 'monkey','39';, prevkey = 'DATE'
red_headerinfo_addlib, head, 'kitten','0.345';, prevkey = prevkey
red_headerinfo_addlib, head, 'hedgehog','1.4.3';, prevkey = prevkey

printarr,head[where(head ne strjoin(replicate(' ',80)))]

end
