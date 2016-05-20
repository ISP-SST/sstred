; docformat = 'rst'

;+
; Clean up a FITS header to make it conform to standards.
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;     J. Lewis Fox, ISP, 2016-05-19
; 
; 
; :Returns:
; 
;     A FITS compatible header from an input header.
; 
; :Params:
; 
;    hdr : in, type=string array
;
;       The header, string array(36) from a FITS file.
;
; :Keywords:
;
; :History:
; 
;   2016-05-19 : JLF. Created.
;
;-
pro red_cleanheader, dirty_hdr

hdr = dirty_hdr

;; keys are in front of the first =
keys = gettok(hdr,'=')

values = hdr
comments = hdr

;; separate values from comments on the /
values = gettok(comments,'/')

;; find tick characters, there could be a correctly formatted string
hdr = strtrim(hdr,1)
firsttick = strpos(hdr,"'")
lasttick = strpos(hdr,"'",/reverse_search)

;; only a pair of ticks is a proper string
pc = where(firsttick ne lasttick,cnt)

if cnt ne 0 then begin
  for i = 0, n_elements(pc)-1 do begin
  
    ;; value is between the ticks
    values[pc[i]] = strmid(hdr[pc[i]],firsttick[pc[i]]+1,lasttick[pc[i]]-1)
    
    ;; comment is after that with or without a slash
    comments[pc[i]] = strmid(hdr[pc[i]],lasttick[pc[i]]+1)
    
    ;; if the first non-blank in the comment is / strip it
    ;; sxaddpar is going to add one.
    if strpos(strtrim(comments[pc[i]],1),'/') eq 0 then $
      comments[pc[i]] = strmid(strtrim(comments[pc[i]],1),1)
  
  endfor
endif

;; this is mostly for boolean values ('T'/'F') but also makes
;; string values prettier. Trailing blanks are not significant
;; according to spec.
values = strtrim(values,0)

for i = 0, 35 do begin
  
  ;; spec is a ' inside a string is represented with '' 
  values[i] = strjoin(strsplit(values[i],"'",/extract),"''")
  
  ;; ignore keywords with empty values. This excludes standard keys
  ;; like END, HISTORY, and COMMENT, but also leaves alone random
  ;; lines of text inserted in the header which don't contain an =
  ;; 
  ;; It's an open question if things like that should be cleaned.
  if values[i] ne '' then begin
  
    ;; try to replace keywords which might not be all-caps with their
    ;; upper case version.
    ;; *** NOT WORKING ***
    if stregex(keys[i],'[a-z]',/bool) then sxdelpar,dirty_hdr,keys[i]
    
    out_value = values[i]
    
    ;; check for numerical values.
    ;; ints are cast to the smallest type that will contain them
    ;; floats and doubles are cast to doubles
    ;; 
    ;; we don't handle complex values yet!
    if valid_num(values[i],val,/int) || $
       valid_num(values[i],val) then out_value = val
    sxaddpar,dirty_hdr,keys[i],out_value,comments[i],after=keys[i-1>0]
  endif
endfor


end 

;; test the code by .run red_cleanheader

blank = string(replicate(32b,80))
test_header = replicate(blank,36)

str='SIMPLE  = T'
str=[str,'BITPIX  = 16']
str=[str,'NAXIS   = 3']
str=[str,'NAXIS1  = 1920']
str=[str,'NAXIS2  = 1200']
str=[str,'NAXIS3  = 100']
str=[str,'ORIGIN  = Dutch Open Telescope']
str=[str,'TELESCOP= Dutch Open Telescope']
str=[str,'INSTRUME= CHROMIS-FOO']
str=[str,'DATE_OBS= 2016-05-07T11:14:21']
str=[str,'DATE    = 2016-05-07T11:15:14']
str=[str,'OBSERVER']
str=[str,'OBJECT  = disk center']
str=[str,'EXPTIME = 0.006006 / [s]']
str=[str,'INTERVAL= 0.020000 / [s]']
str=[str,'GAIN    = 0.000000']
str=[str,'OFFSET  = 0.000000']
str=[str,'COMMENT   For measuring relative intensities']
str=[str,'TESTBOOL=                    F /']
str=[str,"TESTSTRING = 'This 'is / 'a test' Str const with comment."]
str=[str,"TESTSTR2 = 'This is also' / a test' / With a leading slash."]
str=[str,"TESTSTR3 = 'Another test string / with only one tick."]
str=[str,'TESTSTR4 = before first / is the value / after is comment.']
str=[str,'TESTSTR5 = before first = is the key.']
str=[str,"teststr6 = 'before first = is the key,' / proper str"]
str=[str,"TESTSTR7 =        'string constant' / with extra space"]
str=[str,"TESTSTR8 = a ' in / the middle"]
str=[str,'END']

for i = 0, n_elements(str)-1 do begin
    blank = test_header[i]
    strput,blank,str[i]
    test_header[i] = blank
endfor

clean_hdr = test_header
red_cleanheader,clean_hdr

print,clean_hdr

end