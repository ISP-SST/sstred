; docformat = 'rst'

;+
; Helper procedure for red_writedata. Add all the write stuff to the
; main header and create the bintable header.
;
; The bintable header, bdr, comes out in the tabhdu.extname
; sub-structure with key name 'bdr'.
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;    Lewis Fox, ISP
; 
; 
; :Returns:
; 
; 
; :Params:
; 
; 
; :Keywords:
; 
;   
;   
;   
; 
; 
; :History:
;
;   2016-09-16 : JLF. First version.
; 
;   2016-12-07 : MGL. Split into a separate file.
; 
;-
pro red_create_extheaders,tabhdu,extname,head

  top_tags = tag_names(tabhdu)
  top_idx = where(top_tags eq extname,cnt)
  if cnt eq 0 then message,"Can't find "+extname+' in Tab-HDU structure.'

  struct = tabhdu.(top_idx)

  keys = tag_names(struct)

  ;; Extname comment, if it exists
  idx = where(keys eq 'COMMENT',cnt)
  if cnt ne 0 then extcomment = struct.(idx) else extcomment=''

  keys = red_strreplace(keys,'_','-')

  ;; Exclude possible BDR and COMMENT keys
  cidx = where(keys ne 'BDR' and keys ne 'COMMENT',cnt,comp=nk_idx)
  if cnt eq 0 then message,'No tabulated keywords in '+extname

  untabkeys = keys[nk_idx]      ; not tabulated keywords
  keys = keys[cidx]

  ;; Create the bintable header
  fxbhmake,bdr,1,extname,extcomment

  ;; Over-head for fxparpos
  keywd = strmid(head,0,8)
  iend = where(keywd eq 'END     ')

  for i = 0, n_elements(keys)-1 do begin
                                ; the tabulated keyword sub-structure
    tabkeywd = struct.(cidx[i]) 
    if size(tabkeywd,/type) ne 8 then $
       message,'Structure tag '+keys[i]+' is not a structure.'
    
    ;; Update the main header with the tabulated keywords. After is
    ;; 'OBS_SHDU' or keys[i-1] whichever is further down
    if i eq 0 then after = 'OBS_SHDU' else begin
      after = keys[i-1]
      if fxparpos(keywd,iend,after=after) lt fxparpos(keywd,iend,after='OBS_SHDU') then $
         after = 'OBS_SHDU'
    endelse
    
    tags = tag_names(tabkeywd)
    
    ;; VAL tag is the only required tag. make sure it exists before we
    ;; start using it.
    idx = where(tags eq 'VAL',cnt)
    if cnt eq 0 then message,'VAL tag is required in '+keys[i]
    
    ;; Summary value (what goes in the main header keyword value)
    idx = where(tags EQ 'SUMMARY',cnt)
    if cnt ne 0 then summary_value = tabkeywd.(idx) else $
       summary_value = mean(tabkeywd.val)
    
    ;; Comment (main header comment)
    idx = where(tags eq 'COMMENT',cnt)
    if cnt ne 0 then comment = tabkeywd.(idx) else comment = ''
    
    ;; Add the keyword to the header
    sxaddpar,head,keys[i],summary_value,comment,after=after

    ;; Add tabulated values to the bintable
    fxbaddcol,i+1,bdr,tabkeywd.val,keys[i],comment
  endfor

  ;; Add the tabulatd keyword according to Stein Vidar's solarnet
  ;; tab-hdu extension schema.
  tabulatd = extname+'; '+strjoin(keys,', ')
  sxaddpar,head,'TABULATD',tabulatd,after=keys[i-1]

  ;; Add the bintable header to the structure so we can use it later.
  ;; Is there already a bintable header in the structure?
  junk = where(untabkeys eq 'BDR',cnt)
  if cnt ne 0 then struct.bdr = bdr else struct = create_struct('bdr',bdr,struct)

  if n_elements(top_tags) eq 1 then newstruct = create_struct(extname,struct) $
  else begin
    newstruct = top_idx eq 0 ? create_struct(top_tags[0],struct) : $
                create_struct(top_tags[0],tabhdu.(0))
    for i = 1, n_elements(top_tags)-1 do $
       newstruct = i eq top_idx ? create_struct(newstruct,top_tags[i],struct) : $
       create_struct(newstruct,top_tags[i],tabhdu.(i))
  endelse

  tabhdu = newstruct

end
