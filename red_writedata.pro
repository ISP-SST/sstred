; docformat = 'rst'

;+
; Camera-independent write routine.
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;     Tomas Hillberg, ISP, 2016-05-26
; 
; 
; :Returns:
; 
; 
; :Params:
; 
;    fname : in, type=string
;
;       The file name.
;
; :Keywords:
;
;    header : in, type=strarr
;
; 
;    filetype : in, type=string
;
;        The type of file to write. Allowed values are: 
;        "ana" or "fits"
;        If not set auto-detection will be attempted.
;
;    structheader : in, type=flag
;
;        If set, convert the header from struct to strarr
;        before writing.
;
;    status : out, type=signed int
;
;        Output the return status, 0 for success, -1 for failure.
;
;    tabhdu: in, type=struct
;
;	 struct containing the header keywords expressed as binary 
;	 table vectors.
;	 tabhdu.[extnames] = 1 extension per substruct (named [extname])
;	 tabhdu.[extnames].comment = comment for the extname keyword
;	 tabhdu.[extnames].[tabkeys] = 1 tabluated keyword per sub-substruct
;					(named [tabkey], also used for ttypen)
;	 tabhdu.[extnames].[tabkeys].{val,comment,summary}:
;			val = the array to be written to the table
;			comment = main hdr comment and ttypen comment
;			summary = value to put in main hdr for this keyword
;	 sub-structure keys that don't exist are replaced with defaults, blank
;	 for comment, mean(val) for summary
;
; :History:
; 
;   2016-05-26 : THI. First version
;
;   2016-09-16 : JLF. Added tabhdu, binary table writing.
;		 Also fixed a /overwrite bug (it'd overwrite even if you didn't
;		 use /overwrite)
;
;-

; helper procedure to add all the write stuff to the main header
; and create the bintable header.
;
; The bintable header, bdr, comes out in the tabhdu.extname sub-structure
; with key name 'bdr'.
pro red_create_extheaders,tabhdu,extname,head

top_tags = tag_names(tabhdu)
top_idx = where(top_tags eq extname,cnt)
if cnt eq 0 then message,"Can't find "+extname+' in Tab-HDU structure.'

struct = tabhdu.(top_idx)

keys = tag_names(struct)

; extname comment, if it exists
idx = where(keys eq 'COMMENT',cnt)
if cnt ne 0 then extcomment = struct.(idx) else extcomment=''

keys = red_strreplace(keys,'_','-')

;exclude possible BDR and COMMENT keys
cidx = where(keys ne 'BDR' and keys ne 'COMMENT',cnt,comp=nk_idx)
if cnt eq 0 then message,'No tabulated keywords in '+extname

untabkeys = keys[nk_idx] ; not tabulated keywords
keys = keys[cidx]

; create the bintable header
fxbhmake,bdr,1,extname,extcomment

; over-head for fxparpos
keywd = strmid(head,0,8)
iend = where(keywd eq 'END     ')

for i = 0, n_elements(keys)-1 do begin
  ; the tabulated keyword sub-structure
  tabkeywd = struct.(cidx[i]) 
  if size(tabkeywd,/type) ne 8 then $
    message,'Structure tag '+keys[i]+' is not a structure.'
    
  ; update the main header with the tabulated keywords.
  ; after is 'OBS_SHDU' or keys[i-1] whichever is further down
  if i eq 0 then after = 'OBS_SHDU' else begin
    after = keys[i-1]
    if fxparpos(keywd,iend,after=after) lt fxparpos(keywd,iend,after='OBS_SHDU') then $
      after = 'OBS_SHDU'
  endelse
  
  tags = tag_names(tabkeywd)
  
  ; VAL tag is the only required tag. make sure it exists before we start
  ; using it.
  idx = where(tags eq 'VAL',cnt)
  if cnt eq 0 then message,'VAL tag is required in '+keys[i]
  
  ; summary value (what goes in the main header keyword value)
  idx = where(tags EQ 'SUMMARY',cnt)
  if cnt ne 0 then summary_value = tabkeywd.(idx) else $
    summary_value = mean(tabkeywd.val)
  
  ; comment (main header comment)
  idx = where(tags eq 'COMMENT',cnt)
  if cnt ne 0 then comment = tabkeywd.(idx) else comment = ''
  
  ; add the keyword to the header
  sxaddpar,head,keys[i],summary_value,comment,after=after

  ; add tabulated values to the bintable
  fxbaddcol,i+1,bdr,tabkeywd.val,keys[i],comment
endfor

; add the tabulatd keyword according to Stein Vidar's solarnet tab-hdu
; extension schema.
tabulatd = extname+'; '+strjoin(keys,', ')
sxaddpar,head,'TABULATD',tabulatd,after=keys[i-1]

; add the bintable header to the structure so we can use it later.
; is there already a bintable header in the structure?
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

; helper procedure to do the actual writing of the bintable.
;
; some funniness involved with uint data. FXB* adds the TSCALE,TZERO
; keywords for uints, but FXBWRITE doesn't know what to do with the 
; uint data (it doesn't test for type code 12) so you have to scale it
; down yourself.
pro red_write_tabhdu,tabhdu,filename

extnames = tag_names(tabhdu)

for i = 0, n_elements(extnames)-1 do begin
  ; create the extension table in the file.
  fxbcreate,lun,filename,tabhdu.(i).bdr,ext_no
  
  ; get column data keyword names
  keys = tag_names(tabhdu.(i))
  ; COMMENT and BDR aren't data keywords, everything else is.
  idx = where(keys ne 'COMMENT' and keys ne 'BDR')
  
  ; write the data column by column, everything is in row 1
  for j = 0, n_elements(idx)-1 do begin
    data = tabhdu.(i).(idx[j]).val
    ; fxbwrite can't handle uint data!?!
    if size(data,/type) eq 12 then data = fix(data-32768)
    fxbwrite,lun,data,j+1,1 ; FITS tables are unit indexed
  endfor
  ; done with this extension
  fxbfinish,lun
  
endfor

end

pro red_writedata, filename, $
                   data, $
                   header = header, $
                   filetype = filetype, $
                   structheader = structheader, $
                   status = status, $
                   overwrite = overwrite, $
                   compress = compress, $
                   tabhdu = tabhdu

    if ~keyword_set(overwrite) && file_test(filename) then begin
        message, 'File exists: '+filename + ' (use /overwrite to replace)', /info
        status = -1
        return
    endif

    if n_elements(header) ne 0 then begin
        header2 = header
        if keyword_set(structheader) then begin
            header2 = red_paramstostruct(header, /inverse)
        endif
        ; make sure it is a valid fits header.
        check_fits, data, header2, /UPDATE, /SILENT
    endif
    
    if n_elements(filetype) eq 0 then begin

        filetype = rdx_filetype(filename)
  
        if filetype eq '' then begin
            message,'Cannot detect filetype. Pass it manually as', /info
            message, "red_writedata,'"+filename+"',data,filetype='fits'", /info
            status = -1
            return
        endif
        
    endif                         ; filetype


    case strupcase(filetype) of

        'ANA' : begin
            if n_elements(header2) ne 0 then begin
                fzwrite, data, filename, strjoin(temporary(header2),''), compress = compress
            endif else begin
                fzwrite, data, filename, compress = compress
            endelse
        end

        'FITS' : begin
            if n_elements(header2) ne 0 then begin
		if n_elements(tabhdu) ne 0 then begin
		  if ~fxpar(header2,'extend') then $ ; add the EXTEND keyword
		    fxaddpar,header2,'extend','T',' FITS data may contain extensions'
		  extname = tag_names(tabhdu)
		  for i = 0, n_elements(extname)-1 do begin
		    red_create_extheaders,tabhdu,extname[i],header2
		  endfor
		endif
                writefits, filename, data, header2
                ; are there Tab-HDU tabulated keywords?
                if n_elements(tabhdu) ne 0 then begin
		  red_write_tabhdu,tabhdu,filename
                endif
                
            endif else begin
                writefits, filename, data
            endelse
        end

    endcase
  
    status = 0


end

;; fully worked example

datadir = '/storage/sand15n/Incoming/2016.06.21/Flat-tests_000/00:00:00/Chromis-W'

cd,datadir,current=curdir
files = file_search('*')

img = red_readdata(files[1],head=hdr) ; this is a dark, 1 msec exposure time.

print,hdr

dark1 = total(img,3)/500d0
pmm,dark1
print,mean(dark1)

; datamean
dmean = dblarr(1,1,500)
dmin = uintarr(1,1,500)
dmax = dmin
dmed = dmin
for i = 0, 499 do begin 
  dmean[*,*,i] = total(img[*,*,i])/1920d0/1200d0
  dmin[*,*,i] = min(img[*,*,i],max=tmp)
  dmax[*,*,i] = tmp
  dmed[*,*,i] = median(img[*,*,i])
endfor
; date-beg
timestr = strsplit(fxpar(hdr,'date-beg'),'T',/extract)
start = red_time2double(timestr[1])
datebeg = start+fxpar(hdr,'cadence')*dindgen(500)
arr = strarr(2,500)
arr[0,*] = timestr[0]
arr[1,*] = red_time2double(datebeg,/inv)
datebeg = strjoin(temporary(arr),'T')
datebeg = reform(datebeg,[1,1,500])

struct = {tabulations: {date_beg: {val:datebeg,summary:datebeg[0]},$
			datamean: {val:dmean,comment:' [DN] Mean of data'},$
			datamin: {val:dmin,comment:' [DN] Minimum of data',$
				  summary: min(dark1)},$
			datamax: {val:dmax,comment:' [DN] Maximum of data',$
				  summary:max(dark1)},$
			datamedn:{val:dmed,comment:' [DN] Median of data',$
				  summary:median(dark1)},$
			comment: ' For storing tabulated keywords'}}

cd,curdir

red_writedata,(strsplit(red_strreplace(files[1],'.','_'),'.',/extr))[0]+$
  '_reduced_dark.fits',dark1,header=hdr,tabhdu=struct



end