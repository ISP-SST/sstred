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
;   2016-09-16 : JLF. Added tabhdu, binary table writing. Also fixed a
;                /overwrite bug (it'd overwrite even if you didn't use
;                /overwrite)
; 
;   2016-12-07 : MGL. Split binary table helper procedures into
;                separate files.
;
;-
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
          if ~fxpar(header2,'extend') then $ ; Add the EXTEND keyword
             fxaddpar,header2,'extend','T',' FITS data may contain extensions'
          extname = tag_names(tabhdu)
          for i = 0, n_elements(extname)-1 do begin
            red_create_extheaders,tabhdu,extname[i],header2
          endfor
        endif
        writefits, filename, data, header2
        ;; Are there Tab-HDU tabulated keywords?
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

;; Fully worked example

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

struct = {tabulations: {date_beg: {val:datebeg,summary:datebeg[0],unit:'ISODATE (ISO-8601)'},$
                        datamean: {val:dmean,comment:' [DN] Mean of data',unit:'DN'},$
                        datamin: {val:dmin,comment:' [DN] Minimum of data',$
                                  summary: min(dark1)},$
                        datamax: {val:dmax,comment:' [DN] Maximum of data',$
                                  summary:max(dark1)},$
                        datamedn:{val:dmed,comment:' [DN] Median of data',$
                                  summary:median(dark1)},$
                        test:{val:1l,comment:' test scalar'},$
                        testarr:{val:fltarr(1),comment:' test single element array'},$
                        comment: ' For storing tabulated keywords'}}

cd,curdir

red_writedata,(strsplit(red_strreplace(files[1],'.','_'),'.',/extr))[0]+$
              '_reduced_dark.fits',dark1,header=hdr,tabhdu=struct

end
