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
;        before sriting.
;
;    status : out, type=signed int
;
;        Output the return status, 0 for success, -1 for failure.
;
; :History:
; 
;   2016-05-26 : THI. First version
;
;-
pro red_writedata, filename, $
                   data, $
                   header = header, $
                   filetype = filetype, $
                   structheader = structheader, $
                   status = status, $
                   overwrite = overwrite, $
                   compress = compress

    if ~keyword_set(overwrite) && file_test(filename) then begin
        message, 'File exists: '+filename + ' (use /overwrite to replace)', /info
        status = -1
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
                writefits, filename, data, header2
            endif else begin
                writefits, filename, data
            endelse
        end

    endcase
  
    status = 0


end
