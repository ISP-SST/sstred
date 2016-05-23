; docformat = 'rst'

;+
; Read headers from data files in a camera independent way.
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;     J. Lewis Fox, ISP, 2016-05-18
; 
; 
; :Returns:
; 
;     A FITS compatible header (from a data file of any instrument)
; 
; :Params:
; 
;    fname : in, type=string
;
;       The data file.
;
; :Keywords:
;    
;    filetype : in, type=string
;
;	The type of file to read. Allowed values are: 
;		ptgrey-fits
;		fz
;	If not set auto-detection will be attempted.
;
;    structheader : in, type=flag
;
;	If set output the header as a structure.
;
;    status : out, type=signed int
;
;	Output the return status, 0 for success, -1 for failure.
;
; :History:
; 
;   2016-05-18 : JLF. Created.
;
;-
function red_readhead, fname, $
                       filetype = filetype, $
                       structheader = structheader, $
                       status = status

    if( file_test(fname) eq 0 ) then begin
        message, 'File does not exist: '+fname,/info
        status = -1
        return, 0B
    endif

    if( n_elements(filetype) eq 0 ) then begin
    
        filetype = rdx_filetype(fname)

        if( filetype eq '' ) then begin
            message, 'Cannot detect filetype. Pass it manually as', /info
            message, "head = red_readhead('"+fname+"',filetype=fits')", /info
            status = -1
            return, 0B
        endif
        
    endif                         ; filetype

        
    case strupcase(filetype) of

        'ANA' : begin
            ;; Data stored in ANA files.
            anaheader = fzhead(fname)
            if n_elements(anaheader) ne 0 then $           
            ;; Convert ana header to fits header
            header = red_anahdr2fits(anaheader)
        end

        'FITS' : begin
            ;; Data stored in fits files
            red_rdfits, fname, header = header
        end

    endcase

    header = red_filterchromisheaders(header)
    
    if n_elements(header) ne 0 and keyword_set(structheader) then begin
        header = red_paramstostruct(header)
    endif

    status = 0

    return, header

end
