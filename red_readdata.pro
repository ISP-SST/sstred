; docformat = 'rst'

;+
; Camera-independent read routine.
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;     Mats Löfdahl, ISP, 2016-05-16
; 
; 
; :Returns:
; 
;    The image contained in the file fname.
; 
; :Params:
; 
;    fname : in, type=string
;
;       The file name.
;
; :Keywords:
;
;    header : out, type=string array/struct
;
;  	Output the header. 
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
;   2016-05-16 : MGL. First version
; 
;   2016-05-17 : MGL. Default for filetype based on file name.
;
;   2016-05-18 : JLF. Returns the image (like readfits). Added Status flag keyword.
;
;   2016-05-18 : JLF. Use red_mkhdr to create FITS headers from ANA headers.
;-
function red_readdata, fname $
		       , header = header $
                       , filetype = filetype $
                       , structheader = structheader $
                       , status = status

  if file_test(fname) eq 0 then begin

    message, 'File does not exist: '+fname,/info
    status = -1
    return, 0B


  endif

  if n_elements(filetype) eq 0 then begin

     if strmatch(fname, '*.fits') then begin
      ;; Try to find out if it's a PointGrey camera based on file name
      cam = (strsplit(file_basename(fname), '.', /extract))[0]
      if n_elements(cam) ne 0 then caminfo = red_camerainfo(cam)
      if strmatch(caminfo.model,'PointGrey*') then begin
	  filetype = 'ptgrey-fits'
      endif else begin			; PointGrey
	message, 'Cannot detect filetype. Pass it manually as',/info
	message,"head = red_readhead('"+fname+"',filetype='ptgrey-fits')",/info
	statue = -1
	return, 0B
      endelse
        
     endif else begin
        
        filetype = 'fz'         
        
     endelse

  endif                         ; filetype


  case filetype of

     'fz' : begin
        
        ;; Data stored in ANA fz format files.
        
	data = f0(fname)
	if arg_present(header) then anaheader = fzhead(fname)

        if n_elements(anaheader) ne 0 then $           
           ;; Convert ana header to fits header
           header = red_anahdr2fits(anaheader,img=data)

     end

     'ptgrey-fits' : begin

        ;; Data stored in fits files, but in the strange format
        ;; returned by the PointGrey cameras.

	red_rdfits, fname, image = data, header = header, /uint, swap=0
	data = ishft(data, -4) ; 12 bits in 16-bit variables

     end

  endcase

  if n_elements(header) ne 0 and keyword_set(structheader) then begin
     header = red_paramstostruct(header)
  endif
  
  status = 0
  
  return, data

end
