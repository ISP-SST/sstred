; docformat = 'rst'

;+
; Get a frame from an SST crispex LP file.
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
; 
; :Params:
;
;    filename_or_fileassoc : in, type="integer or string"
;
;       If a string, the name of the file from which to get the frame.
;       Otherwise the assoc variable set up to access the main data
;       part of the already opened file.
;
;    frame : out, type=array
;
;       The read image frame.
; 
; 
; :Keywords:
; 
;    iframe : in, optional, type=integer, default="based on ituning, istokes, and iscan"
;
;       The frame index in the data cube seen as a 3D cube. 
;
;    ituning  : in, optional, type=integer, default=0
;
;       The tuning index, used to calculate iframe.
;
;    istokes  : in, optional, type=integer, default=0
;
;       The stokes index, used to calculate iframe.
;
;    iscan : in, optional, type=integer, default=0
;
;       The scan index, used to calculate iframe.   
;   
;    nscans : in, optional, type=integer
;
;       The number of scans in the file. If not given, will try to get
;       it by parsing the file name.
; 
; :History:
; 
;     2018-11-27 : MGL. First version.
; 
;-
pro red_lpcube_getframe, filename_or_fileassoc, frame $
                         , iframe = iframe $
                         , ituning = ituning $
                         , istokes = istokes $
                         , iscan = iscan $
                         , nscans = nscans

  ;; Did we get a file name
  open_and_close = size(filename_or_fileassoc, /tname) eq 'STRING'
  
  if open_and_close then begin
    ;; We have the file name, open the file and set up an assoc
    ;; variable.
    filename = filename_or_fileassoc

    ;; Get the cube size and variable type from the header, but Nt
    ;; is actually Nstokes*Nwav*Nscans.
    red_lp_header, filename, header=header, datatype=datatype, $
                   dims=dims, nx=Nx, ny=Ny, nt=Nt, endian=endian_file

    
    if ((byte(1L, 0, 1))[0] eq 1) then endian = 'l' else endian='b'
    swap_endian = (datatype gt 1) and (endian ne endian_file)
    openr, lun, filename, /get_lun, swap_endian=swap_endian
    case datatype of            ; Header is 512 bytes
      2 : fileassoc = assoc(lun, intarr(Nx,Ny,/nozer), 512)
      4 : fileassoc = assoc(lun, fltarr(Nx,Ny,/nozer), 512)
      else : stop
    end
  endif else begin
    ;; We have an assoc variable, get array dimensions from the file.
    fileassoc = filename_or_fileassoc
    lun = (size(fileassoc,/struc)).file_lun
    fs = fstat(lun)
    filename = fs.name
    ;; Get the cube size and variable type from the header, but Nt
    ;; is actually Nstokes*Nwav*Nscans.
    red_lp_header, filename, header=header, datatype=datatype, $
                   dims=dims, nx=Nx, ny=Ny, nt=Nt, endian=endian_file
  endelse

  ;; Calculate Nstokes, Nwav, and Nscans
  if strmatch(header,'stokes=\[I,Q,U,V]*') then Nstokes = 4 else Nstokes = 1
  if n_elements(Nscans) eq 0 then begin
    ;; Can we get it from the file name?
    scans = (stregex(filename,'scans=([0-9]+-[0-9]+)',/extract,/subexpr))[1]
    if scans eq '' then stop
    s_array = red_expandrange(scans)
    Nscans = n_elements(s_array)
  end
  Ntuning = Nt/(Nstokes*Nscans)

  if n_elements(iframe) eq 0 then begin

    if n_elements(ituning) eq 0 then ituning = 0L
    if n_elements(istokes) eq 0 then istokes = 0L
    if n_elements(iscan)   eq 0 then iscan   = 0L
    
    dimensions = [Nx, Ny, Ntuning, Nstokes, Nscans]
  
    ;; Calculate the frame number
    iframe = long(ituning) + long(istokes)*Ntuning $
             + long(iscan)*Ntuning*Nstokes
  endif

  frame = fileassoc[iframe]

  ;; Close if we opened.
  if open_and_close then free_lun, lun
  
end
