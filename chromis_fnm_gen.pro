; docformat = 'rst'

;+
; Construct IDL strings to generate directory and file names for CHROMIS files
; with aid of 'execute' function. The strings are to be stored in the database 
; instead of full filenames.
; 
; Derived (including comments) from crisp_filename_template.pro by Mats Löfdahl.
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Oleksii Andriienko, Institute for Solar Physics
; 
; 
; :Params:
; 
;    filename : in, type=string
; 
;      The file name sample.
; 
; 
; :Keywords:
; 
;   
;   fnm_gen : out, optional, type=string
;
;      The constructed IDL string to generate (with aid of execute function) the base filename.
;
;   dir_gen, optional, type=string
;   
;      The constructed IDL string to generate the directory name (starting from
;      the date).
; 
; :History:
; 
;    2019-03-12 : MGL. First version of crisp_filename_template.pro.
;    2019-03-15 : OA. Rename the procedure. Switch to IDL strings generation instead of 
;                 abstract templates.
;-
pro chromis_fnm_gen, filename $
                             , fnm_gen = fnm_gen $
                             , dir_gen = dir_gen $
                             , verbose = verbose
 
  timeregex = '[0-2][0-9]:[0-5][0-9]:[0-6][0-9]'
  dateregex = '20[0-2][0-9][.-][01][0-9][.-][0-3][0-9]'

  basename = file_basename(filename)
  dirname = file_dirname(filename)

  dirname_format = 'FORMAT = "('
  dirname_vars = ''
  filename_format = 'FORMAT = "('
  filename_vars = ''

  ;; Construct the template for the file directory
  if dirname ne '' and arg_present(dir_gen) then begin

    ;; Find the date and remove everything before it. That part is
    ;; site dependent.
    date = stregex(dirname, dateregex, /extract) 
    pos = strpos(dirname, date)
    dirname = strmid(dirname, pos)

    dirname_parts = strsplit(dirname, '/', /extract)

    ;; First part is always the date. But are the parts delimited by
    ;; dots of dashes? We don't know since Tomas decided to change
    ;; from dots to dashes and we don't know if other sites have done
    ;; this as well. So it's site dependent and some sites might have
    ;; a mix... We'll construct the template according to the
    ;; format in the input dir name. But programs using this have to
    ;; be aware that they might have to check which form is on disk.
    Ndots   = n_elements(strsplit(dirname_parts[0], '.', /extract))
    Ndashes = n_elements(strsplit(dirname_parts[0], '-', /extract))
    if Ndots gt Ndashes then begin
      ;dirname_format += 'C(CYI4,".",CMOI2,".",CDI2), A1, '
      ;dirname_vars += 'date, "/",'
      dirname_format += 'I4, A1, 2(I02, A1), ' 
      dirname_vars += 'Y, ".", M, ".", D, "/", '
    endif else begin
      ;dirname_format += 'C(CYI4,"-",CMOI2,"-",CDI2), A1, '
      ;dirname_vars += 'date, "/",'  
      dirname_format += 'I4, A1, 2(I02, A1), '
      dirname_vars += 'Y, "-", M, "-", D, "/", '
    endelse  
    
    ;; Second part is most often the type (darks, flats, polcal,
    ;; pinholes, pf-scan, etc.) but it can also be a free-format
    ;; string describing a particular type of science data. So we'll
    ;; make this part of the template a literal string. This means we
    ;; have to store one template per directory in the database! Any
    ;; unidentified type should be stored as "science data" in the
    ;; data base.
    dirname_format += 'A, A, A1, '
    dirname_vars += '"CHROMIS", datatype, "/", '

    ;; Followed by the timestamp and the camera
    ;dirname_format += 'C(CHI2.2, ":", CMI2.2, ":", CSI2.2), A1, A, A1)"'
    ;dirname_vars += '"time, "/", camera, "/"'
    dirname_format += '3(I02,A1), A, A1, '
    dirname_vars += 'H, ":", M, ":", S, "/", camera, "/", '
       
    ll = strlen(dirname_format)
    dirname_format = strmid(dirname_format, 0, ll-2) + ')"'
    
    dir_gen = 'dir = string(' + dirname_vars + dirname_format + ')'

    if keyword_set(verbose) then begin
      print
      print, dirname
      print, dir_gen
    endif
    
  endif 


  ;; Construct the template for the file basename
  if basename ne '' and arg_present(fnm_gen) then begin

    basename_parts = strsplit(basename, '_', /extract)
    ; Find out datatype
    if n_elements(basename_parts) eq 4 then begin ; darks or WB flats for which no simultaneous NB data were collected. 
      filename_format += '3A, I05, A1, I07, A)"'
      filename_vars += '"sst_", detector, "_", scannum, "_", first_frame, ".fits", '
    endif else begin ; science(data), flats, pinholes
      filename_format += '3A, I05, A, I07, A, I05, A, I05, A)"'
      filename_vars += '"sst_", detector, "_", scannum, "_", first_frame, "_wheel", wheel, "_hrz", hrz, ".fits", '  
    endelse
 
    fnm_gen = 'fnm = string(' + filename_vars + filename_format + ')'       

    if keyword_set(verbose) then begin
      print
      print, basename
      print, fnm_gen
    endif
    
  endif

end


paths = '/data/2018/2018-05/2018-05-26/' $
        + [$
        'CHROMIS-darks/15:24:35/Chromis-D/sst_camXXVII_00000_0000000.fits', $
        'CHROMIS-pinholes/15:20:04/Chromis-D/sst_camXXVII_00000_0000000_wheel00005_hrz32122.fits'  $                
          ]

detector='camXXVII'
camera='Chromis-N'
datatype='Darks'
scannum=12
first_frame=200
wheel='wheel00005'
hrz='hrz32122'
Y=2018
M=7
D=3
H=8
M=10
S=1
for i = 0, n_elements(paths)-1 do begin

  print, '--------------------------------'
  chromis_fnm_gen, paths[i] $ ;, /verbose $
                           , fnm_gen = fnm_gen $
                           , dir_gen = dir_gen
  v=execute(fnm_gen)                       
  if v then print, 'generated filename: ', fnm else print, 'Error!'
  v=execute(dir_gen)
  if v then print, 'generated directory name: ', dir else print, 'Error!'
endfor


end
