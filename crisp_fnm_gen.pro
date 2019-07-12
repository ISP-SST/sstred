; docformat = 'rst'

;+
; Construct IDL strings to generate directory and file names for CRISP files
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
pro crisp_fnm_gen, filename $
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
    ;red_append, dirname_template_parts, dirname_parts[1]
    ;;red_append, dirname_template_parts, '<TYPE>'
    dirname_format += 'A, A1, '
    dirname_vars += 'datatype, "/", '

    ;; For polcal data we have to add the prefilter
    if strmatch(dirname_parts[1], 'Polcal') then begin
      ; sometimes (it seems that only for 5173 prefilter)
      ; timestamp is after Polcal and prefilter part should be omitted
      if ~strmatch(dirname_parts[2], timeregex) then begin 
        dirname_format += 'I04, A1, '
        dirname_vars += 'filter1, "/", '
      endif
    endif

    ;; Followed by the timestamp and the camera
    ;dirname_format += 'C(CHI2.2, ":", CMI2.2, ":", CSI2.2), A1, A, A1)"'
    ;dirname_vars += '"time, "/", camera, "/"'
    dirname_format += '3(I02,A1), A, A1, '
    dirname_vars += 'hr, ":", min, ":", sec, "/", camera, "/", '
        
    ll = strlen(dirname_format)
    dirname_format = strmid(dirname_format, 0, ll-2) + ')"'
    ll = strlen(dirname_vars)
    dirname_vars = strmid(dirname_vars, 0, ll-6)
    
    dir_gen = 'dir = string(' + dirname_vars + dirname_format + ')'

    if keyword_set(verbose) then begin
      print
      print, dirname
      print, dir_gen
    endif
    
  endif 


  ;; Construct the template for the file basename
  if basename ne '' and arg_present(fnm_gen) then begin

    ;; Loop through the dot-delimited parts of the file basename and
    ;; try to identify the parts.
    basename_parts = strsplit(basename, '.', /extract)
    for ipart = 0, n_elements(basename_parts)-1 do begin
      case 1 of
        ;; 4-digit prefilter
        strmatch(basename_parts[ipart], '[0-9][0-9][0-9][0-9]') : begin
          filename_format += 'I04, A1, '
          filename_vars += 'filter1, ".", '
        end
        ;; 5-digit scan number
        strmatch(basename_parts[ipart], '[0-9][0-9][0-9][0-9][0-9]') : begin 
          filename_format += 'I05, A1, '
          filename_vars += 'scannum, ".", '
        end
        ;; 7-digit frame number
        strmatch(basename_parts[ipart], '[0-9][0-9][0-9][0-9][0-9][0-9][0-9]') : begin 
          filename_format += 'I07, A1, '
          filename_vars += 'framenum, ".", '
        end
        ;; Camera tag
        strmatch(basename_parts[ipart] , 'cam[XVI]*') : begin 
          filename_format += 'A, A1, '
          filename_vars += 'detector, ".", '
        end
        ;; Liquid crystal state
        strmatch(basename_parts[ipart], 'lc[0-9]') : begin 
          filename_format += 'A2, I1, A1,'
          filename_vars += '"lc", lc_state, ".", '
        end
        ;; Linear polarizer state
        strmatch(basename_parts[ipart], 'LP[0-9][0-9][0-9]') : begin 
          filename_format += 'A2, I03, A1, '
          filename_vars += '"LP", lp_state, ".", '
        end
        ;; Quarterwave plate state
        strmatch(basename_parts[ipart], 'qw[0-9][0-9][0-9]') : begin 
          filename_format += 'A2, I03, A1, '
          filename_vars += '"qw", qw_state, ".", '
        end
        ;; DM focus
        strmatch(basename_parts[ipart], 'f[+-][0-9][0-9][0-9]') : begin 
          filename_format += 'A1, I+03, A1, '
          filename_vars += '"f", focus, ".", '
        end
        ;; Tuning
        strmatch(basename_parts[ipart], '[0-9][0-9][0-9][0-9]_[+-][0-9]*') : begin 
          ;; The tuning info is a four-digit LINE followed by an
          ;; underscore and the finetuning. We need to find out how
          ;; many digits are used for the tuning part. For CRISP it is
          ;; always zero-padded and has a sign.
          Ndigits = strlen(basename_parts[ipart]) - 5
          filename_format += 'I04, A1, I+0' + strtrim(Ndigits, 2) + ', A1, '
          filename_vars += 'line, "_", tuning, ".", '
        end
        ;; Unidentified parts, add as text.
        else : begin
          filename_format += 'A, A1, '
          filename_vars += '"' + basename_parts[ipart] + '", ".", '
        end
      endcase
    endfor                      ; ipart
    
    ll = strlen(filename_format)
    filename_format = strmid(filename_format, 0, ll-2) + ')"' 
    ll = strlen(filename_vars)
    filename_vars = strmid(filename_vars, 0, ll-6)

    fnm_gen = 'fnm = string(' + filename_vars + filename_format + ')'       

    if keyword_set(verbose) then begin
      print
      print, basename
      print, fnm_gen
    endif
    
  endif

end


paths = '/data/2018/2018.07/2018.07.05/' $
        + [$
        'Darks/13:31:51/Crisp-R/camXIX.00004.im.ex.im.0000970', $
        'Flats/12:13:14/Crisp-T/camXXV.00013.im.ex.6173.6173_-090.lc4.im.0208796', $        
        'Pinholes/13:30:26/Crisp-W/camXX.00000.im.ex.8542.8542_+0850.lc4.im.0000998', $        
        'Polcal/8542/15:53:23/Crisp-R/camXIX.00000.im.ex.LP090.qw360.8542.cont.lc3.im.0023894', $        
        'Science/10:17:53/Crisp-T/camXXV.00003.im.ex.8542.8542_+0000.lc0.im.0005248', $        
        'Science/10:17:53/Crisp-T/camXXV.00003.im.ex.6173.6173_-210.lc2.im.0004204' $        
          ]

; Older type CRISP, with free-form science dir name and focus info:
red_append, paths, '/data/2014/2014-05/2014-05-27/SCI-8542+7772/08:10:30/Crisp-W/camXX.00010.im.ex.f-004.8542.8542_-1750.lc1.im.0013666'

detector='camXIX'
camera='Crisp-R'
datatype='Science'
scannum=12
framenum=234
filter1=6173
line=6173
tuning=340
lc_state=1
lp_state=45
qw_state=30
focus=2
Y=2018
M=7
D=3
H=8
M=10
S=1
for i = 0, n_elements(paths)-1 do begin

  print, '--------------------------------'
  crisp_fnm_gen, paths[i], /verbose $
                           , fnm_gen = fnm_gen $
                           , dir_gen = dir_gen
  v=execute(fnm_gen)                       
  if v then print, 'generated filename: ', fnm else print, 'Error!'
;  v=execute(dir_gen)
;  if v then print, 'generated directory name: ', dir else print, 'Error!'
endfor


end
