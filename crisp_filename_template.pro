; docformat = 'rst'

;+
; Construct a name template from a sample CRISP raw file name.
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
;   basename_template : out, optional, type=string
;
;      The constructed template for the base filename.
;
;   dirname_template, optional, type=string
;   
;      The constructed template for the directory name (starting from
;      the date).
; 
; :History:
; 
;    2019-03-12 : MGL. First version.
; 
;-
pro crisp_filename_template, filename $
                             , basename_template = basename_template $
                             , dirname_template = dirname_template $
                             , verbose = verbose
  
  timeregex = '[0-2][0-9]:[0-5][0-9]:[0-6][0-9]'
  dateregex = '20[0-2][0-9][.-][01][0-9][.-][0-3][0-9]'

  basename = file_basename(filename)
  dirname = file_dirname(filename)

  dirname_template = '?'
  basename_template = '?'

  ;; Construct the template for the file directory
  if dirname ne '' and arg_present(dirname_template) then begin

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
      red_append, dirname_template_parts, '<DOTDATE>'
    endif else begin
      red_append, dirname_template_parts, '<ISODATE>'
    endelse
    
    ;; Second part is most often the type (darks, flats, polcal,
    ;; pinholes, pf-scan, etc.) but it can also be a free-format
    ;; string describing a particular type of science data. So we'll
    ;; make this part of the template a literal string. This means we
    ;; have to store one template per directory in the database! Any
    ;; unidentified type should be stored as "science data" in the
    ;; data base.
    red_append, dirname_template_parts, dirname_parts[1]
    ;;red_append, dirname_template_parts, '<TYPE>'

    ;; For polcal data we have to add the prefilter
    if strmatch(dirname_parts[1], 'Polcal') then begin
      red_append, dirname_template_parts, '<FILTER1(04)>'
    endif

    ;; Followed by the timestamp and the camera
    red_append, dirname_template_parts, ['<TIME-OBS>', '<CAMERA>']
    
    dirname_template = strjoin(dirname_template_parts, '/') + '/'

    if keyword_set(verbose) then begin
      print
      print, dirname
      print, dirname_template
    endif
    
  endif 


  ;; Construct the template for the file basename
  if basename ne '' and arg_present(basename_template) then begin

    ;; Loop through the dot-delimited parts of the file basename and
    ;; try to identify the parts.
    basename_parts = strsplit(basename, '.', /extract)
    for ipart = 0, n_elements(basename_parts)-1 do begin
      case 1 of
        ;; 4-digit prefilter
        strmatch(basename_parts[ipart] $
                       , '[0-9][0-9][0-9][0-9]') : $ 
           red_append, basename_template_parts $
                       , '<FILTER1(04)>'
        ;; 5-digit scan number
        strmatch(basename_parts[ipart] $
                       , '[0-9][0-9][0-9][0-9][0-9]') : $ 
           red_append, basename_template_parts $
                       , '<SCANNUM(05)>'
        ;; 7-digit frame number
        strmatch(basename_parts[ipart] $
                       , '[0-9][0-9][0-9][0-9][0-9][0-9][0-9]') : $ 
           red_append, basename_template_parts $
                       , '<FRAMENUM(07)>'
        ;; Camera tag
        strmatch(basename_parts[ipart] $
                       , 'cam[XVI]*') : $ 
           red_append, basename_template_parts $
                       , '<CAMERA>'
        ;; Liquid crystal state
        strmatch(basename_parts[ipart] $
                       , 'lc[0-9]') : $ 
           red_append, basename_template_parts $
                       , 'lc<LCSTATE(1)>'
        ;; Linear polarizer state
        strmatch(basename_parts[ipart] $
                       , 'LP[0-9][0-9][0-9]') : $ 
           red_append, basename_template_parts $
                       , 'LP<LPSTATE(03)>'
        ;; Quarterwave plate state
        strmatch(basename_parts[ipart] $
                       , 'qw[0-9][0-9][0-9]') : $ 
           red_append, basename_template_parts $
                       , 'qw<QWSTATE(03)>'
        ;; DM focus
        strmatch(basename_parts[ipart] $
                       , 'f[+-][0-9][0-9][0-9]') : $ 
           red_append, basename_template_parts $
                       , 'f<FOCUS(+03)>'
        ;; Tuning
        strmatch(basename_parts[ipart] $
           , '[0-9][0-9][0-9][0-9]_[+-][0-9]*') : begin 
          ;; The tuning info is a four-digit LINE followed by an
          ;; underscore and the finetuning. We need to find out how
          ;; many digits are used for the tuning part. For CRISP it is
          ;; always zero-padded and has a sign.
          Ndigits = strlen(basename_parts[ipart]) - 6
          red_append, basename_template_parts $
                      , '<LINE(04)>_<TUNING(+0'+strtrim(Ndigits, 2)+')>'
        end
        ;; Unidentified parts, add as text.
        else : $
           red_append, basename_template_parts $
                       , basename_parts[ipart]
      endcase
    endfor                      ; ipart
    
    
    basename_template = strjoin(basename_template_parts, '.') 

    if keyword_set(verbose) then begin
      print
      print, basename
      print, basename_template
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

;; Older type CRISP, with free-form science dir name and focus info:
red_append, paths, '/data/2014/2014-05/2014-05-27/SCI-8542+7772/08:10:30/Crisp-W/camXX.00010.im.ex.f-004.8542.8542_-1750.lc1.im.0013666'

for i = 0, n_elements(paths)-1 do begin

  print, '-----'
  crisp_filename_template, paths[i], /verbose $
                           , basename_template = basename_template $
                           , dirname_template = dirname_template
  
endfor


end
