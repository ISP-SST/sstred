; docformat = 'rst'

;+
; Collection and basic processing of calibration data from CHROMIS and CRISP
; observations.
;
; :Categories:
;
;   SST pipeline
;
; :Author:
;
;   Andrii V. Sukhorukov, Institute for Solar Physics, asukh@astro.su.se
;
; :Params:
;
;   None
;
; :Keywords:
;
;   date        : in, optional, type = string,
;                 default = 'from out_dir if possible'
;
;   out_dir     : in, optional, type = string,
;                 default = 'current working directory'
;
;   search_dirs : in, optional, type = array of strings,
;                 default = '/storage/sand04/... in Stockholm or /data/disk?/...
;                 on La Palma'.
;
; :History:
;
;   2017-07-17 : AVS. Initial commit.
;   2017-07-24 : AVS. TODO: tested calling IDL from inside, doesn't work as
;                     expected.
;   2017-07-27 : AVS. A line-by-line parser is added.
;   2017-07-30 : AVS. TODO: doesn't work on La Palma under the old CRISPRED
;                     installation. Calling hostname -I is not working too.
;   2017-07-31 : AVS. Keyword and directory handling is completely rewritten.
;   2017-08-01 : AVS. The IP-address logic, searching for the instrumental
;                     directories, and the script generation using a common
;                     call_procedure invocation are rewritten.
;   2017-08-02 : AVS. The script parser is updated.  Succesfully tested on
;                     recent data sets in /storage/sand??/Incoming/
;   2017-08-11 : AVS. Added search_dirs keyword to be able to run in different
;                     locations where old data is stored.
;   2017-08-14 : AVS. Changed cd calls to pushd/popd for a better portability.
;   2017-08-24 : AVS. Logic in internal loops have been moved to
;                     red_setupworkdir.  Instead call red_setupworkdir with the
;                     proper keywords and the doit.pro/config.txt pairs will be
;                     created outside this script.
;-
pro red_getcalibdata,       $
  date        = date,       $
  out_dir     = out_dir,    $
  search_dirs = search_dirs ;

  ; Instruments to process data for.  TRIPPEL and SLITJAW are not included yet.
  instruments = [ 'CHROMIS', 'CRISP' ]

  ; The default calling sequence is to specify the output directory in out_dir
  ; and it should contain the date as the last subdirectory, for example,
  ; /scratch/asukh/2017.07.15.

  ; If out_dir is not given, set it to the current directory.
  if ( n_elements( out_dir ) eq 0 ) then begin
    out_dir = getenv( 'PWD' )
    message, 'out_dir is not given and is set to the current working ' + $
      'directory, ' + out_dir, $
      /informational
  endif

  ; out_dir must end with a slash.
  if ~strmatch( out_dir, '*/' ) then out_dir += '/'

  ; Search if there is a date at the end of out_dir.
  pos = stregex( out_dir,                                  $
    '\/[12][0-9][0-9][0-9][.-][01][0-9][.-][0-3][0-9]\/?$' )
  if ( pos eq -1 ) then begin
    out_dir_date = ''
  endif else begin
    out_dir_date = strmid( out_dir, pos + 1, 10 )
    ; Replace '-' separators to '.'.
    out_dir_date = red_strreplace( out_dir_date, '-', '.', n = 2 )
  endelse

  ; Check if date is given.
  if ( n_elements( date ) eq 0 ) then begin
    date = ''
  endif else begin
    ; Check if the date format is correct.
    isdateformat = stregex( date,                         $
      '^[12][0-9][0-9][0-9][.-][01][0-9][.-][0-3][0-9]$', $
      /boolean                                            )
    if ( isdateformat ) then begin
      ; Replace '-' separators to '.'.
      date = red_strreplace( date, '-', '.', n = 2 )
    endif else begin
      message, 'Date format is wrong, "' + date + '"'
      retall
    endelse
  endelse

  if ( date eq '' ) then begin
    if ( out_dir_date eq '' ) then begin
      message, 'Either date or out_dir must specify the date.'
      retall
    endif else begin
      ; The date is specified by the out_dir only.
      date = out_dir_date
    endelse
  endif else begin
    if ( out_dir_date eq '' ) then begin
      ; Create a sub-directory in out_dir corresponding to date.
      out_dir = out_dir + date + '/'
      file_mkdir, out_dir
    endif else begin
      ; Both out_dir and date specify the date.  Check if it is the same.
      if ( out_dir_date ne date ) then begin
        message, 'out_dir and date specify different dates.'
        retall
      endif
    endelse
  endelse

  ; Date in ISO format and with dots.
  isodate = red_strreplace( date,    '.', '-', n = 2 )
  date    = red_strreplace( isodate, '-', '.', n = 2 )

  ;goto, jump1

  red_setupworkdir, search_dir = search_dirs, out_dir = out_dir, $
    instruments = instruments, date = date, /calibrations_only

  ;return
;jump1:
  ; The current directory is where you run this scrip.  Push it to the stack and
  ; and go the out_dir directory.
  pushd, out_dir ; cd, out_dir

  ; Go to each instrument directory (if exists) and run all generated scripts
  ; inside one by one.
  for i = 0, n_elements( instruments ) - 1 do begin

    instrument_dir = './' + instruments[ i ] + '-calibrations/'

    ; Check if the directory exists
    if file_test( instrument_dir, /directory ) then begin

      pushd, instrument_dir ;cd, instrument_dir
      message, 'At ' + instrument_dir, /informational

      ; Find all scripts.
      doit_files = file_search( 'doit*.pro', count = nfound_scripts )

      if ( nfound_scripts ge 0 ) then begin

        for f = 0, nfound_scripts - 1 do begin

          message, 'Script file: ' + doit_files[ f ], /informational

          ; Open the script file as a plain text...
          openr, flun, doit_files[ f ], /get_lun
          tstr = ""

          while ~eof( flun ) do begin

            ; ...parse it line-by-line and...
            readf, flun, tstr
            tstr = strtrim( tstr, 2 )

            ; ...execute each line but empty ones and end statements.
            if (    ( tstr ne ''    ) $
                 && ( tstr ne 'end' ) ) then begin

              ; Note: execute() function doesn't work in IDL Virtual Machine.
              ; There must be a possible but more complicated solution using
              ; call_procedure, call_function, and call_method tools.
              result = execute( tstr )
              if ~result then begin
                message, 'Failed execute() at line = ' + tstr, /informational
              endif else begin
                message, 'Executed line: ' + tstr, /informational
              endelse

            endif

          endwhile ; eof( flun )

          free_lun, flun

          ; Destroy the reduction class.  It is always named "a" in doit.pro.
          obj_destroy, a

        endfor ; doit_files[ f ]

      endif ; nfound_scripts > 0

      ; Step one level up and proceed to the next instrument's directory.
      popd ; cd, '..'

    endif ; instrument_dir exists

  endfor ; instruments[ i ]

  popd ; Go from the out_dir back to the current script's directory.

  return
end