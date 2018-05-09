; docformat = 'rst'

;+
; Collection and basic processing of calibration data from CHROMIS and
; CRISP observations.
;
; :Categories:
;
;   SST pipeline
;
; :Author:
;
;   Andrii V. Sukhorukov, Institute for Solar Physics, asukh@astro.su.se
;
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
;                 default = '/storage/sand*/... in Stockholm or /data/disk?/...
;                 on La Palma'.
;
;   restart     : in, optional, type = boolean,
;                 default = False
;
; :History:
;
;   2017-07-17 : AVS. Initial commit.
;
;   2017-07-24 : AVS. TODO: tested calling IDL from inside, doesn't work as
;                     expected.
;
;   2017-07-27 : AVS. A line-by-line parser is added.
;
;   2017-07-30 : AVS. TODO: doesn't work on La Palma under the old
;                CRISPRED installation. Calling hostname -I is not
;                working too.
;
;   2017-07-31 : AVS. Keyword and directory handling is completely
;                rewritten. 
;
;   2017-08-01 : AVS. The IP-address logic, searching for the
;                instrumental directories, and the script generation
;                using a common call_procedure invocation are
;                rewritten.
;
;   2017-08-02 : AVS. The script parser is updated. Succesfully tested
;                on recent data sets in /storage/sand??/Incoming/
;
;   2017-08-11 : AVS. Added search_dirs keyword to be able to run in
;                different locations where old data is stored.
;
;   2017-08-14 : AVS. Changed cd calls to pushd/popd for a better
;                portability.
;
;   2017-08-24 : AVS. Logic in internal loops have been moved to
;                red_setupworkdir. Instead call red_setupworkdir with
;                the proper keywords and the doit.pro/config.txt pairs
;                will be created outside this script.
;
;   2017-08-28 : AVS. Logging of done (successfully executed) lines
;                from doit*.pro scripts to done*.log logs is added.
;
;   2017-08-30 : AVS. a) Default values of out_dir are provided. b)
;                Use current date if date or out_dir do not specify
;                it.
;
;   2017-09-11 : AVS. A new keyword, /restart, is added to avoid
;                re-creating existing script/config files in case of a
;                restart.
;
;-
pro red_getcalibdata,         $
   date        = date,        $
   out_dir     = out_dir,     $
   search_dirs = search_dirs, $
   restart     = restart      
  
  ;; Instruments to process data for. TRIPPEL and SLITJAW are not
  ;; included yet.
  instruments = [ 'CHROMIS', 'CRISP' ]
  
  ;; The default calling sequence is to specify the output directory
  ;; in out_dir and it should contain the date as the last
  ;; subdirectory, for example, /scratch/asukh/2017.07.15.

  ;; If out_dir is not given, set it to some commonly-agreed place
  ;; otherwise use the current directory.
  if ( n_elements( out_dir ) eq 0 ) then begin
    case ( get_login_info( ) ).machine_name of
      'camera4' : begin
          ;; This is the default storage for camera4, it is the NFS mounted
          ;; directory of same name on transport1 for easy transfers
        out_dir = '/scratch/calibration-data'
      end
      'transport1' : begin
        ;; This is a local (=fast) SDD drive that should be used to transfer
        ;; processed calibration data from La Palma to Stockholm.
        ;; It is available from outside as royac6::calib/ for rsync  
        out_dir = '/scratch/calibration-data'
      end
      'freija' : begin
        ;; This my temporary folder.
        out_dir = '/scratch/asukh/calibration-data'
      end
      else : begin
        out_dir = getenv( 'PWD' )
        message, 'out_dir is not given and is set to the current working ' + $
                 'directory, ' + out_dir, $
                 /informational
      end
    endcase
  endif

  ;; out_dir must end with a slash.
  if ~strmatch( out_dir, '*/' ) then out_dir += '/'

  ;; Search if there is a date at the end of out_dir.
  pos = stregex( out_dir,                                  $
                 '\/[12][0-9][0-9][0-9][.-][01][0-9][.-][0-3][0-9]\/?$' )
  if ( pos eq -1 ) then begin
    out_dir_date = ''
  endif else begin
    out_dir_date = strmid( out_dir, pos + 1, 10 )
    ;; Replace '-' separators to '.'.
    out_dir_date = red_strreplace( out_dir_date, '-', '.', n = 2 )
  endelse

  ;; Check if date is given.
  if ( n_elements( date ) eq 0 ) then begin
    date = ''
  endif else begin
    ;; Check if the date format is correct.
    isdateformat = stregex( date,                         $
                            '^[12][0-9][0-9][0-9][.-][01][0-9][.-][0-3][0-9]$', $
                            /boolean                                            )
    if ( isdateformat ) then begin
      ;; Replace '-' separators to '.'.
      date = red_strreplace( date, '-', '.', n = 2 )
    endif else begin
      message, 'date format is wrong, "' + date + '"'
      retall
    endelse
  endelse

  if ( date eq '' ) then begin
    if ( out_dir_date eq '' ) then begin
      message, 'neither date nor out_dir specify the date, so, the current ' + $
               'date is taken.', /informational
      date = strmid( red_timestamp( /utc, /iso ), 0, 10 )
    endif else begin
      ;; The date is specified by the out_dir only.
      date = out_dir_date
    endelse
  endif else begin
    if ( out_dir_date eq '' ) then begin
      ;; Create a sub-directory in out_dir corresponding to date.
      out_dir = out_dir + date + '/'
      file_mkdir, out_dir
    endif else begin
      ;; Both out_dir and date specify the date. Check if it is the
      ;; same.
      if ( out_dir_date ne date ) then begin
        message, 'out_dir and date specify different dates.'
        retall
      endif
    endelse
  endelse

  ;; Date in ISO format and with dots.
  isodate = red_strreplace( date,    '.', '-', n = 2 )
  date    = red_strreplace( isodate, '-', '.', n = 2 )
  
  ;; If /restart is set then the instrument directories with
  ;; corresponding script and config files inside will not be updated
  ;; and will be kept as is.
  if ~keyword_set( restart ) then begin
    red_setupworkdir, search_dir = search_dirs, out_dir = out_dir, $
                      instruments = instruments, date = date, /calibrations_only
  endif
  
  ;; The current directory is where you run this script. Push it to
  ;; the directory stack and go the out_dir directory.
  pushd, out_dir

  ;; Go to each instrument directory (if exists) and run all generated
  ;; scripts inside one by one.
  for i = 0, n_elements( instruments ) - 1 do begin

    instrument_dir = './' + instruments[ i ] + '-calibrations/'

    ;; Check if the directory exists
    if file_test( instrument_dir, /directory ) then begin

      pushd, instrument_dir
      message, 'at ' + instrument_dir, /informational
      
      ;; Find all scripts.
      doit_files = file_search( 'doit*.pro', count = nfound_scripts )

      if ( nfound_scripts gt 0 ) then begin

        for f = 0, nfound_scripts - 1 do begin

          doit_file = doit_files[ f ]
          message, 'script file: ' + doit_file, /informational

          ;; Read the script file line by line in an array of strings.
          openr, doit_lun, doit_file, /get_lun
          doit_str = ""
          while ~eof( doit_lun ) do begin
            readf, doit_lun, doit_str
            doit_str = strtrim( doit_str, 2 )
            if ( doit_str ne '' ) then red_append, doit_lines, doit_str
          endwhile
          free_lun, doit_lun

          ;; The done*.log file has the same suffix as doit*.pro.
          dotpos = strpos( doit_file, '.pro', /reverse_search )
          ;; No need to check for match as the doit*.pro name is already matched.
          done_file = doit_file
          strput, done_file, 'done', 0      ; doit*.pro -> done*.pro
          strput, done_file, '.log', dotpos ; done*.pro -> done*.log

          if ( file_test( done_file ) ) then begin
            ;; Read the done*.log line by line into an array of strings.
            openr, done_lun, done_file, /get_lun
            done_str = ""
            while ~eof( done_lun ) do begin
              readf, done_lun, done_str
              done_str = strtrim( done_str, 2 )
              if ( done_str ne '' ) then red_append, done_lines, done_str
            endwhile
            free_lun, done_lun
          endif else begin
            ;; No done*.log file, nothing has been processed before.
            done_lines = [ ]
          endelse
          ;; Create the done*.log file anew by opening it for writing.
          openw, done_lun, done_file, /get_lun
          free_lun, done_lun

          n_doit_lines = n_elements( doit_lines )
          n_done_lines = n_elements( done_lines )

          for iline = 0, n_doit_lines - 1 do begin
            line = doit_lines[ iline ]

            ;; Execute each line but the end statement.
            if ( line ne 'end' ) then begin

              ;; Search if the line has been completed before.
              if ( n_done_lines gt 0 ) then begin
                is_done = -1
                for idoneline = 0, n_done_lines - 1 do begin
                  is_done = strpos( done_lines[ idoneline ], line )
                  if ( is_done ge 0 ) then begin
                    message, 'line "' + line + '" has been completed before.', $
                             /informational
                    ;; Append the old completed line to the new done*.log file.
                    openw, done_lun, done_file, width = 256, /append, /get_lun
                    printf, done_lun, done_lines[ idoneline ]
                    free_lun, done_lun
                                ; Exit the for-loop over done lines.
                    break
                  endif
                endfor                            ; idoneline
                if ( is_done ge 0 ) then continue ; Next doit line.
              endif

              ;; The line has not been executed before. Run it and
              ;; write to the done*.log file with the corresponding
              ;; UTC time-stamp. Warning: execute() function does not
              ;; work in IDL Virtual Machine.
              result = execute( line )
              if ~result then begin
                message, 'execute("' + line + '") failed.', /informational
              endif else begin
                message, 'execute("' + line + '") succeeded.', /informational
                ;; Only method calls such as a->sumdark must be saved.
                is_redclass_method = stregex( line,         $
                                              '^a *-> *sum(dark|flat|pinh|polcal) *,',  $
                                              /boolean ) ;
                if is_redclass_method then begin
                  ;; Append the completed line to the done*.log file.
                  openw, done_lun, done_file, width = 256, /append, /get_lun
                  printf, done_lun, line + ' ; ' + red_timestamp( /iso, /utc )
                  free_lun, done_lun
                endif
              endelse
            endif               ; line ne 'end'
          endfor                ; line = doit_lines[ iline ]

          ;; Destroy the reduction class.  It is always named "a" in doit.pro.
          obj_destroy, a
          ;; Nullify the doit_lines and the done_lines arrays (only in IDL > v8.*).
          doit_lines = [ ]
          done_lines = [ ]

        endfor                  ; doit_file = doit_files[ f ]
        
      endif                     ; nfound_scripts > 0
      
      ;; Step one level up and proceed to the next instrument's
      ;; directory.
      popd                      ; cd, '..'

    endif                       ; instrument_dir exists

  endfor                        ; instruments[ i ]

  popd                          ; Go from the out_dir back to the current script's directory.

  return
end
