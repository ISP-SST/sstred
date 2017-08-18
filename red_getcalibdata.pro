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
;-
pro red_getcalibdata,       $
  date        = date,       $
  out_dir     = out_dir,    $
  search_dirs = search_dirs ;

  ; Instruments to process data for.  TRIPPEL and SLITJAW are not included yet.
  instruments = [ 'CHROMIS', 'CRISP' ]

  ; The telescope geographic location.
  obsgeo_xyz = round( red_obsgeo( 28.759733d, -17.880736d, 2360d ) )

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

  ; If search_dirs are not given, pick them depending on the current computer
  ; location.
  if ( n_elements( search_dirs ) eq 0 ) then begin
    message, 'search_dirs are not given and are set depending on the ' + $
      'current computer location.', $
      /informational

    ; Which computer you are running this at the moment?
    ;
    ; Old 'hostname -I' doesn't work on La Palma under an old OpenSuSE
    ; installation therefore the only supported keyword there is -i that gives
    ; the right IP-address on La Palma but not in Stockholm.  There are more
    ; robust but less trivial ways to get the IPv4-address of the machine.  As
    ; for 31 July 2017, the following works:
    ;   $ ip addr show | grep -Po 'inet \K[\d.]+'
    ; or
    ;   $ /sbin/ifconfig -a | grep -Po 't addr:\K[\d.]+'
    ; Both produce a list of IPv4-addresses for each network interface on La
    ; Palma as well as in Stockholm.  You can narrow the list to only one
    ; address by specifying 'ip addr show eth0' or '/sbin/ifconfig eth0' but it
    ; is better to use full listing as some interfaces might be temporarily off
    ; the network.
    spawn, "ip addr show | grep -Po 'inet \K[\d.]+'", ipv4addresses

    ; Check if any of the addressess matches the local range of IPv4-addresses.
    case 1 of
      max( strmatch( ipv4addresses, '*161.72.15.*' ) ) : begin
        ; SST network address range on La Palma.
        ; Observed data is in any folder mounted at /data/, either disk? or
        ; camera?, or some other drive.  As Pit said on 28 July 2017, we should
        ; search /data/camera? folders only if it is really necessary.  Observed
        ; data must be frist transfered to /data/disk? folders, by default to
        ; disk1.  The inner directory level is either CHROMIS folder for CHROMIS
        ; data or something different, usually the name of the insitute like
        ; ISP, UK, UIO, for CRISP data.  For example,
        ;   /data/disk1/CHROMIS/2017.04.05
        ;   /data/disk1/ISP/2017.04.05
        ;search_dirs = "/data/*/*/" ; Full contents of /data folder.
        search_dirs = "/data/disk?/*/" ; Only /data/disk? storages.
      end
      max( strmatch( ipv4addresses, '*130.237.166.*' ) ) : begin
        ; ISP network range at AlbaNova in Stockholm.
        ; Transfered data is stored in all the sandboxes at /storage/sand*.
        ; The inner directory level is either empty, or Incoming, or
        ; Incoming/Checked.  For example:
        ;   /storage/sand04n/2017.04.05
        ;   /storage/sand05/Incoming/2017.04.05
        search_dirs = '/storage/sand*/' + ['', 'Incoming/', 'Incoming/Checked/']
      end
      else : begin
        message, 'No matching IPv4-address in ' + strjoin( ipv4addresses, ', ' )
        retall
      end
    endcase

  endif ; no search_dirs are given

  for i = 0, n_elements( search_dirs ) - 1 do begin

    ; Each search directory must end with a slash.
    if ~strmatch( search_dirs[ i ], '*/' ) then search_dirs[ i ] += '/'

    ; Include the date directory at the end, if it is not added yet.
    if ( file_basename( search_dirs[ i ] ) ne date ) then begin
      search_dirs[ i ] += date
    endif

  endfor

  ; Search all directories.  This might be slow for non-mounted camera? or
  ; disk? folders.
  found_dirs = file_search( search_dirs, count = nfound_dirs )

  ; Remove duplicates.
  if nfound_dirs gt 1 then begin
    found_dirs  = found_dirs[ uniq( found_dirs, sort( found_dirs ) ) ]
    nfound_dirs = n_elements( found_dirs )
  endif

  if ( nfound_dirs eq 0 ) then begin
    print, 'No data is found for ' + date + '.'
    return
  endif

  ; Loop over all found directories.
  for i = 0, nfound_dirs - 1 do begin

    ; For each found directory, a doit.pro script and a config.txt file are
    ; generated.  They differ by the number offset such as doit02.pro etc.
    found_dir = found_dirs[ i ]

    ; Each found dir must end with a slash.
    if ~strmatch( found_dir, '*/' ) then found_dir += '/'

    ; Instrument directories are named Chromis-D, Chromis-N, and Chromis-W for
    ; the CHROMIS data and Crisp-R, Crisp-T, and Crisp-W for the CRISP data.
    ; They are nested at the third level from the root of each found directory.
    ; For example,
    ;   .../2017.05.02 /CHROMIS-flats /8545     /Chromis-W
    ;   .../2017.05.02 /Darks         /12:09:18 /Crisp-T
    ; where .../2017.05.02 is the current found directory.
    instrument_dirs = file_search( found_dir + '*/*/*', $
      /test_directory,                                  $
      count = ninstrument_dirs                          )

    if ( ninstrument_dirs eq 0 ) then begin
      message, 'No instrument directories found in ' + found_dir, /informational
      continue
    endif else begin
      message, 'Data found in ' + found_dir, /informational
    endelse

    ; Run over all instruments and ...
    for inst = 0, n_elements( instruments ) - 1 do begin

      instrument = instruments[ inst ]

      ; Check if this instrument is present in the instrument directories.  An
      ; instrument diretory is somethin like
      ;   /storage/sand04/Incoming/2017.07.01/CHROMIS-darks/14:24:47/Chromis-D
      ; that ends with /<instrument>-<[DWT],[RTW]>.
      if ( max( strmatch( instrument_dirs, $
                  '*\/' + instrument + '-*', /fold_case ) ) gt 0 ) then begin

        ; A work dir is, for example, .../CRISP-calibrations/.
        work_dir = out_dir + instrument + '-calibrations/'

        ; If work_dir doesn't exist, create it and put the preliminary metadata
        ; inside.
        if ~file_test( work_dir, /directory ) then begin

          file_mkdir, work_dir

          ; Write string metadata
          red_metadata_store, fname = work_dir + '/info/metadata.fits',    $
            [ { keyword : 'OBSRVTRY',                                      $
                value   : 'Observatorio del Roque de los Muchachos (ORM)', $
                comment : 'Name of observatory' },                         $
              { keyword : 'TELESCOP',                                      $
                value   : 'Swedish 1-meter Solar Telescope (SST)',         $
                comment : 'Name of telescope' },                           $
              { keyword : 'OBJECT',                                        $
                value   : 'Sun',                                           $
                comment : '' } ]

          ; Write numerical metadata
          red_metadata_store, fname = work_dir + '/info/metadata.fits',    $
            [ { keyword : 'OBSGEO-Z',                                      $
                value   : obsgeo_xyz[ 2 ],                                 $
                comment : '[m] SST location' },                            $
              { keyword : 'OBSGEO-Y',                                      $
                value   : obsgeo_xyz[ 1 ],                                 $
                comment : '[m] SST location' },                            $
              { keyword : 'OBSGEO-X',                                      $
                value   : obsgeo_xyz[ 0 ],                                 $
                comment : '[m] SST location' } ]

        endif ; if work_dir doesn't exist

        ; Set the config file name and the script file name.  Each name has a
        ; two-figure extension, e.g., doit00.pro to separate scripts for
        ; different found folders.
        config_file = string( 'config', i, '.txt', format = '( a, i02, a )' )
        script_file = string( 'doit',   i, '.pro', format = '( a, i02, a )' )

        ; A generic interfaces to different instrumental subroutines.
        call_procedure, 'red_setupworkdir_' + instrument,         $
          work_dir, found_dir, config_file, script_file, isodate, $
          calibrations_only = 1                                   ;
        print, instrument + ' setup in ' + work_dir

      endif ; instrument in instrument_dirs

    endfor ; instruments[ inst ]

  endfor ; found_dirs[ i ]

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
      doit_files = file_search( 'doit??.pro', count = nfound_scripts )

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