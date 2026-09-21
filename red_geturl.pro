; docformat = 'rst'

;+
; Downloads a file specified by a url. 
; 
; :Categories:
;
;    SST observations
; 
; 
; :Author:
; 
;    Mats Löfdahl, 2013-12-19
; 
; :Returns:
; 
;    Returns a boolean value, true if download worked, false if it didn't.
; 
; :Params:
; 
;    url : in, type=string
; 
;      The url corresponding to the file to be downloaded.
; 
; 
; :Keywords:
;
;    file : in, optional, type=string, default="Get from url"
;
;      The local file name under which to store the downloaded file.
;      (But see also the CONTENTS keyword.)
;
;    link : in, optional, type=string
;
;      If this keyword is present, then make a softlink with this name
;      to the downloaded file.
;
;
;    dir : in, optional, type=string, default=PWD
;
;      The local directory name under which to store the downloaded
;      file. (But see also the CONTENTS keyword.)
;
;    overwrite : in, optional, type=boolean
;
;      Set this to download without checking if the file already
;      exists. 
;
;    path : out, optional, type=string
;
;      The path to the downloaded file (or the empty string if the
;      download failed).
;
;    contents : out, optional,type=strarr
;
;      The contents of the downloaded file. If any of keywords LINK,
;      PATH, DIR, or FILE is present, a download will only happen if
;      the file does not already exist locally or if OVERWRITE is set.
;
; 
; :History:
; 
;     2013-12-19 : MGL. Make it do soft links. Use idl function
;                  parse_url() instead of doing my own parsing. Make
;                  the link also if the file already existed.
;
;     2013-12-20 : MGL. Remove link if download did not work (and file
;                  did not already exist). New keyword "path". Delete
;                  any old versions of the files we try but fail to
;                  download. 
;
;     2014-01-02 : MGL. New keyword "contents". Some new logic for
;                  when to save a file to disk and to avoid reading
;                  the contents into IDL unless necessary.
;
;     2014-01-10 : MGL. Add some error handling.
; 
;     2017-08-18 : THI. Workaround so that paths/filenames with ":"
;                  does not fail. Get rid of warning when download fails.
; 
;     2026-09-21 : MGL. Use SPAWN + curl to bypass IDLnetURL SSL
;                  issues (Error 60) on pre-IDL 9.
;
;-
function red_geturl, url $
                     , file = file $
                     , dir = dir $
                     , overwrite = overwrite $
                     , link = link $
                     , path = path $
                     , contents = contents

  ;; Should disk I/O be involved?
  DiskIO = (n_elements(file) ne 0 $
            or n_elements(link) ne 0 $
            or n_elements(dir) ne 0 $
            or arg_present(path) $
            or ~arg_present(contents) $
           )
  
  urlComponents = parse_url(url)
  
  ; parse_url is broken, it can not handle file/dirnames with ':' as it will be interpreted as a
  ; port-separator. The hack below tries to fix the mess.
  slash_pos = strpos( urlComponents.Host, '/' )
  if slash_pos ne -1 then begin
    urlComponents.Host = strmid( urlComponents.Host, 0, slash_pos )
    port_pos = strpos( urlComponents.Host, ':' )    ; did the original url conain a port?
    if port_pos ne -1 then begin
      host_len = strlen( urlComponents.Host )
      urlComponents.Port = strmid( urlComponents.Host, port_pos+1, host_len-(port_pos+1))
    endif else begin
      urlComponents.Port = '80'     ; re-set the default, as it was most likely set to '', or garbage, by parse_url
    endelse
  endif

  ;; Check OS family to set correct SPAWN flags later
  isWindows = (!VERSION.OS_FAMILY eq 'Windows')

  if DiskIO then begin

     ;; Download to local file if necessary, then read the local file
     ;; into CONTENTS if this keyword is present.
     
     if n_elements(dir) eq 0 then begin
        dir = './'              ; Default dir is PWD
     endif else begin
        file_mkdir, dir         ; Make the directory just in case.
     endelse

     ;; Default file name taken from url.
     if n_elements(file) eq 0 then begin
        file = (strsplit(urlComponents.path,'/',/extract, count = n))[n-1]
     endif

     ;; This is where the file is to be stored/read.
     path = dir+file

     if n_elements(link) ne 0 then begin
        ;; Delete any existing link. If the download works we will create
        ;; a new link, if it doesn't it should not be there.
        file_delete, link, /allow_nonexistent
     endif

     if file_test(path) and ~keyword_set(overwrite) then begin

        print, 'red_geturl : Do not download, file already exists (or use /overwrite):'
        print, '             '+url
        print, '             '+path

        if n_elements(link) ne 0 then begin
           ;; Link anyway
           file_link, path, link
           print, 'red_geturl : Linked to ' + link
        endif

        if arg_present(contents) then begin
           ;; Read the existing file (cross-platform compatible method)
           openr, lun, path, /get_lun
           file_info = file_info(path)
           contents = strarr(file_info.size) ; fallback if text, though image binaries shouldn't use contents
           readf, lun, contents
           free_lun, lun
        endif
        
        ;; Return true (for OK) since the file is there to be used.
        return, 1 

     endif

     ;; At this point we know that the wanted file does not exist
     ;; locally and it needs to be downloaded to disk. 
     
     print, 'red_geturl : Try to download '+url
     
     ;; Download to a temporary file name so we do not unnecessarily
     ;; overwrite an existing version.
     tmpfile = String('tmp_', Bin_Date(SysTime()), format='(A, I4, 5I2.2)')
     
     ;; Build the curl command: 
     ;; -s (silent), -L (follow redirects), -f (fail silently on server errors like 404)
     ;;cmd = 'curl -s -L -f -o "' + tmpfile + '" "' + url + '"'
     cmd = (isWindows ? '' : 'env LD_LIBRARY_PATH="" ') + 'curl -s -L -f -o "' + tmpfile + '" "' + url + '"'

     curl_status = 0
     if (isWindows) then begin
         SPAWN, cmd, /NOSHELL, EXIT_STATUS=curl_status
     endif else begin
         SPAWN, cmd, EXIT_STATUS=curl_status
     endelse

     ;; curl returns 0 upon successful download [1]
     DownloadOK = (curl_status eq 0) and file_test(tmpfile)

     if DownloadOK then begin
        
       file_move, tmpfile, path, /overwrite
       print, 'red_geturl : Downloaded OK to ' + path

       if n_elements(link) ne 0 then begin
         file_link, path, link
         print, 'red_geturl : Linked to ' + link
       endif
        
       if arg_present(contents) then begin
          spawn, (isWindows ? 'type "' : 'cat "') + path + '"', contents
       endif

     endif else begin

        path = ''
        file_delete, tmpfile, /ALLOW_NONEXISTENT
        print, 'red_geturl : Download failed with curl exit status '+strtrim(curl_status, 2)
        
     endelse
     
     return, DownloadOK         ; True if OK

  endif else begin

     ;; Do not involve the disk in any way. Strarr output directly.
     print, 'red_geturl : Try to download '+url
     
     ;;cmd = 'curl -s -L -f "' + url + '"'
     cmd = (isWindows ? '' : 'env LD_LIBRARY_PATH="" ') + 'curl -s -L -f "' + url + '"'
     curl_status = 0

     if (isWindows) then begin
         SPAWN, cmd, contents, /NOSHELL, EXIT_STATUS=curl_status
     endif else begin
         SPAWN, cmd, contents, EXIT_STATUS=curl_status
     endelse

     DownloadOK = (curl_status eq 0)
     if ~DownloadOK then begin
         print, 'red_geturl : Download failed with curl exit status '+strtrim(curl_status, 2)
     endif

     return, DownloadOK         ; True if OK

  endelse

end
