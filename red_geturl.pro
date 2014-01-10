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
;  tmp = (strsplit(url,'/',/extract, count = n))
  
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
           ;; Read the existing file
           spawn, 'cat '+path, contents
        endif
        
        ;; Return true (for OK) since the file is there to be used.
        return, 1 

     endif

     ;; At this point we know that the wanted file does not exist
     ;; locally and it needs to be downloaded to disk. 
     
     print, 'red_geturl : Try to download '+url
     
     oUrl = OBJ_NEW('IDLnetUrl' $
                    , URL_SCHEME = urlComponents.scheme $
                    , URL_HOSTNAME = urlComponents.host $
                    , URL_PATH = urlComponents.path $
                    , URL_PORT = urlComponents.port $
                   ) 

     ;; Download to a temporary file name so we do not unnecessarily
     ;; overwrite an existing version.
     tmpfile = String('tmp_', Bin_Date(SysTime()), format='(A, I4, 5I2.2)')
     
     CATCH, Error_status
     if Error_status ne 0 then begin
        
        print, 'Caught an error'
        CATCH, /CANCEL

     endif else begin
        
        retrievedFilePath = oUrl -> Get(FILENAME=tmpfile) 

     endelse


     oUrl -> GetProperty, RESPONSE_CODE=RespCode ; 200 = OK

     ;; Beware that some web servers return 200 in spite of failure! The
     ;; sst server seems to behave, though.
     DownloadOK = RespCode eq 200 ; True if OK

     if DownloadOK then begin
        
        file_move, tmpfile, path, /overwrite
        print, 'red_geturl : Downloaded OK to ' + path
        if n_elements(link) ne 0 then begin
           file_link, path, link
           print, 'red_geturl : Linked to ' + link
        endif

        if n_elements(link) ne 0 then begin
           ;; Link anyway
           file_link, path, link
           print, 'red_geturl : Linked to ' + link
        endif
        
        if arg_present(contents) then begin
           ;; Read the existing file
           spawn, 'cat '+path, contents
        endif

     endif else begin

        path = ''
        file_delete, tmpfile
        if n_elements(respcode) eq 0 then begin
           print, 'red_geturl : Download failed.'
        endif else begin
           print, 'red_geturl : Download failed with code '+strtrim(respcode, 2)
        endelse
        ;; Should we make a link if there is an already existing file?
        
     endelse
     
     oUrl->CloseConnections 

     OBJ_DESTROY, oUrl 

     return, DownloadOK         ; True if OK

  endif else begin

     ;; Do not involve the disk in any way.
     
     print, 'red_geturl : Try to download '+url
     
     oUrl = OBJ_NEW('IDLnetUrl' $
                    , URL_SCHEME = urlComponents.scheme $
                    , URL_HOSTNAME = urlComponents.host $
                    , URL_PATH = urlComponents.path $
                    , URL_PORT = urlComponents.port $
                   ) 

     CATCH, Error_status
     if Error_status ne 0 then begin
        
        print, 'Caught an error'
        CATCH, /CANCEL

     endif else begin
        
        contents = oUrl -> Get(/STRING_ARRAY) 

     endelse

     oUrl -> GetProperty, RESPONSE_CODE=RespCode ; 200 = OK

     ;; Beware that some web servers return 200 in spite of failure! The
     ;; sst server seems to behave, though.
     DownloadOK = RespCode eq 200 ; True if OK

     oUrl -> CloseConnections 

     OBJ_DESTROY, oUrl 

     return, DownloadOK         ; True if OK

  endelse

end
