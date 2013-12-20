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
;    file : in, optional, type=string
;
;      The local file name under which to store the downloaded file.
;      Default is to get the file name from the url.
;
;    link : in, optional, type=string
;
;      If this keyword is present, then make a softlink with this name
;      to the downloaded file.
;
;
;    dir : in, optional, type=string
;
;      The local directory name under which to store the downloaded
;      file. Default is the current directory.
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
; :History:
; 
;     2013-12-19 : MGL. Make it do soft links. Use idl function
;                  parse_url() instead of doing my own parsing. Make
;                  the link also if the file already existed.
;
;     2013-12-20 : MGL. Remove link if download did not work (and file
;                  did not already exist). New keyword "path".
;
;-
function red_geturl, url, file = file, dir = dir, overwrite = overwrite, link = link, path = path

  if n_elements(dir) eq 0 then dir = './' else file_mkdir, dir ; Make the directory just in case.

  urlComponents = parse_url(url)

  tmp = (strsplit(url,'/',/extract, count = n))
  if n_elements(file) eq 0 then begin
     file = (strsplit(urlComponents.path,'/',/extract, count = n))[n-1]
  endif

  if n_elements(link) ne 0 then begin
     ;; Delete any existing link. If the download works we will create
     ;; a new link, if it doesn't it should not be there.
     file_delete, link, /allow_nonexistent
  endif

  if ~keyword_set(overwrite) and file_test(dir+file) then begin
     print, 'red_geturl : Do not download, file already exists (or use /overwrite):'
     print, '             '+url
     print, '             '+dir+file
     if n_elements(link) ne 0 then begin
        ;; Link anyway
        file_link, dir+file, link
        print, 'red_geturl : Linked to ' + link
     endif
     ;; Return true (for OK) anyway since the file is there to be
     ;; used:
     return, 1 
  endif

  print, 'red_geturl : Try to download '+url
  
  oUrl = OBJ_NEW('IDLnetUrl' $
                 , URL_SCHEME = urlComponents.scheme $
                 , URL_HOSTNAME = urlComponents.host $
                 , URL_PATH = urlComponents.path $
                 , URL_PORT = urlComponents.port $
                ) 

  tmpfile = String('tmp_', Bin_Date(SysTime()), format='(A, I4, 5I2.2)')
  
  retrievedFilePath = oUrl->Get(FILENAME=tmpfile) 
  oUrl->GetProperty, RESPONSE_CODE=RespCode ; 200 = OK

  DownloadOK = RespCode eq 200  ; True if OK

  if DownloadOK then begin
     path = dir+file
     file_move, tmpfile, dir+file, /overwrite
     print, 'red_geturl : Downloaded OK to ' + dir+file
     if n_elements(link) ne 0 then begin
        file_link, dir+file, link
        print, 'red_geturl : Linked to ' + link
     endif
  endif else begin
     path = ''
     file_delete, tmpfile
     print, 'red_geturl : Download failed with code '+respcode
  endelse

  oUrl->CloseConnections 

  OBJ_DESTROY, oUrl 

  return, DownloadOK            ; True if OK

end
