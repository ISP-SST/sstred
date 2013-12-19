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
;      Set this to download without checking if the file already exists
;
; :History:
; 
;     2013-12-19 : MGL. Make it do soft links. Use idl function
;                  parse_url() instead of doing my own parsing.
;
;-
function red_geturl, url, file = file, dir = dir, overwrite = overwrite, link = link

  if n_elements(dir) eq 0 then dir = './' else file_mkdir, dir ; Make the directory just in case.

  urlComponents = parse_url(url)

;  url_scheme = (strsplit(url, ':',/extract))[0]                 ; http, ftp, etc.
  tmp = (strsplit(url,'/',/extract, count = n))
;  url_hostname = tmp[1]
;  url_path = strjoin(tmp[2:*], '/')
  if n_elements(file) eq 0 then begin
     file = (strsplit(urlComponents.path,'/',/extract, count = n))[n-1]
  endif

  if ~keyword_set(overwrite) and file_test(dir+file) then begin
     print, 'red_geturl : Do not download, file already exists (or use /overwrite):'
     print, '             '+url
     print, '             '+dir+file
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
     file_move, tmpfile, dir+file, /overwrite
     print, 'red_geturl : Downloaded OK to ' + dir+file
     if n_elements(link) ne 0 then begin
        file_delete, link, /allow_nonexistent
        file_link, dir+file, link
        print, 'red_geturl : Linked to ' + link
     endif
  endif else begin
     file_delete, tmpfile
     print, 'red_geturl : Download failed with code '+respcode
  endelse

  oUrl->CloseConnections 

  OBJ_DESTROY, oUrl 

  return, DownloadOK            ; True if OK

end
