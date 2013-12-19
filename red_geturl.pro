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
;    dir : in, optional, type=string
;
;      The local directory name under which to store the downloaded
;      file. Default is the current directory.
;
; :History:
; 
; 
; 
;
;
;
;
;
;
;-
function red_geturl, url, file = file, dir = dir

  if n_elements(dir) eq 0 then dir = './' else file_mkdir, dir ; Make the directory just in case.

  url_scheme = (strsplit(url, ':',/extract))[0]                 ; http, ftp, etc.
  tmp = (strsplit(url,'/',/extract, count = n))
  url_hostname = tmp[1]
  url_path = strjoin(tmp[2:*], '/')
  if n_elements(file) eq 0 then begin
     file = tmp[n-1]
  endif
  
  oUrl = OBJ_NEW('IDLnetUrl', URL_SCHEME = url_scheme, URL_HOSTNAME = url_hostname, URL_PATH = url_path) 

  tmpfile = String('tmp_', Bin_Date(SysTime()), format='(A, I4, 5I2.2)')
  
  retrievedFilePath = oUrl->Get(FILENAME=tmpfile) 
  oUrl->GetProperty, RESPONSE_CODE=RespCode ; 200 = OK
  if RespCode eq 200 then begin
     file_move, tmpfile, dir+file, /overwrite
  endif else begin
     file_delete, tmpfile
  endelse

  oUrl->CloseConnections 

  OBJ_DESTROY, oUrl 

  return, RespCode eq 200       ; True if OK

end
