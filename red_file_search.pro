; docformat = 'rst'

;+
; Search for files using the shell "find" command, optionally on a
; remote node.
;
; Can be used instead of file_search() for simple searches. Should be
; faster for directories with many files.
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
; :Returns:
; 
;   Returns all matched filenames in a string array, one file name per
;   array element. If no files exist with names matching the input
;   arguments, a null scalar string is returned instead of a string
;   array.
;
; :Params:
; 
;    dir : in, optional, type=string
; 
;      The (local) directory in which to do the search.
;
;    pattern : in, type=string
;
;      The pattern "find" should use to filter the file names. If this
;      includes the directory, the dir parameter is not needed.
;
; :Keywords:
; 
;   remote_dir : in, optional, type=string, default=dir
; 
;     Set this to the remote node's equivalent of the dir parameter. 
; 
;   remote_login : in, optional, type="boolean or string"
; 
;     If true, then use ssh to do the find operation on a remote node.
;     If a string, use it as the node name (or as user@nodename).
; 
;   verbose : in, optional, type=boolean
; 
;     Be verbose!
; 
; :History:
;  
;    2018-07-10 : MGL. First version, based on code from
;                 red__quicklook movie.pro.
;  
;    2018-07-26 : MGL. Possible to optionally separate directory from
;                 the pattern.
; 
;-
function red_file_search, searchstring, dir $
                          , count = count $
                          , remote_dir = remote_dir $
                          , remote_login = remote_login $
                          , verbose = verbose

  inam = red_subprogram(/low, calling = inam1)

  if n_elements(dir) eq 0 then begin
    dir = file_dirname(searchstring)
    pattern = file_basename(searchstring)
  endif else begin
    pattern = searchstring
  end

  case n_elements(pattern) of
    0    : name_expression = "-name '*'"
    1    : name_expression = "-name '"+pattern+"'"
    else : name_expression = "\( " $
                             + strjoin("-name '"+pattern+"'", ' -o ') $
                             + " \)"
  endcase

  print, inam + ' : ' + name_expression
  
  if keyword_set(remote_login) then begin
    
    if size(remote_login, /tname) ne 'STRING' then remote_login = 'root@transport1'
    if n_elements(remote_dir) eq 0 then remote_dir = dir
    
    cmd = 'ssh ' + remote_login $
          + ' "find ' + remote_dir + '/ -type f ' $
          + name_expression+'"'
    
    if keyword_set(verbose) then print, cmd
    spawn, cmd, files

    if remote_dir ne dir then files = dir + '/' + file_basename(files)
    
  endif else begin

    spawn, 'find ' + dir + '/ -type f ' + name_expression, files

  endelse

  if size(files, /n_dim) eq 0 then count = 0 else count = n_elements(files)
  
  return, files
  
end


;    ;; Case 1, called as /ssh_find: transport1, so identical names
;    if n_elements(ssh_find) eq 1 then begin
;      spawn, 'ssh root@transport1 find '+dirs[iset]+'/'+cam+' -type f -name '+pattern, files
;    endif else begin
;      ;; Case 2, called as  ssh_find=['user@host','<local data dir>']
;      if keyword_set(verbose) then $
;         print, 'ssh '+ssh_find(0)+' find '+ssh_find(1)+'/'+cam+' -type f -name '+pattern
;      spawn, 'ssh '+ssh_find(0)+' find '+ssh_find(1)+'/'+cam+' -type f -name '+pattern, files
;      files = strmid(files, strlen(ssh_find(1)))
;      files = dir+files
;    endelse
