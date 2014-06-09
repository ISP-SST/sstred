; docformat = 'rst'

;+
; Find directories that match a regular expression.
;
; When a directory matches, do not search in its subdirectories.
; 
; :Categories:
;
;    SST observations
; 
; 
; :author:
; 
;    Mats Löfdahl, 2014-06-09
;
; :Keywords:
;
;
;
; :History:
; 
;     2014-06-09 : MGL. First version.
; 
;-
function red_find_matching_dirs, regex, rootdir = rootdir, maxdepth = maxdepth, count = count, verbose = verbose

  if n_elements(rootdir) eq 0 then rootdir = './'
  if n_elements(maxdepth) eq 0 then maxdepth = 10000

  if strmid(rootdir,strlen(rootdir)-1) ne '/' then rootdir += '/'

  if keyword_set(verbose) then print, 'Rootdir = ', rootdir
  if keyword_set(verbose) then print, 'regex = ', regex

  if strmatch(rootdir, '*'+regex+'*') then begin
     count = 1
     return, rootdir
  endif

  if maxdepth gt 0 then begin

     ;; Find all subdirectories
     dnames = file_search(rootdir+'*', count = Nd, /test_dir)
     
     if keyword_set(verbose) then print, 'dnames : ', dnames
     
     founddirs = ['']
  
     for i = 0, Nd-1 do begin
        if keyword_set(verbose) then print
;        if strmatch(dnames[i], '*'+regex+'*') then begin
;           if keyword_set(verbose) then print, 'Match : '+dnames[i]
;           founddirs = [founddirs, dnames[i]]
;        endif else begin
           if keyword_set(verbose) then print, 'Descend into :'+dnames[i]
           snames = red_find_matching_dirs(regex, rootdir = dnames[i], maxdepth = maxdepth-1 $
                                           , count = Nsubdirs, verbose = keyword_set(verbose))
           if Nsubdirs gt 0 then founddirs = [founddirs, snames]
;        endelse
        
     endfor                     ; i
  endif                         ; maxdepth

  count = n_elements(founddirs)
 
  if count gt 1 then begin
     count += -1
     return, founddirs[1:*]
  endif else begin
     count = 0
     return, ''
  endelse

end


