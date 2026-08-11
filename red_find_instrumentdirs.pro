; docformat = 'rst'

;+
; Look recursively for instrument-specific subdirectory below a given
; directory. 
;
; Subdirectories will be searched for by using regular expressions of
; the following type:
; "topdir+'/'+data_prefix+'*/'+levels+instrument_prefix+'*'", where
; levels consist of an increasing number of repetitions of "*/". Once
; matching subdirectories are found at a certain level, the search
; will not go deeper.
; 
; :Categories:
;
;    SST observations
; 
; 
; :Author:
; 
;    Mats Löfdahl, ISP, 2016-05-05
; 
; 
; :Returns:
; 
;    A string array with directory names.
; 
; 
; :Params:
; 
;    topdir : in, type=string
;
;       Look for directories below this one.
;
;    instrument_prefix : in, type=string
;
;       A string matching the instrument specific directories, like
;       "crisp" or "chromis".
;
;    data_prefix : in, type=string
;   
;       A string matching the data type directory, like "dark" or
;       "flat". 
; 
; 
; 
; 
; :Keywords:
; 
;
;    sublevels : in, optional, type=integer, default=3
;
;       The max number of levels of subdirectories to search.
;
;    fold : out, optional, type=boolean, default=true
;
;       Whether to fold cases when searching.
; 
;    count : out, optional, type=integer
;
;       The number of found subdirectories.
;
;
; :History:
; 
;    2016-05-05 : MGL. Inital version.
;
;
;
;-
function red_find_instrumentdirs, topdir, instrument_prefix, data_prefix $
                                  , sublevels = sublevels $
                                  , count = count $
                                  , fold = fold

  if n_elements(sublevels) eq 0 then sublevels = 3
  if n_elements(fold) eq 0 then fold = 1

  levels = ''
  for i = 0, sublevels-1 do begin
    print, 'level:', i
    instrumentdirs = file_search(topdir+'/'+data_prefix+'*/'+levels+instrument_prefix+'*', count = count, fold = fold)
    if count gt 0 then return, instrumentdirs
    levels += '*/'
  endfor
  
  count = 0
  return, ''
  
end
