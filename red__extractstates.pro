; docformat = 'rst'

;+
; Extract states information from an array of strings (typically file
; names). 
;
; Substitutes original crisp__extractstates.pro and
; chromis__extractstates.pro that were copied to
; crisp__extractstates_nondb.pro and chromis__extractstates_nondb.pro,
; respectively. This procedure is a wrapper around db and nondb
; versions.
; 
; :Categories:
;
;    SSTRED
; 
; 
; :Author:
; 
;     Oleksii Andriienko, ISP
; 
; 
; :Returns:
; 
; 
; :Params:
; 
;    strings : in, type=strarr
;   
;        A list of strings from which to extract the states
;        information. 
;   
;    states : out, optional, type=array(struct)
;
;        An array of structs, containing (partially filled) state
;        information. 
; 
; 
; :Keywords:
; 
;     datasets : in, optional, type = strarr
;
;        List of datasets timestamps to be used instead of list of
;        filenames. Can be used only if sst_db is installed.
; 
;     force : in, optional, type=boolean
;
;        Do not use cached states.
; 
;     nondb : in, optional, type=boolean
; 
; 
;     polcal : in, optional, type=boolean
; 
;        Set this to add polcal-specific items in the states, qw and
;        lp.
; 
;     strip_settings : in, optional, type=boolean
;
;        Exclude exposure/gain information from the fullstate entries.
;
; 
; :History:
; 
;   2019-07-23 : OA. Created.
; 
;   2022-07-29 : MGL. Change from a CRISP:: method to a RED:: method
;                to prepare for CRISP camera upgrade.
;
;-
pro red::extractstates, strings, states $
                        , cam = cam $
                        , datasets = datasets $ 
                        , force = force $
                        , nondb = nondb $
                        , polcal = polcal $
                        , strip_settings = strip_settings
  
  if keyword_set(datasets) then begin ; if we use datasets then we should use the database
    if n_elements(datasets) eq 0 then return
    self->extractstates_db, strings, states, datasets = datasets
    return
  endif
  Nstrings = n_elements(strings)
  
  if Nstrings eq 0 then return

  ;; Check for raw data directories in 'strings'.
  raw_data_dirs = [self -> raw_data_dirs(),'*/data/*']
  
  is_raw = 0B
  for i=0,n_elements(raw_data_dirs)-1 do begin
        bb = strmatch(strings,raw_data_dirs[i])
        ss = where(bb eq 1)
        if n_elements(ss) gt 1 then begin
      is_raw = 1B
      break  ; we should not have raw and processed data files in one call
    endif else if n_elements(ss) eq 1 and ss ne -1 then begin
      is_raw = 1B
      break 
    endif
  endfor
  if is_raw and self.db_present and ~keyword_set(nondb) then begin ; we can use sst_db only with raw data
    self->extractstates_db, strings, states, cam = cam
  endif else begin
    self->extractstates_nondb, strings, states $
                               , force = force $
                               , strip_settings = strip_settings $
                               , polcal = polcal
  endelse

  return

end


