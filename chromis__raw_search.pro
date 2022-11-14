; docformat = 'rst'

;+
; Search for raw files.
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
;    A string array with file names. 
; 
; :Params:
; 
;    dir : in, type=string
;   
;      The directory in which to search.
; 
; 
; :Keywords:
; 
;    count : out, optional, type=integer
;   
;       The number of found files.
; 
;    fpi_states : in, optional, type=strarr
; 
;       Limit the search to files with these fpi_states as part of
;       their file names.
; 
;    prefilters : in, optional, type=strarr
; 
;       Limit the search to files with these prefilters as part of
;       their file names.
; 
;    scannos : in, optional, type="array of strings or integers"
; 
;       Limit the search to files with these scannumbers as part of
;       their file names. Can be the numbers in integer or string
;       form, or a comma and dash separated list.
; 
;    tunings : in, optional, type=strarr
; 
;       Limit the search to files with these tunings as part of
;       their file names.
; 
; :History:
; 
;  2019-07-05 : MGL. First version.
; 
;  2021-08-23 : MGL. Adapt searchstrings to WB flats file names.
; 
;  2021-08-23 : MGL. Make it into a method.
; 
;-
function chromis::raw_search, dir $
                              , count = count $
                              , fpi_states = fpi_states $
                              , prefilters = pref $
                              , scannos = scannos_in $
                              , tunings = tunings

  ;; Name of this subprogram
  inam = red_subprogram(/low, calling = inam1)

  iswb = 'W' eq strupcase((strsplit(file_basename(dir), '-', /extract))[1])
  isflats = strmatch(dir,'*[Ff]lat*')

  ;; Massage the scannos
  case n_elements(scannos_in) of

    0 : scannos = '[0-9][0-9][0-9][0-9][0-9]' ; Proper wildcard

    1 : if size(scannos_in, /tname) eq 'STRING' then begin
      scannos = string(rdx_str2ints(scannos_in), format = '(I05)') ; Possibly a comma and dash separated string
    endif else begin
      scannos = string(round(scannos_in), format = '(I05)') ; A number
    endelse

    else : scannos = string(round(scannos_in), format = '(I05)') ; An array, assume of numbers.
    
  endcase
  Nscans = n_elements(scannos)
  
  ;; Prefilter
  case n_elements(prefilters) of
    0 : prefilters = '[0-9][0-9][0-9][0-9]'
    else : prefilters = strtrim(pref, 2)
  endcase
  Npref = n_elements(prefilters)
  
  ;; Tuning
  case n_elements(tunings) of
    0 : tunings = '[0-9][0-9][0-9][0-9]_[+-][0-9]*'
    else : tunings = strtrim(tunings, 2)
  endcase
  Ntuning = n_elements(tunings)

  ;; Fpi_state
  case n_elements(fpi_states) of
    0 : fpi_states = 'wheel0000[0-9]_hrz[0-9][0-9][0-9][0-9][0-9]' ;prefilters+'_'+tunings+'_lc?'
    else : fpi_states = strtrim(fpi_states, 2)
  endcase
  Nstates = n_elements(fpi_states)

  ;; Construct search strings
  ;; Typical file name: sst_camXXX_00013_0009375_wheel00002_hrz33621.fits
  Nstrings = Nscans * Nstates
  searchstrings = strarr(Nstrings)
  istring = 0
  for iscan = 0, Nscans-1 do begin
    for istate = 0, Nstates-1 do begin
      if self.isodate ge '2022-11-03' then begin
        if iswb and isflats then begin
          ;; Typical name: sst_camXXXI_00000_0001250_5896.fits
          searchstrings[istring] = 'sst_cam*_' + scannos[iscan] $
                                   + '_[0-9][0-9][0-9][0-9][0-9][0-9][0-9]_' $
                                   + prefilters $
                                   + '.fits'
        endif else begin
          ;; Typical name: sst_camXXXI_00002_0002076_5896_5896_+0092_lc4.fits
          searchstrings[istring] = 'sst_cam*_' + scannos[iscan] $
                                   + '_[0-9][0-9][0-9][0-9][0-9][0-9][0-9]_' $
                                   + prefilters $
                                   + '_*.fits'
        endelse
      endif else begin 
        if iswb and isflats then begin
          searchstrings[istring] = 'sst_cam*_' + scannos[iscan] $
                                   + '_[0-9][0-9][0-9][0-9][0-9][0-9][0-9]_' $
                                   + 'wheel0000[0-9]' $ ; Will have to change when new file naming scheme is implemented
                                   + '.fits'
        endif else begin
          searchstrings[istring] = 'sst_cam*_' + scannos[iscan] $
                                   + '_[0-9][0-9][0-9][0-9][0-9][0-9][0-9]_' $
                                   + fpi_states[istate] $
                                   + '.fits'
        endelse
      endelse
      istring++
    endfor
  endfor 

  stop
  
  files = red_file_search(searchstrings, dir, count = count)
  
  return, files

end
