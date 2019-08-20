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
;    instrument : in, optional, type=string, default='Get from dir'
; 
;       Specify the instrument.
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
;-
function red_raw_search, dir $
                         , count = count $
                         , fpi_states = fpi_states $
                         , instrument = instrument $
                         , prefilters = prefilters $
                         , scannos = scannos_in $
                         , tunings = tunings

  ;; Name of this subprogram
  inam = red_subprogram(/low, calling = inam1)

  if n_elements(instrument) eq 0 then begin
    ;; Find out from dir
    instrument = strupcase((strsplit(file_basename(dir), '-', /extract))[0])
  endif

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
    else : prefilters = strtrim(prefilters, 2)
  endcase
  Npref = n_elements(prefilters)
  
  ;; Tuning
  case n_elements(tunings) of
    0 : tunings = '[0-9][0-9][0-9][0-9]_[+-][0-9]*'
    else : tunings = strtrim(tunings, 2)
  endcase
  Ntuning= n_elements(tunings)

  ;; Fpi_state
  case n_elements(fpi_states) of
    0 : case strupcase(instrument) of
      'CHROMIS': fpi_states = 'wheel0000[0-9]_hrz[0-9][0-9][0-9][0-9][0-9]'
      'CRISP' : fpi_states = prefilters+'.'+tunings+'.lc?'
      else: stop
    endcase
    else : fpi_states = strtrim(fpi_states, 2)
  endcase
  Nstates = n_elements(fpi_states)
  
  
  ;; Construct search strings
  case strupcase(instrument) of
    
    'CRISP' : begin
      
      Nstrings = Nscans * Nstates
      searchstrings = strarr(Nstrings)
      istring = 0
      for iscan = 0, Nscans-1 do begin
        for istate = 0, Nstates-1 do begin
          searchstrings[istring] = 'cam*.' + scannos[iscan] + '.*.' $
                                   + fpi_states[istate] $
                                   + '.im.[0-9][0-9][0-9][0-9][0-9][0-9][0-9]'
          istring++
        endfor
      endfor 

    end

    'CHROMIS' : begin

      ;; Typical file name: sst_camXXX_00013_0009375_wheel00002_hrz33621.fits
      Nstrings = Nscans * Nstates
      searchstrings = strarr(Nstrings)
      istring = 0
      for iscan = 0, Nscans-1 do begin
        for istate = 0, Nstates-1 do begin
          searchstrings[istring] = 'sst_cam*_' + scannos[iscan] $
                                   + '_[0-9][0-9][0-9][0-9][0-9][0-9][0-9]_' $
                                   + fpi_states[istate] $
                                   + '.fits'
          istring++
        endfor
      endfor 

    end

    else : begin

      print, inam+' : This instrument is not implemented (yet): ', instrument
      
    end
    
  endcase 

;  files = file_search(dir+'/'+searchstrings, count = count)
  files = red_file_search(searchstrings, dir, count = count)
;  stop

  return, files

end

dir1 = '/data/2019/2019.04/2019.04.14/CHROMIS-data/08:46:09/Chromis-W/'
files1 = red_raw_search(dir1, scannos = '4,6-7,10' $
                        , fpi_states = ['wheel00006_*', 'wheel00002_hrz33578'], count = cnt1)


dir2 = '/data/2019/2019.04/2019.04.14/CHROMIS-data/08:46:09/Chromis-N/'
files2 = red_raw_search(dir2, scannos = '4,6-7,10', count = cnt2)

stop

dir3 = '/data/2019/2019.04/2019.04.14/Science/08:09:35/Crisp-W/'
files3 = red_raw_search(dir3, scannos = '4,6-7,10', pref = '8542', tunings = ['8542_-575', '8542_+700'], count = cnt3)
dir4 = '/data/2019/2019.04/2019.04.14/Science/08:09:35/Crisp-T/'
files4 = red_raw_search(dir4, scannos = '4,6-7,10', count = cnt4)

end

