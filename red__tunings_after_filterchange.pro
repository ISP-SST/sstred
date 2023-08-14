; docformat = 'rst'

;+
; Find tunings right after a prefilter change.
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
;    A strarr with the tunings.
; 
; :Params:
; 
;    dir : in, type=string
; 
;       The path to a directory with narrowband data in fits files.
; 
; 
; :Keywords:
; 
;   count : out, optional
;   
;      The number of found tunings.
; 
; 
; :History:
; 
;   2023-05-29 : MGL. First version.
; 
;-
function red::tunings_after_filterchange, dir, count = count

  files = file_search(dir+'/*_00000_*.fits', count = Nfiles)
  
  ;; Are there any filter changes?
  prefs = red_fitsgetkeyword_multifile(files, 'FILTER1')
  upref = prefs[uniq(prefs,sort(prefs))]
  Nprefs = n_elements(upref)

  if Nprefs le 1 then begin
    ;; No prefilter changes, so return nothing
    count = 0
    return, 0
  endif

  ;; The number of prefilter changes is >= the number of prefilters,
  ;; since we have to change back to the first filter at the beginning
  ;; of the scan, and we may have changed between the filters multiple
  ;; times during the scan.

  ;; So we sort the files in temporal order and then detect filter
  ;; changes in the sorted array. The extractstates method
  ;; doesn't return the time of collecting the files so we'll read the
  ;; DATE-BEG keyword directory from the file headers.
  
  date_beg = red_fitsgetkeyword_multifile(files, 'DATE-BEG')
  sindx = sort(date_beg)

  files    = files[sindx]
  date_beg = date_beg[sindx]
  prefs    = prefs[sindx]
  
  tindx = [0, where(red_differential(long(prefs)) ne 0)]
  count = n_elements(tindx)

  ;; Old CHROMIS files have the tuning in scrambled form in the file
  ;; headers so we'll extract the states of the relevant files and
  ;; return the tuning from them.
  
  files = files[tindx]
  self -> extractstates, files, states
  
  return, states.tuning

end

if 1 then begin
  cd, '/scratch/mats/test-crisp2/2022-11-05/CRISP/'
  a = crisp2red(/dev, /no)
  dir = 'data/08:33:52/Crisp-R/'
endif else begin
  cd, '/scratch/mats/2016.09.19/CHROMIS-jan19'
  a = chromisred(/dev, /no)
  dir = 'data/10:42:01/Chromis-N/'
endelse

tuns = a -> tunings_after_filterchange(dir, count = count)
if count gt 0 then hprint, tuns

end
