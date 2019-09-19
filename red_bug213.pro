; Check momfbd config files for bug 213 (see
; https://dubshen.astro.su.se/redmine/issues/213). 
;
; Change directory to a CRISP workdir. Start IDL and do
;
; IDL> red_bug213
;
; This program will examine one config file in all
; momfbd*/??:??:??/????/cfg/ directories and report in the text file
; bug213.out if it finds at least one case where the ALIGN_MAPs for a
; tuning state differ between lc states. (The bug should not have any
; dramatic effects for non-polarimetric data.)
;
; For such directories, look in the corresponding NB cube for
; alignment artifacts in the differential Stokes components. In
; particular for tuning states that there are pinhole array data for,
;
; The solution to the problem is to do a red_update to get a version
; of prepmomfbd where this is fixed (not available yet, still being
; tested). And then, for the data sets where you find that the
; artifacts matter, delete the momfbd*/??:??:??/????/cfg/results/
; directories (including any stokes*/ subdirectories therein), rerun
; prepmomfbd, rerun momfbd restoration, and rerun whatever steps you
; do after that like cube making.
;
;
; Below please find sample output for a config file where the problem
; was detected. Note that the first ALIGN_MAP line differs from the
; three following lines. This is the symptom of the bug.
;
; Examining momfbd_nopd/10:33:29/6173/cfg/momfbd_reduc_6173_00000.cfg
; Alignment maps differ for at least one tuning.
; Tuning 6173_-175 : 
; ALIGN_MAP=-1.00166,-0.00866148,1001.91,-0.00916016,0.998210,8.01989,-2.38452e-06,-1.07953e-06,1.00000
; ALIGN_MAP=-1.00173,-0.00864516,1001.97,-0.00912338,0.998334,8.07150,-2.39547e-06,-1.04110e-06,1.00000
; ALIGN_MAP=-1.00173,-0.00864516,1001.97,-0.00912338,0.998334,8.07150,-2.39547e-06,-1.04110e-06,1.00000
; ALIGN_MAP=-1.00173,-0.00864516,1001.97,-0.00912338,0.998334,8.07150,-2.39547e-06,-1.04110e-06,1.00000
; 
; Mats Löfdahl 2019-09-19.
;
pro red_bug213

  restore, 'calib/alignments.sav'

  ;; Select a detector
  detector = alignments[0].state2.detector

  ;; Find and loop over all momfbd cfg directories
  cfgdirs = file_search('momfbd*/??:??:??/????/cfg', count = Ndirs)

  if Ndirs eq 0 then return

  openw, lun, 'bug213.out', /get_lun
  
  for idir = 0, Ndirs-1 do begin

    ;; Check one of the cfg files
    cfgfiles = file_search(cfgdirs[idir]+'/*.cfg', count = Nfiles)

    if Nfiles eq 0 then continue

    printf, lun
    printf, lun, 'Examining '+cfgfiles[0]

    
    spawn, 'cat '+cfgfiles[0], cfg

    
    
    ;; Check if polarimetric
    if total(strmatch(cfg,'*lc2*')) eq 0 then continue ; non-polarimetric
    
    ;; Determine prefilter
    prefilter = (strsplit(cfgdirs[idir], '/', /extract))[2]
    
    ;; Find tunings with alignments
    indx = where(alignments.state2.prefilter eq prefilter, Ntunings)
    tunings = alignments[indx].state2.tuning
    
    for ituning = 0, Ntunings-1 do begin

      ;; Read maps for such tunings, all LC states
      indx = where(strmatch(cfg, '*FILENAME_TEMPLATE='+detector+'*'+tunings[ituning]+'_*'), Nmatch)
      if Nmatch eq 0 then continue
      
      for i = indx[0], indx[1]-1 do if strmatch(cfg[i], '*ALIGN_MAP*') then break
      indx += (i-indx[0])
      maps = cfg[indx]

      if round(total(maps[0] ne maps[1:*])) gt 0 then begin

        ;; If maps differ, then report in output file
        printf, lun, 'Alignment maps differ for at least one tuning.'
        printf, lun, 'Tuning '+tunings[ituning]+' : '
        printf, lun, maps
        break
        
      endif

      
    endfor                      ; ituning
    
  endfor                        ; idir

  free_lun, lun
  
end
