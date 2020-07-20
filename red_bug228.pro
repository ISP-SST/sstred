; Check momfbd config files for bug 228 (see
; https://dubshen.astro.su.se/redmine/issues/228). 
;
; This is a fix to a problem that can affect users who want to make
; fitscubes in workdirs where prepmomfbd was run prior to pulling the
; mid July update of SSTRED.
;
; In that update, I changed some metadata FITS header keywords so they
; would conform with the SOLARNET recommendations. I thought that
; would just propagate through all the processing steps to the final
; fitscubes but I had forgotten that one of those keywords is actually
; accessed by make_nb_cube, and if it (with its new value) is not
; found, some information in the wideband cube is not found.
;
; The symptom is that make_nb_cube stops with an error message similar to: 
;
; % Compiled module: RED::HEADERINFO_COPYSTEP.
; % Variable is undefined: STEPNUMS.
; % Execution halted at: RED::HEADERINFO_COPYSTEP 89 /home/user/idl/bin/sstred/red__headerinfo_copystep.pro
; % CHROMIS::MAKE_NB_CUBE 357 /home/user/idl/bin/sstred/chromis__make_nb_cube.pro
;
; The line numbers can vary due to later edits (and are definitely
; different for the crisp and chromis versions) but the problem
; appears when make_nb_cube calls headerinfo_copystep.
;
; To fix this, change directory to the CRISP or CHROMIS workdir. Start
; IDL and do
;
; IDL> red_bug228
;
; This program will loop through all momfbd*/.../results/*.fitsheader
; files, find any instances of "PRSTEP1 = 'MOMFBD image restoration'"
; and change them to "PRSTEP1 = 'MOMFBD '".
;
; After this has completed, please make new WB files before making any
; new NB files.
;
; Do *not* run fitscube_finalize on the wideband cube!
; 
; Mats Löfdahl 2020-07-20.
;
pro red_bug228

  files = file_search('momfbd*/*/*/cfg/results/*.fitsheader', count = Nfiles)

  print, 'Found '+strtrim(Nfiles, 2)+' fitsheader files.'
  
  for ifile = 0, Nfiles-1 do begin

    red_progressbar, ifile, Nfiles, files[ifile]
    
    h = headfits(files[ifile])

    value = strtrim(fxpar(h, 'PRSTEP1', count = Nkey), 2)

    if Nkey eq 0 then continue
    if value eq 'MOMFBD' then continue

    fxaddpar, h, 'PRSTEP1', 'MOMFBD'
    modfits, files[ifile], 0, h
    
  end
  
end
