; docformat = 'rst'

;+
; Copy summed calibrations to workdir (if they exist).
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
; 
; :Params:
; 
;   sourcedir : in, type=string
; 
;     Path to but not including the subdirectory to be copied.
;
;   subdir : in, type=string
; 
;      The subdirectory to be copied.
; 
;   workdir : in, type=string
; 
;      Where to copy the files.
; 
; :History:
; 
;   2025-04-25 : MGL. First version.
; 
;-
pro red_setupworkdir_copy, sourcedir, subdir, workdir 
  
  if n_elements(sourcedir) gt 0 then begin
  
    targetdir = workdir+'/'+subdir+'/'

    ;; Check subdirectories
    sss = file_search(sourcedir+'/'+subdir+'/','*.fits', count = Nsearch)
    tmpdirs = red_uniquify(file_dirname(sss))
    if subdir eq 'polcal_sums' then tmpdirs = red_uniquify(file_dirname(tmpdirs))
    ;; Make a selection list
    sel_list = strarr(n_elements(tmpdirs))
    for ilist = 0, n_elements(tmpdirs)-1 do begin
      sel_list[ilist] = red_strreplace(tmpdirs[ilist], sourcedir, '') 
      if subdir eq 'polcal_sums' then begin
        fnames = file_search(tmpdirs[ilist]+'/*/cam*fits', count = Nfiles)
      endif else begin
        fnames = file_search(tmpdirs[ilist]+'/cam*fits', count = Nfiles)
      endelse
      red_extractstates, fnames, /basename, cam = cams, pref = prefilters, settings = settings
      prefilters = red_uniquify(prefilters)
      cams = red_uniquify(cams)
      settings = red_uniquify(settings)
      if n_elements(cams) gt 0 && cams[0] ne '' then sel_list[ilist] += ' - ' + strjoin(cams, ',')
      if n_elements(settings) gt 0 && settings[0] ne '' then sel_list[ilist] += ' - ' + strjoin(settings, ',')
      if n_elements(prefilters) gt 0 && prefilters[0] ne '' then sel_list[ilist] += ' - ' + strjoin(prefilters, ',')
    endfor                      ; ilist
    if n_elements(tmpdirs) gt 1 then begin
      print
      tmp = red_select_subset(sel_list $
                              , qstring = 'What '+subdir+' directories to import from '+sourcedir+'?' $
                              , indx = indx $
                             )
      searchdir = tmpdirs[indx]
    endif else searchdir = tmpdirs[0]

    file_mkdir, targetdir

    for idir = 0, n_elements(searchdir)-1 do begin

      if subdir eq 'polcal_sums' then begin

        ;; Copy whole subdirectories for polcal sums
        subdirs = file_search(searchdir[idir]+'/*')
        print, 'Copy subdirectories '+strjoin(file_basename(subdirs), ',')+' from '+searchdir[idir]
        file_copy, subdirs, targetdir, /overwrite, /recursive
        
      endif else begin
        ;; Search for summed files and corresponding "_discarded.txt"
        ;; files.
        files = file_search(searchdir[idir]+'/cam*', count = Nfiles)
        
        if Nfiles gt 0 then begin
          print, 'Copy '+strtrim(Nfiles, 2)+' files from '+red_uniquify(file_dirname(files)) + '...' 
          file_copy, files, targetdir, /overwrite
        endif
        
      endelse

    endfor                      ; idir
    
  endif

end
