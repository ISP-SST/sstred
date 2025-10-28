; docformat = 'rst'

;+
; Link (or optionally copy) summed calibrations to workdir (if they
; exist).
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
;   copy : in, optional, type=boolean
; 
;     Copy instead of making soft links.
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
;   2025-09-22 : MGL. New keyword link.
;   
;   2025-10-01 : MGL. New keyword copy, remove keyword link.
; 
;-
pro red_setupworkdir_copy, sourcedir, subdir, workdir, copy = copy

  if n_elements(sourcedir) eq 0 then return
  if file_test(sourcedir + '/' + subdir) eq 0 then return
  
  targetdir = workdir+'/'+subdir+'/'

  ;; Check subdirectories
  sss = file_search(sourcedir+'/'+subdir+'/','*.fits', count = Nsearch)
  tmpdirs = red_uniquify(file_dirname(sss))
  if subdir eq 'flats' then begin
    ;; Avoid copying output from fitgains
    tmpdirs = tmpdirs[where(~strmatch(tmpdirs,'*spectral_*'))]
  endif
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
  endfor                        ; ilist
  if n_elements(tmpdirs) gt 1 then begin
    print
    ;; If there are subdirectories with as well as without timestamps,
    ;; then default to selecting the ones with timestamps. (The ones
    ;; without are probably automatically generated when making
    ;; quicklooks.)
    default = where(strmatch(sel_list,'*??:??:??*'),Nwhere)
    if Nwhere eq n_elements(sel_list) or Nwhere eq 0 then default='*' else default = rdx_ints2str(default)
    tmp = red_select_subset(sel_list $
                            , qstring = 'What '+subdir+' directories to import from '+sourcedir+'?' $
                            , indx = indx $
                            , default = default $
                           )
    searchdir = tmpdirs[indx]
  endif else searchdir = tmpdirs[0]

  file_mkdir, targetdir
  
  for idir = 0, n_elements(searchdir)-1 do begin

    case subdir of

      'polcal_sums' :  begin
        ;; Copy/link files in subdirectories for polcal sums
        subdirs = file_search(searchdir[idir]+'/*', count = Nsubdirs)
        for isubdir = 0, Nsubdirs-1 do begin
          files = file_search(subdirs[isubdir]+'/cam*', count = Nfiles)
          if size(files,/n_dim) gt 0 then begin
            targetsubdir = targetdir + file_basename(red_uniquify(file_dirname(files)))+'/'
            file_mkdir, targetsubdir
            if ~keyword_set(copy) then begin
              ;; print, 'Link subdirectories '+strjoin(file_basename(subdirs), ',')+' from '+searchdir[idir]
              print, 'Link '+strtrim(Nfiles, 2)+' files from '+red_uniquify(file_dirname(files)) + '...' 
              file_link, files, targetsubdir, /allow_same
            endif else begin
              ;; print, 'Copy subdirectories '+strjoin(file_basename(subdirs), ',')+' from '+searchdir[idir]
              print, 'Copy '+strtrim(Nfiles, 2)+' files from '+red_uniquify(file_dirname(files)) + '...' 
              file_copy, files, targetsubdir, /overwrite
            endelse
          endif
        endfor                  ; isubdir
      end

      else : begin
        ;; Copy/link files for other types of summed calibration data
        files = file_search(searchdir[idir]+'/cam*', count = Nfiles)
        if size(files,/n_dim) gt 0 then begin
          if subdir eq 'flats' then begin
            ;; Avoid copying output from fitgains
            files = files[where(~strmatch(files, '*cavityfree*'), Nfiles)]
          endif
          if Nfiles gt 0 then begin
            if ~keyword_set(copy) then begin
              print, 'Link '+strtrim(Nfiles, 2)+' files from '+red_uniquify(file_dirname(files)) + '...' 
              file_link, files, targetdir, /allow_same
            endif else begin
              print, 'Copy '+strtrim(Nfiles, 2)+' files from '+red_uniquify(file_dirname(files)) + '...' 
              file_copy, files, targetdir, /overwrite
            endelse
          endif
        endif
      end
      
    endcase

;    
;
;
;    
;    if subdir eq 'polcal_sums' then begin
;
;      ;; Copy whole subdirectories for polcal sums
;      subdirs = file_search(searchdir[idir]+'/*')
;      if keyword_set(link) then begin
;        print, 'Link subdirectories '+strjoin(file_basename(subdirs), ',')+' from '+searchdir[idir]
;        file_link, subdirs, targetdir, /overwrite
;      endif else begin
;        print, 'Copy subdirectories '+strjoin(file_basename(subdirs), ',')+' from '+searchdir[idir]
;        file_copy, subdirs, targetdir, /overwrite, /recursive
;      endelse
;      
;    endif else begin
;      ;; Search for summed files and corresponding "_discarded.txt"
;      ;; files.
;      files = file_search(searchdir[idir]+'/cam*', count = Nfiles)
;      
;      if Nfiles gt 0 then begin
;        print, 'Copy '+strtrim(Nfiles, 2)+' files from '+red_uniquify(file_dirname(files)) + '...' 
;        file_copy, files, targetdir, /overwrite
;      endif
;      
;    endelse

  endfor                        ; idir
  

end
