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
  
    if subdir eq 'polcal_sums' then begin
      searchdir = file_search(sourcedir+'/'+subdir+'/*', count = Ndirs)
      if Ndirs eq 0 then return
      targetdir = workdir+'/'+subdir+'/' + file_basename(searchdir)+'/'
    endif else begin
      searchdir = sourcedir+'/'+subdir+'/'
      targetdir = workdir+'/'+subdir+'/'
    endelse
    
    for idir = 0, n_elements(searchdir)-1 do begin
      
      ;; Search for summed files and corresponding "_discarded.txt"
      ;; files.
      files = file_search(searchdir[idir]+'/cam*', count = Nfiles)

      file_mkdir, targetdir[idir]

      if Nfiles gt 0 then begin
        str = 'Copy '+strtrim(Nfiles, 2)+' files from '+red_uniquify(file_dirname(files)) + '...'
        print, string(13B), str, FORMAT = '(A,A,$)'
        file_copy, files, targetdir[idir], /overwrite
        print, string(13B), str + 'Done!', FORMAT = '(A,A,$)'
        print
      endif
      
    endfor                      ; idir
    
  endif

end
