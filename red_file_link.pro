; docformat = 'rst'

;+
; Like IDL's file_link but allow overwriting of existing links (but
; not regular files).
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
; :Params:
; 
;   SourcePath : in, type="string or strarr"
; 
;     The target files. See IDL's file_link command.
; 
;   DestPath
; 
;     Where to make the links. See IDL's file_link command.
; 
; 
; :Keywords:
; 
;   overwrite : in, optional, type=boolean
;   
;     Allow overwriting existing links.
; 
;   _ref_extra : in, optional
;
;     See IDL's file_link for other supported keywords.
;
; :History:
; 
;    2026-07-02 : MGL. First version.
; 
;-
pro red_file_link, SourcePath, DestPath $
                   , overwrite = overwrite $
                   , VERBOSE = VERBOSE, _REF_EXTRA=extra

  COMPILE_OPT idl2
  
  case !true of

    ~keyword_set(overwrite) : begin
      ;; If we don't want to overwrite, let IDL's file_link do its thing.
      if keyword_set(verbose) then red_message, 'No overwrite.'
    end
 
    isa(DestPath, /scalar) : begin

      if keyword_set(verbose) then red_message, 'DestPath is a scalar.'
      
      if file_test(DestPath) then begin
        
        ;; Scenario A: DestPath is a directory
        if file_test(DestPath, /dir) then begin
          if keyword_set(verbose) then $
             red_message, 'DestPath is a directory. Expanding SourcePath and checking for conflicting symlinks inside.'
          
          ;; Expandera SourcePath om det innehåller wildcards (t.ex. '*.dat')
          ;; Om SourcePath redan är en vanlig array eller filnamn, returnerar FILE_SEARCH det som det är.
          resolved_sources = file_search(SourcePath, count = Nsearch)
          
          ;; If SourceDest expands to something then check if there
          ;; are matching files in DestPath that need to be deleted.
          ;; Otherwise, let file_link deal with it.
          if Nsearch gt 0 then begin
            base_names = file_basename(resolved_sources)
            
            ;; Loop through the files
            for i = 0L, n_elements(base_names) - 1L do begin
              check_file = DestPath + '/' + base_names[i]
              if file_test(check_file, /symlink) and ~file_test(check_file, /dir) then begin
                if keyword_set(verbose) then $
                   red_message, 'Deleting existing symlink inside directory: ' + check_file
                file_delete, check_file, /quiet
              endif
            endfor
          endif
          
        ;; Scenario B: DestPath is not a directory
        endif else begin
          if file_test(DestPath, /symlink) then begin
            if keyword_set(verbose) then $
               red_message, 'DestPath is a symlink and not a directory, so delete it.'
            file_delete, DestPath, /quiet
          endif
        endelse
        
      endif
    end

    isa(SourcePath, /array) : begin
      ;; Both are arrays. They'd better be of the same length.
      if keyword_set(verbose) then $
         red_message, 'SourcePath and DestPath are both arrays.'
      if n_elements(SourcePath) ne n_elements(DestPath) then begin
        red_message, 'SourcePath and DestPath are different lengths. Do nothing.'
        return
      endif
      ;; Each element of SourcePath can potentially be a regular
      ;; expression. If it matches several files, the corresponding
      ;; element of DestPath has to be an existing directory. To avoid
      ;; dealing with too many complicated cases here, use recursive
      ;; calls to figure each element out.
      if keyword_set(verbose) then $
         red_message, 'Call recursively for each element.'
      for ipath = 0, n_elements(SourcePath)-1 do begin
        red_file_link, SourcePath[ipath], DestPath[ipath], /overwrite $
                   , VERBOSE = VERBOSE, _STRICT_EXTRA=extra    
      endfor                    ; ipath
      return
    end

    else : begin
      if keyword_set(verbose) then $
         red_message, 'SourcePath is a scalar and DestPath is an array. This is not allowed.'
      return
    end
    
  endcase

  if keyword_set(verbose) then $
         red_message, "Let IDL's file_link do it's thing."
  file_link, SourcePath, DestPath, VERBOSE = VERBOSE, _STRICT_EXTRA=extra


end


