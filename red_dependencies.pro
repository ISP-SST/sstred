;; Use this program to locate subprograms that are not in the crispred
;; directory, IDL distribution, astrolib, mpfit or the coyote library.
;;
;; IDL> .run red_dependencies
;;
;; Two output files are generated. 
;;
;;  * dependencies.out : if verbose>0: lists all dependencies, else it only lists the unresolved ones
;;
;;  * dependencies.grep : lists the crispred files where the found subprograms are called.
;;
;; Before making local copies of the found subprograms and changing
;; their names to the red_* namespace, make sure this isn't
;; already done so you only have to change how they are called.
;;
;; Finding everything is a recursive process. Run again after fixing
;; what was found.
;;
;; Mats Löfdahl 2013-07-24.
;;
;;   2015-08-11 : Turn into procedure in the red namespace to avoid potential name collisions.
;;                Automagically resolve valid directories (crispred, astrolib, coyote, mpfit).
;;                Ignore some unresolved, but unused, routines.

pro prnt, text, lun
    print, text
    printf, lun, text
end

pro rcheck, rname, lun, is_function=is_function, okpaths=okpaths, skip_routines=skip_routines, verbose=verbose

    for i=0, (size(skip_routines, /dim))[0]-1 do begin
        if strmatch(rname, skip_routines[i], /fold) then return
    endfor

    rinfo = routine_info(rname,/source,functions=is_function)

    isOK = 0
    if rinfo.path eq '' or FILE_DIRNAME(rinfo.path) eq '' then begin
        prnt, 'Routine not found: '+rname, lun
    endif else begin
        for i = 0, (size(okpaths, /dim))[0]-1 do begin
            isOK = isOK or strmatch(rinfo.path, okpaths[i], /fold)
        endfor
        if isOK then begin
            if keyword_set(verbose) then prnt, 'Located "'+rname+'"  --->  '+rinfo.path, lun
        endif else begin
            prnt, 'Routine "'+rname+'" found outside accepted scope: '+rinfo.path, lun
       endelse
    endelse
    
    if ~isOK then spawn, 'grep -i '+rname+' '+FILE_DIRNAME(ROUTINE_FILEPATH("prnt"))+PATH_SEP()+'*.pro >> dependencies.grep'
    
end

pro red_dependencies, verbose=verbose

    spawn, 'rm -f dependencies.grep'
    openw, lun, 'dependencies.out', /get_lun

    crispred_dir = FILE_DIRNAME(ROUTINE_FILEPATH("prnt"))   ;; get location of this file (i.e. crispred)
    resolve_routine, 'ZENPOS'                               ;; some random routine in the astro lib
    astrolib_dir = FILE_DIRNAME(FILE_DIRNAME(ROUTINE_FILEPATH("ZENPOS")))
    resolve_routine, 'cgAxis'                               ;; some random routine in the coyote lib
    coyote_dir = FILE_DIRNAME(ROUTINE_FILEPATH("cgAxis"))
    resolve_routine,'binomial',is_function=1                ;; some random function in the idl lib
    idl_dir = FILE_DIRNAME(ROUTINE_FILEPATH("binomial",/is_function))
    resolve_routine,'mpfit',is_function=1                   ;; some random function in the mpfit lib
    mpfit_dir = FILE_DIRNAME(ROUTINE_FILEPATH("mpfit",/is_function))

    okpaths = [crispred_dir, idl_dir, coyote_dir, astrolib_dir, mpfit_dir]+PATH_SEP()+'*'
    skip_routines=[ '$main$', 'boolean', 'retrun', 'trnlog', 'dellog', 'setlog']
    ;; $main$:   meta-routine that has no source path and will always fail our check.
    ;; boolean:  added in IDL 8.4 and will not be called on older versions
    ;; retrun:   bug in idl 7.0, lib/parse_url.pro has a type-o "retrun" should be "return" :-)
    ;; the last 3 are obsoleted but still referenced since they are used in the doc_library debugging

    if keyword_set(verbose) then begin
        prnt, okpaths, lun
        prnt, skip_routines, lun
    endif

    ;; compile all methods and their dependencies
    resolve_all, class='pol', skip_routines=skip_routines, /CONTINUE_ON_ERROR, /QUIET
    resolve_all, class='red', skip_routines=skip_routines, /CONTINUE_ON_ERROR, /QUIET
    resolve_all, class='chromis', skip_routines=skip_routines, /CONTINUE_ON_ERROR, /QUIET
    resolve_all, class='crisp', skip_routines=skip_routines, /CONTINUE_ON_ERROR, /QUIET
    
    scripts = file_basename(file_search( crispred_dir+PATH_SEP()+'*.pro' ), '.pro')     ; all scripts in crispred
    idx = where( strmatch(scripts, '*__*', /fold) EQ 1, compl=cidx )                    ; only the non class-methods
    scripts = scripts[cidx]
    for i=0, n_elements(scripts)-1 do begin
        resolve_all, resolve_either=scripts[i], skip_routines=skip_routines, /CONTINUE_ON_ERROR, /QUIET
    endfor

    subrnames = strlowcase(routine_info(functions=0))
    funcnames = strlowcase(routine_info(functions=1))

    for i=0, (size(subrnames, /dim))[0]-1 do begin
        rcheck, subrnames[i], lun, okpaths=okpaths, skip_routines=skip_routines, verbose=verbose
    endfor

    for i=0, (size(funcnames, /dim))[0]-1 do begin
        rcheck, funcnames[i], lun, okpaths=okpaths, skip_routines=skip_routines, verbose=verbose, /is_function
    endfor


    free_lun, lun


end
