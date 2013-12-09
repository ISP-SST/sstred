; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :author:
; 
; 
; 
; 
; :returns:
;    path to an .so library
; 
; :Params:
; 
;   
; 
; :Keywords:
; 
; 
; :history:
;
;   2013-12-09  initial version

FUNCTION red_libfile, libname

      ;;; first look in DLM_PATH
    libfile = file_search(strsplit(!DLM_PATH, ':', /extr), libname)
    IF libfile EQ '' THEN BEGIN
      ;;; try IDL_PATH, too
        libfile = file_search(strsplit(!PATH, ':', /extr), libname)
        IF libfile EQ '' THEN $
          message, 'Could not locate library file '+libname+'; Exiting'
    ENDIF
    libfile = libfile(0)
    return, libfile

END
