; docformat = 'rst'

;+
; Generate documentation web pages on github.
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, Institute for Solar Physics
;       with help from google's AI.
; 
; :History:
; 
;   2026-08-11 : MGL. First version.
; 
;-
PRO red_generate_docs
  ; 1. Add IDLdoc to the GDL path
  !PATH = './idldoc/src/:' + !PATH
  
  ; 2. Explicitly create the output directory
  FILE_MKDIR, 'out-docs'
  
  ; 3. Run IDLdoc (Removed the invalid RECURSIVE keyword)
  ;    We set ROOT to the current directory '.'
  idldoc, ROOT='.', $
          OUTPUT='out-docs/', $
          TITLE='RED Project Documentation'
          
  ; 4. Fallback: If IDLdoc created an html/ underfolder, move everything up
  IF FILE_TEST('out-docs/html/index.html') THEN BEGIN
    FILE_COPY, 'out-docs/html/*', 'out-docs/', /OVERWRITE
  ENDIF
  
  ; 5. Fallback: Ensure a standard index.html exists in out-docs/
  IF FILE_TEST('out-docs/index.pro.html') AND (NOT FILE_TEST('out-docs/index.html')) THEN BEGIN
    FILE_COPY, 'out-docs/index.pro.html', 'out-docs/index.html'
  ENDIF
          
  EXIT
END
