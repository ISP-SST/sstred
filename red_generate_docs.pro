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
  
  ; 3. Get the absolute path to the current directory
  CD, CURRENT=current_dir
  PRINT, '--- IDLdoc Root Directory: ', current_dir
  
  ; 4. Run IDLdoc using the absolute path for ROOT
  idldoc, ROOT=current_dir, $
          OUTPUT='out-docs/', $
          TITLE='RED Project Documentation'
          
  ; 5. Fallback: If IDLdoc placed files in an html/ subfolder, move them up
  IF FILE_TEST('out-docs/html/index.html') THEN BEGIN
    FILE_COPY, 'out-docs/html/*', 'out-docs/', /OVERWRITE
  ENDIF
          
  EXIT
END
