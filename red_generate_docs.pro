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
  
  ; 2. Explicitly create the output directory so GDL doesn't crash
  FILE_MKDIR, 'out-docs'
  
  ; 3. Run IDLdoc over the entire project
  idldoc, ROOT='.', $
          OUTPUT='out-docs/', $
          TITLE='RED Project Documentation', $
          RECURSIVE=1
          
  ; 4. Explicitly exit GDL
  EXIT
END
