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
  ; 1. Add IDLdoc to the GDL path so the command can be found
  !PATH = './idldoc/src/:' + !PATH
  
  ; 2. Run IDLdoc over the entire project
  ;    ROOT='.' starts the search from the repository root.
  ;    RECURSIVE=1 ensures all subdirectories are scanned.
  ;    OUTPUT='out-docs/' defines where the HTML pages will be saved.
  ;    EXCLUDE prevents IDLdoc from documenting itself or the output folder.
  idldoc, ROOT='.', $
          OUTPUT='out-docs/', $
          TITLE='RED Project Documentation', $
          RECURSIVE=1, $
          EXCLUDE=['idldoc', 'out-docs']
          
  ; 3. Exit GDL when the documentation generation is complete
  EXIT
END
