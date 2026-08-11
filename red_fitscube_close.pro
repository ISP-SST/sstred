; docformat = 'rst'

;+
; Close a fitscube file that was opened with red_fitscube_open and
; optionally write a modified header.
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
;   fileassoc : in, type=assoc
;
;       The assoc variable of the opened file.
; 
;   fitscube_info : in, out, optional, type=struct
;
;       Info about the opened fitscube file, including the header, the
;       cube dimensions, and the logical unit. Undefined upon return.
; 
; :Keywords:
;
;    newheader : in, optional, type=strarr
;   
;       A (modified) header that should be written to the file.
; 
; 
; :History:
; 
;    2019-04-03 : MGL. First version.
; 
;    2019-09-13 : MGL. Use red_fitscube_newheader.
; 
;-
pro red_fitscube_close, fileassoc, fitscube_info $
                        , newheader = newheader

  inam = red_subprogram(/low, calling = inam1)

  lun = (size(fileassoc,/struc)).file_lun
 

  ;; Get the file name in case it's needed.
  fs = fstat(lun)
  filename = fs.name

  ;; Close the file
  free_lun, lun
  undefine, fitscube_info
  
  if n_elements(newheader) gt 0 then begin

;   ;; TODO? Check if the header length has changed enough that the
;   ;; data past has to be moved. If so, if might be faster to make a
;   ;; new copy of the file, with the new header, and then delete the
;   ;; original file. 
;   
;    print, inam + ' : Writing a modified header to a fitscube file. This may be very time consuming if the file is large and the data part has to be moved on disk.'
;
;    tic
;    modfits, filename, 0, newheader
;    toc

    red_fitscube_newheader, filename, newheader
    
  endif
  

end
