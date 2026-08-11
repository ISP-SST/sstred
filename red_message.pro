; docformat = 'rst'

;+
; Print a message to the terminal window, prepended with the name of
; the calling subprogram.
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
;   msg : in, type=strarr
; 
;     The message to be printed.
; 
; 
; :Keywords:
; 
;   
;   
;   
; 
; 
; :History:
; 
;   2025-02-19 : MGL. First version.
; 
;-
pro red_message, msg

  ;; Get the name of this subprogram and in particular of where it's
  ;; called from.
  inam = red_subprogram(/low, calling = inam1)              

  ;; Then print the message
  case n_elements(msg) of
    0 : return
    1 : red_strflow, inam1 + ' : ' + msg
    else : begin
      ;; If a strarr, we don't want to prepend the subprogram name to
      ;; each line.
      msg1 = msg
      msg1[0] = inam1 + ' : ' + msg1[0]
      red_strflow, msg1
    end
  endcase

end

