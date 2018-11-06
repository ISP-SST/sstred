; docformat = 'rst'

;+
; Get information about the memory use on the local host.
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
;   An anonymous struct representing the contents of /proc/meminfo.
;   The memory amounts are given in bytes.
; 
; :History:
; 
;   2018-11-05 : MGL. First version.
; 
;-
function red_meminfo

  spawn, 'cat /proc/meminfo', meminfo

  Nlines = n_elements(meminfo)

  for iline = 0, Nlines-1 do begin

    docontinue = 0
    
    items = strsplit(meminfo[iline], /extract)
    
    case n_elements(items) of
      2: if items[1] eq '0' then begin ; If no unit, then value should be 0.
        value = 0d 
        tag = strsplit(items[0], ':', /extract)
      endif
      3: begin
        tag = strsplit(items[0], ':', /extract)
        case items[2] of
          'kB' : value = double(items[1])*1d3
          'MB' : value = double(items[1])*1d6
          'GB' : value = double(items[1])*1d9
          'TB' : value = double(items[1])*1d12
          else : stop
        endcase
      end
      else: docontinue = 1
    endcase

    if docontinue then continue ; Can't do this from within the case statement.
    
    tag = red_strreplace(tag, '(', '_')
    tag = red_strreplace(tag, ')', '')
    
    if n_elements(memstruct) eq 0 then begin
      memstruct = create_struct(tag, value)
    endif else begin
      memstruct = create_struct(memstruct, tag, value)
    endelse
    
  endfor                        ; iline

  return, memstruct
  
end
