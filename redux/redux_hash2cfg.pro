; docformat = 'rst'

;+
; Convert a hash representation of cfg information to a string array
; and optionally store as a cfg file.
; 
; :Categories:
;
;    MFBD methods.
; 
; 
; :Author:
; 
;    Mats Löfdahl, Institute for Solar Physics
; 
; 
; :Returns:
; 
; 
; :Params:
; 
;    hcfg : in, type=orderedhash
; 
;       The cfg information to be converted.
; 
; 
; :Keywords:
; 
;    cfg : out, optional, type=strarr
;   
;       The converted output.
;
;    filename : in, optional, type=string
;
;       Where to store the output.
; 
; :History:
; 
;    2020-02-19 : MGL. First version.
; 
;-
pro redux_hash2cfg, hcfg, cfg = cfg, filename = filename
  
  keys = hcfg.keys()
  Nkeys = n_elements(keys)
  undefine, cfg
  
  for ikey = 0, Nkeys-1 do begin
    
;     print, keys[ikey]
    
    case 1 of
      
      strmatch(keys[ikey], 'OBJECT*') : begin
        ;; Recursion for object
        redux_hash2cfg, hcfg[keys[ikey]], cfg = ocfg
        red_append, cfg, 'object{'
        red_append, cfg, '  '+ocfg
        red_append, cfg, '}'
      end
      
      strmatch(keys[ikey], 'CHANNEL*') : begin
        ;; Recursion for channel
        redux_hash2cfg, hcfg[keys[ikey]], cfg = ccfg
        red_append, cfg, 'channel{'
        red_append, cfg, '  '+ccfg
        red_append, cfg, '}'
      end
        
      else : begin
        if isa(hcfg[keys[ikey]], /boolean) then begin
          ;; Treat booleans separately
          if hcfg[keys[ikey]] then begin
            ;; Just the keyword without a value means TRUE
            red_append, cfg, keys[ikey]
          endif else begin
            ;; If we'd ever need to specify a FALSE, that's a zero?
            red_append, cfg, keys[ikey]+'=0'
          endelse
        endif else begin
          red_append, cfg, keys[ikey]+'='+strtrim(hcfg[keys[ikey]], 2)
        endelse
      end
      
    endcase

  endfor                        ; ikey
  
  if n_elements(filename) gt 0 then begin
    openw, lun, filename, /get_lun
    printf, lun, cfg, format = '(a0)'
    free_lun, lun
  endif

end


dir = '/scratch/mats/mfbd_tail/CRISP/momfbd_remap_nopd/09:30:20/6563/cfg/'
hcfg = redux_cfg2hash(dir+'momfbd_reduc_6563_00018.cfg')

redux_hash2cfg, hcfg, filename = 'tmp.cfg', cfg = cfg


end
