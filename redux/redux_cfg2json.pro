; docformat = 'rst'

;+
; Convert cfg file contents to a json string.
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
;    The cfg file contents as a json string.
; 
; :Params:
; 
;    cfg_or_cfgfile : in, type="string or strarr"
; 
;      A cfg file name or its contents as a string array.
; 
; :History:
; 
;     2020-02-17 : MGL. First version.
; 
;     2020-02-18 : MGL. Make all keys upper case.
; 
;-
function redux_cfg2json, cfg_or_cfgfile
  
  case n_elements(cfg_or_cfgfile) of
    0: stop
    1 : spawn, 'cat '+cfg_or_cfgfile, cfg
    else: cfg = cfg_or_cfgfile
  endcase
  
  onumber = 0
  
  Nlines = n_elements(cfg)

  epos = strpos(cfg, '=')


  for iline = 0, Nlines-1 do begin
    
    cfgline = cfg[iline]
    
    if epos[iline] eq -1 then begin
      case 1 of
        strmatch(cfgline, '*object{*')  : begin
          ;; Begin object
          cnumber = 0
          cfgline = '"OBJECT'+strtrim(onumber, 2)+'" : {'
          onumber += 1
        end
        strmatch(cfgline, '*channel{*') : begin
          ;; Begin channel
          cfgline = '"CHANNEL'+strtrim(cnumber, 2)+'" : {'
          cnumber += 1
        end
        strmatch(cfgline, '*}*') : begin
          ;; End of object or channel
          cfgline = '},'
          ;; Remove trailing comma on previous line, if any,
          if strmid(cfg[iline-1], 0, /reverse) eq ',' then begin
            cfg[iline-1] = strmid(cfg[iline-1], 0, strlen(cfg[iline-1])-1)
          endif
        end
        else : begin
          ;; Boolean TRUE keyword
          cfgline = '"'+strtrim(cfgline, 2)+'": true,'
        end
      endcase
    endif else begin
      ;; Keyword-value pair on this line 
      keyword = strtrim(strmid(cfgline, 0, epos[iline]), 2)
      value   = strtrim(strmid(cfgline, epos[iline]+1), 2)
      cfgline = '"'+keyword+'": "'+value+'",'
    endelse
    
    cfg[iline] = cfgline
    
  endfor                        ; iline
  
  ;; Remove trailing comma on previous line, if any,
  if strmid(cfg[iline-1], 0, /reverse) eq ',' then begin
    cfg[iline-1] = strmid(cfg[iline-1], 0, strlen(cfg[iline-1])-1)
  endif
  
  return, strjoin(['{', cfg, '}'])
  
end

dir = '/scratch/mats/mfbd_tail/CRISP/momfbd_remap_nopd/09:30:20/6563/cfg/'
jstruct = redux_cfg2json(dir+'momfbd_reduc_6563_00018.cfg')

help, jstruct



end
