; docformat = 'rst'

;+
; Convert a cfg file to an ordered hash.
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
;    The cfg file contents as an ordered hash.
; 
; :Params:
; 
;    cfg_or_cfgfile : in, type="string or strarr"
; 
;      A cfg file name or its contents as a string array.
; 
; :History:
; 
;     2020-02-18 : MGL. First version.
; 
;-
function redux_cfg2hash, cfg_or_cfgfile

  return, json_parse(redux_cfg2json(cfg_or_cfgfile)) 
     
end

dir = '/scratch/mats/mfbd_tail/CRISP/momfbd_remap_nopd/09:30:20/6563/cfg/'
hcfg = redux_cfg2hash(dir+'momfbd_reduc_6563_00018.cfg')

help, hcfg


end
