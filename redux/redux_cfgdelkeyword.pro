; docformat = 'rst'

;+
; Delete a keyword from a redux cfg. Will quietly pretend to remove
; keywords that do not exist.
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
; :Params:
; 
;     cfginfo : in, type="string, strarr, or struct"
; 
;        If strarr or hash: The cfg file contents, to be updated. If
;        string: The cfg file name, to be rewritten.
; 
;     name : in, type=string
;     
;        The name of the keyword. If ending in *, return multiple
;        keywords as a struct, where the * represents a number.
; 
; :Keywords:
;
;     hcfg : out, optional, type=struct
;
;        The output cfg as an orderedhash for when cfginfo is a file
;        name.
;     
; 
; :History:
; 
;    2020-02-21 : MGL. First version.
; 
;-
pro redux_cfgdelkeyword, cfginfo, name $
                         , hcfg = hcfg 
  
  if isa(cfginfo,'HASH')  then begin
    hcfg = cfginfo
  endif else begin
    hcfg = redux_cfg2hash(cfginfo)
  endelse
  
  if ~strmatch(name, '*.*') then begin

    ;; Not compound, just remove the keyword (if it exists)
    if hcfg.HasKey(name) then hcfg.remove, name
    
  endif else begin
    
    ;; Compund keyword, use recursion
    
    splt = strsplit(name, '.', /extract)
    name_part1 = splt[0]
    name_part2 = strjoin(splt[1:*], '.')

    ;; The sub-hash to be updated
    hcfg_part1 = hcfg[name_part1]

    ;; Recursive call to update the subhash 
    redux_cfgdelkeyword, hcfg_part1, name_part2
    
  endelse 


  ;; Convert to input format or write to file as needed
  
  case 1 of

    isa(cfginfo, /string, /array) : cfginfo = redux_hash2cfg(hcfg)

    isa(cfginfo, /string, /scalar) : tmp = redux_hash2cfg(hcfg, filename = cfginfo)

    else : cfginfo = hcfg
    
  endcase
  
end

fname = '/scratch/mats/mfbd_tail/CRISP/momfbd_remap_nopd/09:30:20/6563/cfg/momfbd_reduc_6563_00018.cfg'

hcfg = redux_cfg2hash(fname)


mls = redux_cfggetkeyword(hcfg, 'MAX_LOCAL_SHIFT', count = count)
print, 'Before : MAX_LOCAL_SHIFT = ', mls
redux_cfgdelkeyword, hcfg, 'MAX_LOCAL_SHIFT'
mls = redux_cfggetkeyword(hcfg, 'MAX_LOCAL_SHIFT', count = count)
print, 'After  : MAX_LOCAL_SHIFT = ', mls

stop

redux_cfgdelkeyword, hcfg, 'OBJECT33.CHANNEL1.DARK_NUM'

print, redux_cfggetkeyword(hcfg, 'OBJECT33.CHANNEL1.WAVELENGTH')

stop

end

