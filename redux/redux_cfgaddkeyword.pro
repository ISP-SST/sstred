; docformat = 'rst'

;+
; Add a keyword to a redux cfg. 
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
;     value : in, type=varies
;
;        The value to be assigned to the keyword.
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
;    2020-02-19 : MGL. First version.
; 
;    2020-03-07 : MGL. Create hash if needed.
; 
;-
pro redux_cfgaddkeyword, cfginfo, name, value $
                         , hcfg = hcfg 
  
  if n_elements(cfginfo) eq 0 then begin
    ;; Undefined cfg, make one
    hcfg = orderedhash([], []) 
  endif else begin
    if isa(cfginfo,'HASH')  then begin
      hcfg = cfginfo
    endif else begin
      hcfg = redux_cfg2hash(cfginfo)
    endelse
  endelse
  
  if ~strmatch(name, '*.*') then begin

    ;; Not compound, just insert the value

    case 1 of
      n_elements(value) eq 0 : hcfg[name] = !true
      isa(value, /array) : hcfg[name] = strjoin(strtrim(value, 2), ',')
      else : hcfg[name] = strtrim(value, 2)
    endcase
    
  endif else begin
    
    ;; Compund keyword, use recursion
    
    splt = strsplit(name, '.', /extract)
    name_part1 = splt[0]
    name_part2 = strjoin(splt[1:*], '.')

    ;; Existing sub-hash to be updated
    if hcfg.haskey(name_part1) then hcfg_part1 = hcfg[name_part1]

    ;; Recursive call to update/create the subhash 
    redux_cfgaddkeyword, hcfg_part1, name_part2, value 

    ;; Assign the updated subhash as the value to the first part of
    ;; the name 
    hcfg[name_part1] = hcfg_part1
    
  endelse 


  ;; Convert to input format or write to file as needed
  
  case 1 of

    isa(cfginfo, /string, /array) : cfginfo = redux_hash2cfg(hcfg)

    isa(cfginfo, /string, /scalar) : tmp = redux_hash2cfg(hcfg, filename = cfginfo)

    else : cfginfo = hcfg
    
  endcase
  
end

;; Test creating with dot notation
redux_cfgaddkeyword, hhh, 'TOP.MID.BOTTOM', '45'



;; Testing
fname = '/scratch/mats/mfbd_tail/CRISP/momfbd_remap_nopd/09:30:20/6563/cfg/momfbd_reduc_6563_00018.cfg'

hcfg = redux_cfg2hash(fname)


mls = redux_cfggetkeyword(hcfg, 'MAX_LOCAL_SHIFT', count = count)
print, 'Before : MAX_LOCAL_SHIFT = ', mls
redux_cfgaddkeyword, hcfg, 'MAX_LOCAL_SHIFT', mls*2
mls = redux_cfggetkeyword(hcfg, 'MAX_LOCAL_SHIFT', count = count)
print, 'After  : MAX_LOCAL_SHIFT = ', mls


redux_cfgaddkeyword, hcfg, 'OBJECT33.CHANNEL1.TJO', '1,2,34'
redux_cfgaddkeyword, hcfg, 'OBJECT33.CHANNEL1.TJO', 3.14

print, redux_cfggetkeyword(hcfg, 'OBJECT33.CHANNEL1.TJO')

stop

end

