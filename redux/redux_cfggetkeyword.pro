; docformat = 'rst'

;+
; Get a keyword from a redux cfg file.
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
;    The keyword value.
; 
; :Params:
; 
;     cfginfo : in, type="string, strarr, or struct"
; 
;        If strarr: The cfg file contents. If string: The
;        cfg file name. If struct: a struct of the type returned by
;        redux_cfg2struct. 
; 
;     name : in, type=string
;     
;        The name of the keyword. If ending in *, return multiple
;        keywords as a struct, where the * represents a number.
; 
; 
; :Keywords:
;
;     cfghash : out, optional, type=struct
;
;        The input cfginfo as an orderedhash. 
;     
;     count, out, optional, type=integer
; 
;        The number of values returned.
; 
; 
; :History:
; 
;    2020-02-17 : MGL. First version.
; 
;    2020-02-18 : MGL. Represent the cfg contents in the form of an
;                 orderedhash rather than a struct.
; 
;-
function redux_cfggetkeyword, cfginfo, name $
                              , cfghash = cfgshash $
                              , count = count
  
  if isa(cfginfo,'HASH')  then begin
    cfghash = cfginfo
  endif else begin
    cfghash = redux_cfg2hash(cfginfo)
  endelse

  if strmatch(name, '*.*') then begin
    ;; Compund keyword, use recursion
    splt = strsplit(name, '.', /extract)
    name_part1 = splt[0]
    name_part2 = strjoin(splt[1:*], '.')

    cfghash_part1 = redux_cfggetkeyword(cfghash, name_part1, count = count)
    
    case count of
      0 : return, []
      1 : return, redux_cfggetkeyword(cfghash_part1, name_part2, count = count)           
      else : begin
        ;; Multiple first parts, return array or list depending
        ;; on second part.    
        keys = cfghash_part1.keys()
        Nkeys = n_elements(keys)
        no_values_count = 0
        repeat begin
          value = redux_cfggetkeyword(cfghash_part1[keys[0]], name_part2, count = count_part2)
          if count_part2 eq 0 then no_values_count += 1
        endrep until count_part2 ne 0
        if isa(value, 'HASH') then begin
          if no_values_count gt 0 then begin
            values = list([])
            for i = 1, no_values_count-1 do values.add, []
            values.add, value
          endif else values = list(value)
          for i = no_values_count+1, count-1 do begin
            value = redux_cfggetkeyword(cfghash_part1[keys[i]], name_part2, count = count_part2)
            if count_part2 eq 0 then values.add, [] else values.add, value
          endfor                ; i
        endif else begin
          if no_values_count gt 0 then values = replicate('', no_values_count)
          red_append, values, value  
          for i = no_values_count+1, count-1 do begin
            value = redux_cfggetkeyword(cfghash_part1[keys[i]], name_part2, count = count_part2)
            if count_part2 eq 0 then red_append, values, '' else red_append, values, value  
          endfor                ; i
        endelse
        return, values
      end
    endcase
  endif


  keys = cfghash.keys()
  indx = where(strmatch(keys.ToArray(), name, /fold), count)
  case count of
    
    0 : begin
      return, []
    end

    1 : begin
      return, cfghash[keys[indx[0]]]
    end

    else : begin
      return, cfghash[keys[indx]]
    end
    
  endcase
  
end

dir = '/scratch/mats/mfbd_tail/CRISP/momfbd_remap_nopd/09:30:20/6563/cfg/'

value = redux_cfggetkeyword(dir+'momfbd_reduc_6563_00018.cfg', 'OBJECT*.CHANNEL1', count = count)
help, count, value

stop
value = redux_cfggetkeyword(dir+'momfbd_reduc_6563_00018.cfg', 'OBJECT*.CHANNEL1.FILENAME_TEMPLATE', count = count)
help, count, value

stop

value = redux_cfggetkeyword(dir+'momfbd_reduc_6563_00018.cfg', 'OBJECT*.WAVELENGTH', count = count)
help, count, value


stop

value = redux_cfggetkeyword(dir+'momfbd_reduc_6563_00018.cfg', 'OBJECT2.WAVELENGTH', count = count)
help, count, value
;
value = redux_cfggetkeyword(dir+'momfbd_reduc_6563_00018.cfg', 'OBJECT2.CHANNEL1.FILENAME_TEMPLATE', count = count)
help, count, value

value = redux_cfggetkeyword(dir+'momfbd_reduc_6563_00018.cfg', 'OBJECT2.CHANNEL1', count = count)
help, count, value

stop

value = redux_cfggetkeyword(dir+'momfbd_reduc_6563_00018.cfg', 'MAX_LOCAL_SHIFT', count = count)
help, count, value

stop

value = redux_cfggetkeyword(dir+'momfbd_reduc_6563_00018.cfg', 'OBJECT2', count = count)
help, count, value

value = redux_cfggetkeyword(dir+'momfbd_reduc_6563_00018.cfg', 'OBJECT*', count = count)
help, count, value


help, count, values

end

