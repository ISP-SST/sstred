; docformat = 'rst'

;+
; Sort the keywords in a redux cfg. 
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
;    2020-03-09 : MGL. First version.
; 
;-
pro redux_cfgsortkeywords, cfginfo, hcfg = hcfg 

  if n_elements(cfginfo) eq 0 then begin
    return                      ; Undefined cfg, nothing to sort
  endif else begin
    if isa(cfginfo,'HASH')  then begin
      hcfg = cfginfo            ; Input is already a hash
    endif else begin
      hcfg = redux_cfg2hash(cfginfo) ; Input is a string or strarr
    endelse
  endelse

  keys = hcfg.keys()

  ;; Objects or channels?
  okeys = keys.filter('strmatch', 'OBJECT[0-9]*')
  ckeys = keys.filter('strmatch', 'CHANNEL[0-9]*')
  Nobjects  = n_elements(okeys)
  Nchannels = n_elements(ckeys)

  if Nobjects gt 0 and Nchannels gt 0 then begin
    print, 'OBJECT and CHANNEL keywords should not be on the same level'
  endif

  case 1 of

    Nchannels gt 0 : begin

      ;; Recursive calls to sort within the channels
      for ikey = 0, Nchannels-1 do begin
        ccfg = hcfg[ckeys[ikey]]
        redux_cfgsortkeywords, ccfg
        hcfg[ckeys[ikey]] = ccfg
      endfor
      
      ;; Sort the channel keys in numerical order
      ckeys = ckeys.sort(compare_function = $
                         Lambda(a,b:(long(strmid(a, 7))).Compare(long(strmid(b, 7)))) $
                        )
      
      ;; Find non-channel keywords
      indx = keys.where(ckeys.toarray(), complement=cindx, Ncomplement = Ncomplement)

      if Ncomplement eq 0 then begin
        ;; The channel keys are the only keys
        keys = ckeys
      endif else begin
        ;; Sort the non-channel keys and concatenate.
        keys = (keys[cindx]).sort() + ckeys
      endelse

    end
    
    Nobjects gt 0 : begin

      ;; Recursive calls to sort within the objects
      for ikey = 0, Nobjects-1 do begin
        ocfg = hcfg[okeys[ikey]]
        redux_cfgsortkeywords, ocfg
        hcfg[okeys[ikey]] = ocfg
      endfor

      ;; Sort the object keys in numerical order
      okeys = okeys.sort(compare_function = $
                         Lambda(a,b:(long(strmid(a, 6))).Compare(long(strmid(b, 6)))) $
                        )
      
      ;; Find non-object keywords
      indx = keys.where(okeys.toarray(), complement=cindx, Ncomplement = Ncomplement)

      
      if Ncomplement eq 0 then begin
        ;; The object keys are the only keys
        keys = okeys
      endif else begin
        ;; Sort the non-object keys and concatenate.
        keys = okeys + (keys[cindx]).sort()
      endelse

    end

    else : begin

      ;; Only ordinary keywords
      keys = keys.sort()

    end
    
  endcase

  ;; Reorder the hash
  hcfg = hcfg[keys]

  ;; Convert to input format or write to file as needed
  
  case 1 of

    isa(cfginfo, /string, /array) : cfginfo = redux_hash2cfg(hcfg)

    isa(cfginfo, /string, /scalar) : tmp = redux_hash2cfg(hcfg, filename = cfginfo)

    else : cfginfo = hcfg
    
  endcase
  
  
end
