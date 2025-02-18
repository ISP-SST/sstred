; docformat = 'rst'

;+
; Copy FITS info about a processing step from one FITS header or file
; to another.
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
; 
; :Params:
; 
;    filename_or_header : in, type=strarr
;
;      The FITS header in which to write the step info. If only one
;      element or a scalar string, interpreted as the name of a file
;      from which to read the header. Note: needs to be a filename if
;      any of the step keywords are are variable keywords.
; 
;    old_filename_or_header : in, type=strarr
;
;      The FITS header from which to read the step info. If only one
;      element or a scalar string, interpreted as the name of a file
;      from which to read the header. Note: needs to be a filename if
;      any of the step keywords are are variable keywords.
; 
; :Keywords:
; 
;    all : in, optional, type=boolean
;   
;       Copy info for all available steps.
; 
;    anchor  : in, optional, type=string
;   
;       Position the step info after the specified keyword.
; 
;    stepnum  : in, optional, type=integer
;   
;       Copy info for the step with this (these) number(s). 
; 
;    last  : in, optional, type=boolean
;   
;       Copy info for the last available step.
; 
;    prstep  : in, optional, type=string
;   
;       Copy info for steps matching this string.
; 
; 
; :History:
; 
;    2018-10-31 : MGL. First version.
; 
;    2020-07-01 : MGL. Reimplement.
; 
;-
pro red::headerinfo_copystep, filename_or_header, old_filename_or_header $
                              , all = all $
                              , anchor = anchor $
                              , last = last $
                              , prstep = prstep $
                              , stepnum = stepnum 

  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)

  case n_elements(filename_or_header) of
    
    0 : begin
      print, inam + ' : No header or filename provided.'
      stop
      ;; return
    end
    
    1 : begin
      ;; Is filename_or_header a file name?
      if file_test(filename_or_header[0]) then begin
        filename = filename_or_header[0]
        hdr = headfits(filename)
      endif else begin
        ;; Nope
        print, inam + ' : filename_or_header is not a file name.'
        help, filename_or_header[0]
        stop
      end
    end
    
    else : hdr = filename_or_header
    
  endcase
  
 
  case n_elements(old_filename_or_header) of
    
    0 : begin
      print, inam + ' : No old header or filename provided.'
      return
    end
    
    1 : begin
      ;; Is old_filename_or_header a file name?
      if file_test(old_filename_or_header[0]) then begin
        old_filename = old_filename_or_header[0]
        oldhdr = headfits(old_filename)
      endif else begin
        ;; Nope
        print, inam + ' : old_filename_or_header is not a file name.'
        help, old_filename_or_header[0]
        return
      end
    end
    
    else : oldhdr = old_filename_or_header
    
  endcase

  ;; Available steps
  prsteps_available = strtrim(fxpar(oldhdr,'PRSTEP*', count = Navailable), 2)
  
  ;; Any steps to copy?
  if Navailable eq 0 then return

  ;; Make step name changes backward compatible
  for istep = 0, Navailable-1 do begin
    if prsteps_available[istep] eq 'Demodulate' then $
       prsteps_available[istep] = 'DEMODULATION'
    if prsteps_available[istep] eq 'MOMFBD image restoration' then $
       prsteps_available[istep] = 'MOMFBD'
  endfor

  ;; Existing steps
  prsteps_existing = strtrim(fxpar(hdr,'PRSTEP*', count = Nexisting), 2)
  
  ;; List of defined PR* keywords.
  prkeys = red_headerinfo_prkeys(count = Nkeys) 
  
  ;; Combine steps indicated by the keywords
  if keyword_set(all) then red_append, stepnums, indgen(Navailable)+1
  if keyword_set(last) then red_append, stepnums, Navailable
  if n_elements(stepnum) gt 0 then red_append, stepnums, stepnum
  for istep = 0, n_elements(prstep)-1 do begin
    ii = where(prsteps_available eq prstep[istep], Nmatch)
    if Nmatch gt 0 then red_append, stepnums, ii+1
  endfor                        ; istep

  ;; Some step could be indicated by multiple keywords, make a
  ;; uniqueified list
  stepnums = stepnums(uniq(stepnums, sort(stepnums))) 
   
  ;; Find the next output stepnumber and the anchor.
  red_headerinfo_anchor, hdr, anchor, next_stepnumber = outstep

;  outstep = Nexisting
;  if Nexisting gt 0 then begin
;
;    ;; Look for last PR* keyword to use as anchor
;    pos = 0
;    ;;keywrd = strmid(oldhdr, 0, 8)
;    keywrd = strmid(hdr, 0, 8)
;    iend = n_elements(keywrd)
;    for ikey = 0, n_elements(prkeys)-1 do begin
;      for istep = 0, n_elements(stepnums)-1 do begin
;        foreach letter, ['', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J'] do begin
;          name = prkeys[ikey]+strtrim(stepnums[istep], 2)+letter
;          thispos = fxparpos(keywrd, iend, before = name)
;          
;          if thispos eq iend then break
;          if thispos gt pos then begin
;            pos = thispos
;            anchor = name
;          endif
;        end
;      endfor                    ; istep
;    endfor                      ; ikey
;    
;  endif else begin
;    ;; Given anchor overridden if there are steps already
;    if n_elements(anchor) eq 0 then anchor = 'OBS_HDU' 
;  endelse

  ;; Now do the copying
  for istep = 0, n_elements(stepnums)-1 do begin
    for ikey = 0, n_elements(prkeys)-1 do begin
      foreach letter, ['', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J'] do begin

        name = prkeys[ikey]+strtrim(stepnums[istep], 2)+letter

;        print, name

        ;; Read keyword
        undefine, variable_values
        if n_elements(old_filename) ne 0 then begin
          value = red_fitsgetkeyword(old_filename, name $
                                     , comment = comment, count = count $
                                     , variable_values = variable_values)
        endif else begin
          value = red_fitsgetkeyword(oldhdr, name $
                                     , comment = comment, count = count)
        endelse

        if count eq 0 then break

        ;; Write keyword
        outname = prkeys[ikey]+strtrim(outstep, 2)+letter
        print, outname
        if n_elements(variable_values) gt 0 then begin

          if n_elements(filename) eq 0 then stop
          ;;output[strtrim(keys[ikey], 2)] = {value:value, comment:comment, variable_values:variable_values}

          ;; fitscube_addvarkeyword reads and writes the header so we have
          ;; to write what we have changed so far...
          red_fitscube_newheader, filename, hdr

          axis_numbers = where(size(variable_values.values, /dim) gt 1) + 1
          
          self -> fitscube_addvarkeyword, filename, outname, reform(variable_values.values) $
                                          , anchor = anchor $
                                          , comment = comment $
                                          , keyword_value = value $
                                          , axis_numbers = axis_numbers

          ;; ...and then read it again.
          hdr = headfits(filename)

        endif else begin
          red_fitsaddkeyword, anchor = anchor, hdr, outname, value, comment
        endelse

      end                       ; letter 
    endfor                      ; ikey
    outstep++                   ; Step number in the new header
  endfor                        ; istep

  if n_elements(filename) gt 0 then begin
    ;; Write the new header to the file
    red_fitscube_newheader, filename, hdr
  endif
  
end


a = crispred(/dev)

hdr = headfits('/scratch/mats/2016.09.19/CRISP-aftersummer/cubes_wb/wb_6302_2016-09-19T09:30:20_scans=2-8_corrected_im.fits')

mkhdr, newhdr, fltarr(10, 10)

hgrep, hdr, 'PR'

print

a -> headerinfo_copystep, newhdr, hdr, stepnum = 1
hprint, newhdr

a -> headerinfo_copystep, newhdr, hdr, stepnum = 2
hprint, newhdr

end
