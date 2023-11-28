; docformat = 'rst'

;+
; Combine info in the last PRSTEPn and associated keywords from
; multiple scans.  
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
;    hdr : out, type=strarr
; 
;       The header of cubefile after adding the combined step.
; 
;    cubefile : in, type=string
; 
;       A fitscube file.
; 
;    files : in, type=strarr
; 
;       An array of file names, from which to combine the (last) step
;       info and put in cubefile.
; 
;    stepnumber : in, out, type=integer
; 
;       The stepnumber to write. Incremented by one upon return.
; 
; 
; :Keywords:
; 
;    anchor : in, optional, type=string
;     
;       Anchor for writing FITS keywords.
; 
; 
; :History:
; 
;    2022-12-09 : MGL. First version.
; 
;-
pro red::headerinfo_combinestep, hdr, cubefile, files, stepnumber, anchor = anchor
  
  ;; Name of this subprogram
  inam = red_subprogram(/low, calling = inam1)

  outstp = strtrim(stepnumber, 2) ; The step to write to cubefile

  prkeys = red_headerinfo_prkeys(count = Nkeys)
  
  Nfiles = n_elements(files)
  
  ;; Available steps from files, we want the last (and hopefully only)
  ;; step. (The method has to be rewritten if we want to do this for
  ;; more than a single step from the files or if the step number
  ;; could vary in the different files.)
  hdr1 = red_readhead(files[0])
  prsteps_available = strtrim(fxpar(hdr1,'PRSTEP*', count = fstep), 2)
  instp = strtrim(stepnumber, 2) ; The step to read from files

  hdr = headfits(cubefile)
  ;; Assume that this step is not already in hdr or cubefile!
    

  values = strarr(Nfiles)
  comments = strarr(Nfiles)
  for ikey = 0, Nkeys-1 do begin

    inkey = prkeys[ikey]+instp
    outkey = prkeys[ikey]+outstp

    ;; Read key values from all files
    nokey = 0
    for ifile = 0, Nfiles-1 do begin

      file = files[ifile]
      if rdx_filetype(file) eq 'MOMFBD' then file = red_strreplace(file, '.momfbd', '.fitsheader')
      
      val = red_fitsgetkeyword(file, inkey, comment = comment, count = cnt)

      if cnt eq 0 then begin
        nokey = 1               ; This keyword does not exist in file
      endif else begin
        values[ifile] = val
        comments[ifile] = comment
      endelse
      
    endfor                      ; ifile

    if nokey then continue      ; Skip to next keyword
    
    if inkey eq 'PRSTEP'+instp then begin
      ;; Some old momfbd output could have an old version of the
      ;; PRSTEP1 keyword. Repair that here.
      indx = where(values eq 'MOMFBD image restoration', cnt)
      if cnt gt 0 then values[indx] = 'MOMFBD' 
    endif

    
    if stddev(total(float(byte(values)),1)) lt 1e-10 then begin
      ;; If key values are all the same, then add just the scalar
      ;; keyword.
      print, 'Add '+outkey+' as a regular keyword.'
      red_fitsaddkeyword, hdr, outkey, values[0], comments[0], anchor = anchor
    endif else begin

      ;; Just take the first value for this
      scalar_value = values[0]
      
      ;; We could think of a way to combine also the comments. Or pick
      ;; the first one?
      scalar_comment = 'Variable keyword'
      unique_comments = comments[uniq(comments, sort(comments))]
      scalar_comment = strjoin(unique_comments, ',')
      
      ;; Add variable keyword to the file.

      print, 'Add '+outkey+' as a VARIABLE keyword.'

      ;; fitscube_addvarkeyword reads and writes the header so we have
      ;; to write what we have changed so far...
      red_fitscube_newheader, cubefile, hdr

      self -> fitscube_addvarkeyword, cubefile, outkey, values $
                                      , anchor = anchor $
                                      , comment = scalar_comment $
                                      , keyword_value = scalar_value $
                                      , axis_numbers = [5]

      ;; ...and then read it again.
      hdr = headfits(cubefile)

    endelse
    
  endfor                        ; ikey
  
  red_fitscube_newheader, cubefile, hdr

  if fstep gt 0 then stepnumber++
  
end

