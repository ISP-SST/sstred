; docformat = 'rst'

;+
; Read the same header keyword from multiple files.
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
;   A string array with keyword values.
; 
; :Params:
; 
;   filenames : in, type=strarr
; 
;     Paths to the files to read the header keywords from.
; 
;   name : in, type=string
; 
;     The keyword name.
; 
; :Keywords:
; 
;   comments : out, optional, type=strarr
;   
;     The keyword comments.
; 
;   counts : out, optional, type=lonarr
;   
;     The number of keywords found for each file.
; 
; :History:
; 
;    2023-05-22 : MGL. First version.
; 
;-
function red_fitsgetkeyword_multifile, filenames, name $
                                       , comments = comments $
                                       , counts = counts $
                                       , _ref_extra = extra

  Nfiles = n_elements(filenames)

  comments = strarr(Nfiles)
  counts = lonarr(Nfiles)
  
  for ifile = 0, Nfiles-1 do begin

    hdr = red_readhead(filenames[ifile])
    value = red_fitsgetkeyword(hdr, name $
                               , comment = comment $
                               , count = count $
                               , _strict_extra = extra)
    
    red_append, values, value   ; Use red_append because we don't know the type
    counts[ifile] = count
    if count gt 0 then comments[ifile] = comment
    
  endfor                        ; ifile

  return, values
  
end
