; docformat = 'rst'

;+
; Return date when something changed in the setup or data collection.
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
;    A date in ISO format YYYY-MM-DD or tags for such a date.
; 
; :Params:
; 
; 
; :Keywords:
; 
;    count : out, optional, type=string
;   
;      The number of returned dates or tags.
; 
;    date : in, out, optional, type=string
;   
;      The date for which to return tags.
; 
;    explanation : out, optional, type=string
;   
;      A longer explanation of the relvant event corresponding to the
;      date and tag.
; 
;    tag : in, out, optional, type=string
;   
;      The tag for which to return a date.
;   
; 
; 
; :History:
; 
;   2025-02-24 : MGL. First version.
; 
;-
function red_dates $
   , count = count $
   , date = date $
   , explanation = explanation $
   , tag = tag

  ;; Name of this subprogram
  inam = red_subprogram(/low, calling = inam1)  

  ;; Define dates, tags, and explanations

  ;; Note: You can add tags but existing tags should never be changed.
  ;; If you do it anyway, you need to make sure to change it also in
  ;; any subprogram that uses it. (Also potentially in other branches.
  ;; So don't do it. Really!)
  
  red_append, tags, 'AO KL'
  red_append, dates, '2013-01-01'
  red_append, explanations, 'AO uses KL modes from '+dates[-1]

  red_append, tags, 'AO 8x8'
  red_append, dates, '2013-10-28'
  red_append, explanations, 'AO log has 8x8 measurements from '+dates[-1]

  red_append, tags, 'CRISP Ximea'
  red_append, dates, '2022-08-19'
  red_append, explanations, 'Installed new CRISP Ximea cameras from '+dates[-1]
  
  red_append, tags, 'CHROMIS tuning metadata'
  red_append, dates, '2022-11-03'
  red_append, explanations, 'CHROMIS tuning metadata avaialble from '+dates[-1]
  
  red_append, tags, 'pinhole array with L'
  red_append, dates, '2023-01-01'
  red_append, explanations, 'Pinhole array with larger pinhole in L configuration installed after '+dates[-1]

  red_append, tags, 'polcal flats'
  red_append, dates, '2023-10-12'
  red_append, explanations, 'Collection of flats for the polcal wavelengths routinely started on '+dates[-1]

  red_append, tags, 'CHROMIS Ximea'
  red_append, dates, '2025-04-01'
  red_append, explanations, 'Installed new CHROMIS Ximea cameras from '+dates[-1]

  
  if n_elements(tag) gt 0 && n_elements(date) gt 0 then begin

    print, inam + ' : Please use only one of tag or date keywords'
    count = 0
    return, ''
    
  endif
  
  if n_elements(tag) gt 0 then begin

    ;; Return matching dates and explanations
    indx = where(tag eq tags, count)
    if count eq 1 then begin
      explanation = explanations[indx]
      return, dates[indx]
    endif else begin
      print, inam + ' : Unknown or ambiguous tag: '+tag
      count = 0
      return, ''
    endelse
  
  endif

  if n_elements(date) gt 0 then begin

    ;; Return matching tags and explanations
    indx = where(date eq dates, count)
    if count gt 0 then begin
      explanation = explanations[indx]
      return, tags[indx]
    endif else begin
      print, inam + ' : Unknown date: '+date
      count = 0
      return, ''
    endelse

  endif

  print, inam + ' : Please use one of tag or date keyword'
  count = 0
  return, ''  

end
