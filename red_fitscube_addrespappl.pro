; docformat = 'rst'

;+
; Add or uppdate RESPAPPL (APPLied RESPonse function) variable keyword.
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
;   filename : in, type=string
; 
;      The name of the file in which to add or update the RESPAPPL
;      variable keyword.
; 
;   respappl : in, type=array
; 
;      The applied response function for one of or both the tuning and
;      scans dimensions.
; 
; 
; :Keywords:
;
;   anchor : in, out, optional, type=string
;
;      Place the RESPAPPL header keyword after this keyword.
; 
;   axis_numbers : out, optional, type=array
;   
;      The fitscube pixel axis numbers that the response function was
;      applied to.
; 
;   scans_dimension : in, optional, type=boolean
; 
;      Whether the scans dimension is represented in the response
;      function.
;
;   tuning_dimension : in, optional, type=boolean
; 
;      Whether the tuning dimension is represented in the response
;      function.
; 
; :History:
; 
;   2021-02-26 : MGL. First version.
; 
;-
pro red_fitscube_addrespappl, filename, respappl $
                              , anchor = anchor  $
                              , axis_numbers = axis_numbers $
                              , scans_dimension = scans_dimension $
                              , tuning_dimension = tuning_dimension $
                              , update = update

  ;; Read existing RESPAPPL keyword, if any.
  old_respappl = red_fitscube_getrespappl(filename $
                                          , count = count $
                                          , scans_dimension = old_scans_dimension $
                                          , tuning_dimension = old_tuning_dimension)

  ;; Get dimensions from existing respappl
  old_dims = size(old_respappl, /dim)
  if n_elements(old_dims) eq 1 then begin
    old_respappl = reform(old_respappl, old_dims, 1)
    old_dims = [old_dims, 1]    
  endif
  Ntun = old_dims[0]
  Nscans = old_dims[1]

;  help, count, scans_dimension, tuning_dimension, old_scans_dimension, old_tuning_dimension
;  print, old_respappl
;  print, respappl
  
  ;; Reform to Ntun by Nscans if necessary and multiply with old RESPAPPL.
  if ~keyword_set(scans_dimension) then begin
    if n_elements(respappl) ne old_dims[0] then stop
    new_respappl = old_respappl * rebin( reform(respappl $
                                                , n_elements(respappl), 1) $
                                         , Ntun, Nscans, /sample)
  endif else if ~keyword_set(tuning_dimension) then begin
    if n_elements(respappl) ne old_dims[1] then stop
    new_respappl = old_respappl * rebin( reform(respappl $
                                                , 1, n_elements(respappl)) $
                                         , Ntun, Nscans, /sample)    
  endif else new_respappl = old_respappl * respappl
 
  ;; Write the new keyword.
  red_fitscube_addvarkeyword, filename, 'RESPAPPL', new_respappl $ 
                              , anchor = anchor $
                              , comment = 'Mean of applied response function' $
                              , keyword_method = 'mean' $
                              , axis_numbers = [3, 5] $
                              , update = update

end
