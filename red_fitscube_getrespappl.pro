; docformat = 'rst'

;+
; Read RESPAPPL (APPLied RESPonse function) variable keyword.
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
;      The name of the file in which to look for the RESPAPPL variable
;      keyword.  
; 
; 
; :Keywords:
; 
;   axis_numbers : out, optional, type=array
;   
;      The fitscube pixel axis numbers that the response function was
;      applied to.
; 
;   count : out, optional, type=integer
; 
;      The number of instances of the RESPAPPL keyword that was found.
; 
;   scans_dimension : out, optional, type=boolean
; 
;      Whether the scans dimension is represented in the returned
;      response function.
;
;   tuning_dimension : out, optional, type=boolean
; 
;      Whether the tuning dimension is represented in the returned
;      response function.
; 
; 
; :History:
; 
;   2021-02-26 : MGL. First version.
; 
;-
function red_fitscube_getrespappl, filename $
                                   , axis_numbers = axis_numbers $
                                   , count = count $
                                   , scans_dimension = scans_dimension $
                                   , tuning_dimension = tuning_dimension

  tmp = red_fitsgetkeyword(filename, 'RESPAPPL' $
                           , count = count $
                           , variable_values = variable_values )
  
  ;; Probably always response correction only for tuning and scans
  ;; dimensions. But define this just in case it changes in the
  ;; future.
  axis_numbers = [3, 5]
  scans_dimension = !true
  tuning_dimension = !true
  
  ;; Find the dimensions of the data cube
  h = headfits(filename)
  Ntun   = fxpar(h, 'NAXIS3')
  Nscans = fxpar(h, 'NAXIS5')

  if count eq 0 then begin
    ;; No RESPAPPL keyword in the file but return array with 1 as a
    ;; convenience. 
    return, replicate(1d, Ntun, Nscans)
  endif

  
  ;; What dimensions were stored?
  dims0 = size(variable_values.values, /dim)
  Ntun0 = dims0[3-2]
  if n_elements(dims0) gt 5-2 then Nscans0 = dims0[5-2] else Nscans0 = 1

  if Ntun eq Ntun0 and Nscans eq Nscans0 then begin
    ;; Dimensions match, just return the values array
    respappl = reform(variable_values.values, Ntun, Nscans)
  endif else begin
    ;; Reform and rebin to the cube [Ntun,Nscans] dimensions
    respappl = rebin( reform(variable_values.values, Ntun0, Nscans0) $
                      , Ntun, Nscans, /sample)
  endelse

  return, respappl
  
end
