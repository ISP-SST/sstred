; docformat = 'rst'

;+
; 
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
; 
; 
; 
; 
; 
; :Keywords:
; 
;   
;   
;   
; 
; 
; :History:
; 
;    2017-04-23 : MGL. First version.
; 
; 
; 
; 
;-
function red_turret_select_pointing, turretdata, pointingtype, time = time

  ;; Return positions in array with length matched to the time array.
  if n_elements(time) eq 0 then begin
    time = turretdata.time
    interpolate_time = 0
  endif else begin
    interpolate_time = 1
  endelse
  Ntimes = n_elements(time)
  pointing = fltarr(2, Ntimes)

  ;; Find data of the specified type
  plength = strlen(pointingtype)
  indx = where( (strmid(turretdata.pointingtype1, 0, plength) eq pointingtype) $
                and (strmid(turretdata.pointingtype2, 0, plength) eq pointingtype) $
                , count, complement=cindx, ncomplement=ccount)

  if count eq 0 then begin
    ;; No coordinates of the specified type available. Return NaNs to
    ;; signal no data available.
    pointing[*, *] = !Values.F_NaN
    return, pointing
  endif
  
  if interpolate_time then begin
    
    ;; Interpolate but avoid gaps in data
    pointing[0, *] = red_interpol_nogaps(turretdata[indx].pointing1, turretdata[indx].time, time)
    pointing[1, *] = red_interpol_nogaps(turretdata[indx].pointing2, turretdata[indx].time, time)

  endif else begin
    
    ;; No interpolation, just use the available data
    pointing[0, indx] = turretdata[indx].pointing1
    pointing[1, indx] = turretdata[indx].pointing2

  endelse
  
  ;; Return NaNs where no data is available
  if ccount gt 0 then begin
    pointing[0, cindx] = !Values.F_NaN
    pointing[1, cindx] = !Values.F_NaN
  endif

  return, pointing
  
end
