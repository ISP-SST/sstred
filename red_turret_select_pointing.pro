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
;    An array with pointing coordinates.
; 
; :Params:
; 
;   isodate : in, type=string
; 
;      The date in ISO format. (Needed only if we need to convert
;      Stonyhurst coordinates from old turret logfiles.)
; 
;   pointingtype : in, type=string
; 
;      The pointing type tag to select, one of "Stonyhurst", "Disk",
;      "Flat", or "Pointing". If "Pointing", information from the
;      former three types are combined and returned as cartesian disk
;      coordinates. 
; 
;   turretdata : in, type=structarr
; 
;      A struct array with info from the turret logfile. The turret
;      program stores pointing information of different types
;      depending on what it's doing, these types are tagged in
;      different ways.
; 
; 
; :Keywords:
; 
;   time : in, out, optional, type=dblarr
;   
;      This is an array of time coordinates in s after midnight. If
;      given as input, the returned pointing coordinates will be
;      interpolated to those points in time. If not, it is simply
;      returned as the time coordinate in turretdata.
; 
; 
; :History:
; 
;    2017-04-23 : MGL. First version.
; 
;    2021-12-16 : MGL. Implemented combined pointingtype "Pointing".
; 
; 
;-
function red_turret_select_pointing, turretdata, pointingtype, isodate, time = time

  ;; Return positions in array with length matched to the time array.
  if n_elements(time) eq 0 then begin
    time = turretdata.time
    interpolate_time = 0
  endif else begin
    interpolate_time = 1
  endelse
  Ntimes = n_elements(time)

  pointing = replicate(!Values.F_NaN, 2, Ntimes) ; Return NaNs where no data is available  
  
  ;; Find data of the specified type
  if pointingtype ne 'Pointing' then begin
    
    plength = strlen(pointingtype)
    indx = where( (strmid(turretdata.pointingtype1, 0, plength) eq pointingtype) $
                  and (strmid(turretdata.pointingtype2, 0, plength) eq pointingtype) $
                  , count, complement=cindx, ncomplement=ccount)
    
    
    if count eq 0 then begin
      ;; No coordinates of the specified type available. Return NaNs to
      ;; signal no data available.
      return, pointing
    end
    
    if interpolate_time then begin
      
      ;; Interpolate but avoid gaps in data
      pointing[0, *] = red_interpol_nogaps(turretdata[indx].pointing1, turretdata[indx].time, time)
      pointing[1, *] = red_interpol_nogaps(turretdata[indx].pointing2, turretdata[indx].time, time)

    endif else begin
      
      ;; No interpolation, just use the available data
      pointing[0, indx] = turretdata[indx].pointing1
      pointing[1, indx] = turretdata[indx].pointing2
      
    endelse
    
    return, pointing
    
  endif

  ;; Return the combination of "Stonyhurst", "Disk", and "Flat"

  pointing_all = replicate(!Values.F_NaN, 2, n_elements(turretdata)) 

  ;; The Cartesian disk coordinates are exactly what we want
  indx_disk = where( (strmid(turretdata.pointingtype1, 0, 4) eq 'Disk') $
                     and (strmid(turretdata.pointingtype2, 0, 4) eq 'Disk') $
                     , count_disk, complement=cindx_disk, ncomplement=ccount_disk)
  if count_disk gt 0 then begin
    pointing_all[0, indx_disk] = turretdata[indx_disk].pointing1
    pointing_all[1, indx_disk] = turretdata[indx_disk].pointing2
    red_append, indx, indx_disk
  endif
  
  ;; The Flats mode data are also Cartesian, only with a bit worse
  ;; accuracy due to them being the wanted coordinates rather than the
  ;; actual coordinates.
  indx_flat = where( (strmid(turretdata.pointingtype1, 0, 4) eq 'Flat') $
                     and (strmid(turretdata.pointingtype2, 0, 4) eq 'Flat') $
                     , count_flat, complement=cindx_flat, ncomplement=ccount_flat)
  if count_flat gt 0 then begin
    pointing_all[0, indx_flat] = turretdata[indx_flat].pointing1
    pointing_all[1, indx_flat] = turretdata[indx_flat].pointing2
    red_append, indx, indx_flat
  end
  
  ;; The Stonyhurst coordinates have to be converted to Cartesian disk
  ;; coordinates. 
  indx_stony = where( (strmid(turretdata.pointingtype1, 0, 5) eq 'Stony') $
                      and (strmid(turretdata.pointingtype2, 0, 5) eq 'Stony') $
                      , count_stony, complement=cindx_stony, ncomplement=ccount_stony)
  if count_stony gt 0 then begin
    helio = replicate(!Values.F_NaN, 2, n_elements(indx_stony))
    helio[0, *] = turretdata[indx_stony].pointing2 ; 2nd value is Longitude
    helio[1, *] = turretdata[indx_stony].pointing1 ; 1st value is Latitude
    pointing_all[*, indx_stony] = red_convert_stonyhurst(helio, isodate, turretdata.time)
    red_append, indx, indx_stony
  end
    
  if n_elements(indx) eq 0 then stop ; Should return something

  ;; Sort indices, potentially concatenation of indx_disk, indx_flat,
  ;; and indx_stony, in no order.
  sindx = sort(turretdata[indx].time)
  indx = indx[sindx]

  
  ;; If no interpolation then return all the data points
  if ~interpolate_time then return, pointing_all[indx]
  
  ;; Interpolate but avoid gaps in data
  pointing[0, *] = red_interpol_nogaps(pointing_all[0, indx], turretdata[indx].time, time)
  pointing[1, *] = red_interpol_nogaps(pointing_all[1, indx], turretdata[indx].time, time)

  return, pointing
  
end
