; docformat = 'rst'

;+
; Download the FPI calibration values (linedef.py).
;
; Code from the red::download method.
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
;   isodate, in, type=string
; 
;     The ISO formatted date.
; 
;   all_dirs, in, type=strarr
;
;     The directories in which to find scans for which the line defs
;     are needed.
;
;   workdir, in, type=string
; 
;     The directory under which the download directory is to be
;     found or created. 
; 
; 
; :History:
; 
;   2018-06-08 : MGL. First version.
; 
;-
pro red_download_linedefs, isodate, all_dirs, workdir

  dir = workdir+'/downloads/'  

  dotdate = red_strreplace(isodate, '-', '.', n = 2)
  ldpath = 'http://www.royac.iac.es/Logfiles/CHROMIS/linedef/'
  downloadOK = red_geturl(ldpath , contents=ldfiles )
  if downloadOK then begin
    todays_linedefs = ldfiles[where(strmatch(ldfiles, 'linedef.py-'+dotdate+'*', /FOLD_CASE) EQ 1)]
    todays_linedef_times = todays_linedefs
    for i=0,n_elements(todays_linedef_times)-1 do begin
      todays_linedef_times[i] = red_time_conv( strmid(todays_linedef_times[i],strlen(todays_linedef_times[i])-8,8) )
    endfor
    ntlt = n_elements( todays_linedef_times )
    if ntlt gt 0 then begin 
      for i=0,n_elements(all_dirs)-1 do begin
        red_append, all_times, red_Time_conv( strmid(all_dirs[i],strlen(all_dirs[i])-8,8) )
      endfor
      data_times = all_times[uniq(all_times, sort(all_times))]
      for i=0,n_elements(data_times)-1 do begin
        lastidx = max(where( todays_linedef_times lt data_times[i] ))
        if lastidx lt 0 then begin ; TODO: if no calib exists, get one from previous day?
          print, 'red::download : There is not calibration done before the data: ' + $
                 all_dirs[ where(all_times eq data_times[i]) ]
          print, 'red::download : It is safe to ignore this warning if this dataset will not be used for science ' + $
                 ' (it might be an early-morning test or calibration).'
          continue
        endif else begin
          red_append, used_linedefs, todays_linedef_times[lastidx]
        endelse
      endfor
      
      used_linedefs = used_linedefs[uniq(used_linedefs, sort(used_linedefs))]
      todays_linedefs = todays_linedefs[ where(todays_linedef_times eq used_linedefs) ]
      for i=0,n_elements(todays_linedefs)-1 do begin
        downloadOK = red_geturl( ldpath+todays_linedefs[i], dir=dir, overwrite = overwrite )
      endfor
    endif
  endif
  
end
