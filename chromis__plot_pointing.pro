; docformat = 'rst'

;+
; Plot pointing during a day.
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
;    file : in, optional, type=varies
;   
;       If a string, then plot to a file with this name, else if TRUE
;       then print to the default file "dir-analysis/pointing.png".
; 
; 
; :History:
; 
;    2017-10-10 : MGL. First version.
; 
; 
; 
; 
;-
pro chromis::plot_pointing, file = file

  if n_elements(file) gt 0 then begin
    if size(1,/tname) ne 'STRING' && keyword_set(file) then begin
      cgPS_Open, 'dir-analysis/pointing.png', /decomposed
    endif else begin
      cgPS_Open, file, /decomposed
    endelse
  endif
  
  psym_calib = 16
  psym_data = 16

  symsize_calib = 0.3
  symsize_data = .5

  ;; Parameters to convert from time to wavelengths to plot colors.
  col_tbeg = 7                  ; hours
  col_tend = 19                 ; hours
  col_wbeg = 380.               ; nm
  col_wend = 700.               ; nm
  
  ;; Get pointing data
  red_logdata, self.isodate, time_pig, pig = metadata_pig, rsun = rsun

  red_plot_diskcoordinates, metadata_pig, rsun=rsun, rplot=1.03, color = 'gray' $
                            , title = 'SST '+self.isodate + ' pointing'

  dirs = file_search('dir-analysis/'+self.isodate+'/*', count = Ndirs)

  for idir = 0, Ndirs-1 do begin

    subdirs = file_search(dirs[idir]+'/??:??:??', count = Nsubdirs)
    
    if Nsubdirs gt 0 then begin
    
      for isubdir = 0, Nsubdirs-1 do begin

        openr, lun, /get_lun, subdirs[isubdir]+'/interval.txt'
        readf, lun, tbeg
        readf, lun, tend
        free_lun, lun

        indx = where(time_pig ge tbeg*3600 and time_pig le tend*3600, count)

        if count gt 0 then begin

          
          if file_basename(dirs[idir]) eq 'CHROMIS-data' $
             or file_basename(dirs[idir]) eq 'Science' then begin
            psym = psym_data
            symsize = symsize_data
            col_wavelengths = (time_pig[indx]/3600.-col_tbeg)/(col_tend-col_tbeg) $
                              * (col_wend-col_wbeg) + col_wbeg
            colors = red_WavelengthToRGB(col_wavelengths, /num)
          endif else begin
            psym = psym_calib
            symsize = symsize_calib
            colors = 0
          endelse
        endif
        red_plot_diskcoordinates, /over, metadata_pig, indx = indx, rsun=rsun, rplot=1.01 $
                                  , color = colors, psym = psym, symsize = symsize
      endfor
      
    endif
  endfor                ; idir
  
  if n_elements(file) gt 0 then cgPS_Close

end
