; docformat = 'rst'

;+
; Compare pointing according to pig and turret.
;
; Change directory to a CRISP2 workdir and then do .rn red_compare_pigturret.
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
; :History:
; 
;    2023-07-03: MGL. First version.
; 
;-

pro red_compare_pigturret, dum

  spawn, 'pwd', pwd
  date = (file_basename(file_dirname(pwd)))[0]
  
  ;; Cordinates from the turret logfile
  red_logdata, date, time, turret = turret_coord, Rsun = Rsun

  ;; Cordinates from the pig logfile, using turret's time coordinates
  red_logdata, date, time, pig = pig_coord

  ;; Create a spectrum color scale based on the time

  start_time = 8. * 3600.
  stop_time = 13 * 3600.
  indx = where(time ge start_time and time le stop_time, Nwhere)

  start_lambda = 380
  stop_lambda = 700

  if Nwhere eq 0 then stop
  
  lambda =  start_lambda + (stop_lambda - start_lambda) * (time[indx] - start_time) / (stop_time - start_time) 
  colors = red_WavelengthToRGB(lambda, /num)
  
                                ; cgwindow
  cgwindow, 'red_plot_diskcoordinates', turret_coord[*, indx], color = colors $
            , title = date + ' : Turret coordinates', rsun = rsun, psym = 16, rplot = 1.2
  cgcontrol, output = date+'_coord_turret.pdf'
  
                                ; cgwindow
  cgwindow, 'red_plot_diskcoordinates', pig_coord[*, indx], color = colors $
            , title = date + ' : PIG coordinates', rsun = rsun, psym = 16, rplot = 1.2
  cgcontrol, output = date+'_coord_pig.pdf'
  
                                ; cgwindow
  cgwindow, 'red_plot_diskcoordinates', turret_coord[*, indx] - pig_coord[*, indx], color = colors $
            , title = date + ' : (Turret - PIG) diff coordinates', rsun = 1, psym = 16, rplot = 200, unit = '1"'
  cgcontrol, output = date+'_coord_diff.pdf'


end

a = crisp2red(/dev, /no)        ; Change to chromisred or crispred if you want to use a non-crisp2 workdir.
a -> download
red_compare_pigturret

end
