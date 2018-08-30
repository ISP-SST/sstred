; docformat = 'rst'

;+
; Plot r0 statistics.
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
; 
; 
; 
; 
; 
;-
pro red_plot_r0_stats,  states, pname = pname

  r0max = 0.3

  ;; Find first and last files
  tmp = max(states.framenumber, mxl, subscript_min=mnl)


  ;; From header
  h0 = red_readhead(states[mnl].filename)
  date_obs = strsplit(fxpar(h0, 'DATE-OBS'), 'T', /extract)
  isodate = date_obs[0]
  timestamp = date_obs[1]
  instrument = fxpar(h0, 'INSTRUME')

  ;; From states
  upref = states(uniq(states.prefilter, sort(states.prefilter))).prefilter
  Npref = n_elements(upref)
  uscan = states[uniq(states.scannumber,sort(states.scannumber))].scannumber
  Nscans = n_elements(uscan)
  
  ;; Download the r0 log file 
  red_logdata, isodate, r0time, r0 = r0data, /use_r0_time
  
  title = instrument + ' ' + strjoin(upref, ',') + ' ' + isodate + ' ' + timestamp

  ;; Arrays to hold the stats
  tmin = dblarr(Nscans)
  tmax = dblarr(Nscans)
  r0_8x8_min    = fltarr(Nscans)
  r0_8x8_mean   = fltarr(Nscans)
  r0_8x8_median = fltarr(Nscans)
  r0_8x8_max    = fltarr(Nscans)
  r0_24x24_min    = fltarr(Nscans)
  r0_24x24_mean   = fltarr(Nscans)
  r0_24x24_median = fltarr(Nscans)
  r0_24x24_max    = fltarr(Nscans)

  for iscan = 0L, Nscans-1 do begin
     
    ;; Time
    indx = where(states.scannumber eq uscan[iscan])
    tmp = max(states[indx].framenumber, mxl, subscript_min=mnl)

    hmn = red_readhead(states[indx[mnl]].filename)
    tmin[iscan] = red_time2double((strsplit(fxpar(hmn, 'DATE-BEG'), 'T', /extract))[1])
    hmx = red_readhead(states[indx[mxl]].filename)
    tmax[iscan] = red_time2double((strsplit(fxpar(hmx, 'DATE-END'), 'T', /extract))[1])

    ;; Statistics
    indx = where(r0time ge tmin[iscan] and r0time le tmax[iscan], Nmatch)
    if Nmatch gt 0 then begin
      if size(r0data, /n_dim) eq 2 then begin
        r0_24x24_min[iscan]    = min(r0data[0, indx])    
        r0_24x24_mean[iscan]   = mean(r0data[0, indx])   
        r0_24x24_median[iscan] = median(r0data[0, indx]) 
        r0_24x24_max[iscan]    = max(r0data[0, indx])    
        r0_8x8_min[iscan]      = min(r0data[1, indx])    
        r0_8x8_mean[iscan]     = mean(r0data[1, indx])   
        r0_8x8_median[iscan]   = median(r0data[1, indx]) 
        r0_8x8_max[iscan]      = max(r0data[1, indx])    
      endif else begin
        r0_24x24_min[iscan]    = min(r0data[indx])    
        r0_24x24_mean[iscan]   = mean(r0data[indx])   
        r0_24x24_median[iscan] = median(r0data[indx]) 
        r0_24x24_max[iscan]    = max(r0data[indx])    
      endelse
    endif
    
  endfor                        ; iscan

  cgwindow, 'cgplot', /nodata, [0], [0] $
            , xrange =  [min(states.scannumber), max(states.scannumber)] $
            , yrange = [0, r0max] $
            , xtitle = 'scan #', ytitle = 'r$\sub0$ / 1 m' $
            , xstyle = 8 

  ;; Add a time axis on top
  
  L = red_gen_timeaxis([min(tmin), max(tmax)])

  cgwindow, /add, 'cgaxis', xaxis = 1 $
            , /xstyle $
            , xrange = [min(tmin), max(tmax)] $
            , XTICKV=L.tickv, XTICKS=L.ticks, XMIN=L.minor, XTICKNAM=L.name

  cgtext, mean(!x.crange), !y.crange[1]*.03, title, /data, align = 0.5,/add

  statscolors = ['green', 'red', 'blue', 'cyan']
  statslines = [0, 1]
  
  xpoly = [uscan, reverse(uscan), uscan[0]]

  ypoly = [r0_24x24_min, reverse(r0_24x24_max), r0_24x24_min[0]]
  cgwindow, /add, 'cgcolorfill', xpoly, ypoly, color = 'Light Gray'

  ypoly = [r0_8x8_min, reverse(r0_8x8_max), r0_8x8_min[0]]
  cgwindow, /add, 'cgcolorfill', xpoly, ypoly, color = 'Pale Goldenrod'

  cgwindow, /add, /over, 'cgplot', uscan, r0_8x8_mean $
            , linestyle = statslines[0], color = statscolors[1]
  cgwindow, /add, /over, 'cgplot', uscan, r0_8x8_median $
            , linestyle = statslines[0], color = statscolors[2]
  cgwindow, /add, /over, 'cgplot', uscan, r0_24x24_mean $
            , linestyle = statslines[1], color = statscolors[1]
  cgwindow, /add, /over, 'cgplot', uscan, r0_24x24_median $
            , linestyle = statslines[1], color = statscolors[2]


  delta = 5
  for ii = delta, uscan[-1]-1, delta do begin
    cgwindow, /add, /over, 'cgplot', [ii, ii], !y.crange, color = 'gray', linestyle = 1
  endfor

  cglegend, /add, align = 0, location = [.15, .85], linestyle = 0 $
            , title = ['mean', 'median'], colors = statscolors[[1, 2]]

  if n_elements(pname) eq 0 then begin
    pname = 'dir-analysis/' + 'r0stats_'+ strjoin(upref, ',')+ '_' + isodate + '_' $
            + timestamp $
            + '.pdf'
  endif

  dir = file_dirname(pname)
  if dir ne '' then file_mkdir, dir
  
  cgcontrol, output = pname

end
