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
;   ismos : in, optional, type=boolean
;   
;     This is a mosaic directory, plot for tiles rather than scans.
; 
;   pname : in, optional, type=string
;
;     File name in which to save the plot.
; 
; :History:
; 
;   2023-08-28 : MGL. New keyword ismos.
; 
;-
pro red_plot_r0_stats, states $
                       , ismos = ismos $
                       , pname = pname

  inam = red_subprogram(/low, calling = inam1)

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

  if keyword_set(ismos) then begin
    Nscans = 1
    amos = reform((stregex(states.filename, 'mos([0-9][0-9])', /extract, /subexpr))[1, *])
    umos = amos[uniq(amos,sort(amos))]
    Nmos = n_elements(umos)
    Npoints = Nmos
  endif else begin
    uscan = states[uniq(states.scannumber,sort(states.scannumber))].scannumber
    Nscans = n_elements(uscan)
    Npoints = Nscans
  endelse
  
  ;; Download the r0 log file 
  red_logdata, isodate, r0time, r0 = r0data, /use_r0_time

  if n_elements(r0data) eq 0 then begin
    ;; Log file probably does not exist on web server (yet)
    print, inam+' : No r0 log data available, no r0 plot will be generated.'
    return                    
  endif
  
  title = instrument + ' ' + strjoin(upref, ',') + ' ' + isodate + ' ' + timestamp

  
  
  ;; Arrays to hold the stats
  tmin = dblarr(Npoints)
  tmax = dblarr(Npoints)
  r0_8x8_min    = fltarr(Npoints)
  r0_8x8_mean   = fltarr(Npoints)
  r0_8x8_median = fltarr(Npoints)
  r0_8x8_max    = fltarr(Npoints)
  r0_24x24_min    = fltarr(Npoints)
  r0_24x24_mean   = fltarr(Npoints)
  r0_24x24_median = fltarr(Npoints)
  r0_24x24_max    = fltarr(Npoints)

  for ipoint = 0L, Npoints-1 do begin
    
    ;; Time
    if keyword_set(ismos) then begin
      indx = where(amos eq umos[ipoint])
    endif else begin
      indx = where(states.scannumber eq uscan[ipoint])
    end
    tmp = max(states[indx].framenumber, mxl, subscript_min=mnl)

    hmn = red_readhead(states[indx[mnl]].filename)
    tmin[ipoint] = red_time2double((strsplit(fxpar(hmn, 'DATE-BEG'), 'T', /extract))[1])
    hmx = red_readhead(states[indx[mxl]].filename)
    tmax[ipoint] = red_time2double((strsplit(fxpar(hmx, 'DATE-END'), 'T', /extract))[1])

    ;; Statistics
    indx = where(r0time ge tmin[ipoint] and r0time le tmax[ipoint], Nmatch)
    if Nmatch eq 0 then begin
      ;; There might not be any r0 values during the exposure of one
      ;; tile. In that case, pic the closest one.
      tt = (tmax[ipoint]+tmin[ipoint])/2.
      diff = abs(r0time - tt)
      tmp = min(diff, indx)
    endif 
    if size(r0data, /n_dim) eq 2 then begin
      r0_24x24_min[ipoint]    = min(r0data[0, [indx]])    
      r0_24x24_mean[ipoint]   = mean(r0data[0, [indx]])   
      r0_24x24_median[ipoint] = median(r0data[0, [indx]]) 
      r0_24x24_max[ipoint]    = max(r0data[0, [indx]])    
      r0_8x8_min[ipoint]      = min(r0data[1, [indx]])    
      r0_8x8_mean[ipoint]     = mean(r0data[1, [indx]])   
      r0_8x8_median[ipoint]   = median(r0data[1, [indx]]) 
      r0_8x8_max[ipoint]      = max(r0data[1, [indx]])    
    endif else begin
      r0_24x24_min[ipoint]    = min(r0data[[indx]])    
      r0_24x24_mean[ipoint]   = mean(r0data[[indx]])   
      r0_24x24_median[ipoint] = median(r0data[[indx]]) 
      r0_24x24_max[ipoint]    = max(r0data[[indx]])    
    endelse
    
    
  endfor                        ; ipoint

  if keyword_set(ismos) then begin
    xtitle = 'tile #'
    xrange = [-0.1, Nmos-1+0.1]
  endif else begin
    xtitle = 'scan #'
    xrange =  [min(states.scannumber), max(states.scannumber)]
  endelse 
  cgwindow, 'cgplot', /nodata, [0], [0] $
            , xrange = xrange $
            , yrange = [0, r0max] $
            , xtitle = xtitle, ytitle = 'r$\sub0$ / 1 m' $
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

  if keyword_set(ismos) then begin
    xpoints = long(umos)
    xpoly = long([umos, reverse(umos), umos[0]])
  endif else begin
    xpoints = uscan
    xpoly = [uscan, reverse(uscan), uscan[0]]
  endelse
  
  ypoly = [r0_24x24_min, reverse(r0_24x24_max), r0_24x24_min[0]]
  cgwindow, /add, 'cgcolorfill', xpoly, ypoly, color = 'Light Gray'

  ypoly = [r0_8x8_min, reverse(r0_8x8_max), r0_8x8_min[0]]
  cgwindow, /add, 'cgcolorfill', xpoly, ypoly, color = 'Pale Goldenrod'

  cgwindow, /add, /over, 'cgplot', xpoints, r0_8x8_mean $
            , linestyle = statslines[0], color = statscolors[1]
  cgwindow, /add, /over, 'cgplot', xpoints, r0_8x8_median $
            , linestyle = statslines[0], color = statscolors[2]
  cgwindow, /add, /over, 'cgplot', xpoints, r0_24x24_mean $
            , linestyle = statslines[1], color = statscolors[1]
  cgwindow, /add, /over, 'cgplot', xpoints, r0_24x24_median $
            , linestyle = statslines[1], color = statscolors[2]

  if ~ismos then begin
    delta = 5
    for ii = delta, uscan[-1]-1, delta do begin
      cgwindow, /add, /over, 'cgplot', [ii, ii], !y.crange, color = 'gray', linestyle = 1
    endfor
  endif
  
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
