; docformat = 'rst'

;+
; Plots r0 for SST data directories.
; 
; :Categories:
;
;    SST observations
; 
; 
; :author:
; 
;    Mats Löfdahl, 2014-05-28
;
; :Keywords:
; 
;    dir : in, optional, type=string
; 
;      The directory where the data is.
; 
;    date : in, optional
;   
;      The date of observations. Any type if input understood by
;      date_conv.pro plus pseudo-ISO ('YYYY.MM.DD').
;   
;    today : in, optional, type=boolean
;   
;      Sets date to today. And tries to download the latest r0 file.
;   
;    lapalma : in, optional, type=boolean 
; 
;      If you set this and do not specify a dir, then the program will
;      look for data daily directories in places where they can
;      usually be found at the SST site in La Palma,
;      "/data/store[1234]/" and "/data/disk[1234]/*/"
;
;
;    r0max : in, optional, type=scalar, default=0.40
;
;      Upper plot r0 limit in m.
;
;    extension : in, optional, type=string, default='jpg'
;
;      Set this to an extension for the plot files, like 'png', 'pdf'.
;
;    plotvertical : in, optional, type=boolean
;
;      Set this to plot vertical lines between scans. Such lines will
;      be plotted by default if all scans are with the same prefilter.
;
;    plot24 : in, optional, type=boolean
;
;       Set this to plot the full-day, large-FOV r0 value.
;
;    plot8  : in, optional, type=boolean
;
;       Set this to plot the full-day, small-FOV r0 value.
;
;    markdata : in, optional, type=boolean
;
;       Set this to mark when data were collected in the full-day r0
;       plot. Requires that the directory where data from the day is
;       known or can be deduced.
;
;    tmin : in, out, optional, type=scalar
;
;       Lower time limit for plot. If given as a variable, on output
;       it contains the real limit, e.g, as detemined by keyword
;       onlydata. 
;
;    tmax : in, out, optional, type=scalar
;
;       Upper time limit for plot. If given as a variable, on output
;       it contains the real limit, e.g, as detemined by keyword
;       onlydata. 
;
;    onlydata : in, optional, type=boolean
;
;       Limit time interval to when data were collected.
;
;    dt_mean : in, optional, type=scalar, default=1
;
;       Time window for running mean in minutes. 
;
;
; :History:
; 
;     2014-05-28 : MGL. Start working on it, based loosely on earlier
;                  routines sst_plotonedir and sst_plotdirs. It can
;                  now plot r0 for an entire day and mark time
;                  intervals when data of various kinds were collected
;                  (red and blue, but also darks, flats, polcal).
;
;     2014-06-09 : MGL. Speed up directory search by using new routine
;                  red_find_matching_dirs.
; 
;-
pro red_plot_r0, dir = dir, today = today, date = date $
                 , lapalma = lapalma $
                 , local = local, r0max = r0max, extension = extension $
                 , tmin = tmin, tmax = tmax $
                 , plotvertical = plotvertical $
                 , plot24 = plot24, plot8 = plot8, markdata = markdata, onlydata = onlydata $
                 , dt_mean = dt_mean


  if n_elements(r0max) eq 0 then r0max = 0.40             
  if n_elements(tmax) eq 0 then tmax = 19
  if n_elements(tmin) eq 0 then tmin =  7
  if n_elements(extension) eq 0 then extension = 'jpg'  
  if n_elements(plotvertical) eq 0 then plotvertical = 0
  if n_elements(dt_mean) eq 0 then dt_mean = 1.

  ;; Colors to use for plotting
  color_8       = 'black'
  color_24      = 'goldenrod'
  color_dark    = 'black'
  color_flat    = 'gray'
  color_polcal  = 'magenta'
  color_pinhole = 'green'
  color_blue    = 'blue'
  color_red     = 'red'


  ;; Did we specify a directory somehow?
  if keyword_set(local) then dir = './' ; "local" keyword highest prio
  if n_elements(dir) eq 0 then begin    ; Given "dir" second prio
     if keyword_set(lapalma) then begin ; La Palma data tree
        dir = '/data/'                  ;
     endif else begin                   ; Root
        dir = '/'                       ;
     endelse                            ;
  endif

  ;; Find the directory/directories to plot.
  
  ;; Does the directory name contain a time-stamp?
  regex = '/[0-9][0-9]:[0-9][0-9]:[0-9][0-9]'
  pos = stregex(dir,regex)      
  if pos ne -1 then begin

     ;; Yes, this is a time-stamped directory.

     dnames = [dir]

  endif else begin

     ;; No, this is not a time-stamped directory.
     ;; Does it contain time-stamped directories?

     dnames = file_search(dir+'/[0-9][0-9]:[0-9][0-9]:[0-9][0-9]', count = Ndirs)
     if Ndirs eq 0 then begin

        ;; No there are no time stamped directories here

        print
        print, 'Please navigate to and select a directory.'
        print, 'If this is a time-stamped directory, it will be plotted.'
        print, 'If not, you will be asked to select one or several time stamped directories'
        print, 'in the selected directory.'
        dialog_title = 'Please select a directory.'
        pickdir = dialog_pickfile(/dir, path = dir, /must_exist, title = dialog_title)
        if strlen(pickdir) gt 0 then dir = pickdir

        ;; Does this directory name contain a time-stamp?
        pos = stregex(dir,regex)      
        if pos ne -1 then begin

           ;; Yes, this is a time-stamped directory.

           dnames = [dir]

        endif else begin

           ;; No, this is still not a time-stamped directory.
           ;; But does it contain time-stamped directories this time?

           dnames = file_search(dir+'/[0-9][0-9]:[0-9][0-9]:[0-9][0-9]', count = Ndirs)
       endelse
        
     endif

  endelse



  ;; Do we have a date?
  if keyword_set(today) then date = systime(/UTC,/JULIAN)
  if n_elements(date) eq 0 then begin
     if n_elements(dir) gt 0 then begin
        print, 'red_plot_r0 : Try to get date from the selected directory.'
        pos = stregex(dir,'/[0-9][0-9][0-9][0-9][.-][0-9][0-9][.-][0-9][0-9]')
        if pos ne -1 then begin
           date = strmid(dir, pos+1, 10)
           datedir = strmid(dir, 0, pos+11)
        endif
     endif
  endif
  if n_elements(date) eq 0 then begin
     print, 'red_plot_r0 : Try to extract a date from the current directory.'
     pwd = getenv('PWD')
     pos = stregex(pwd,'/[0-9][0-9][0-9][0-9][.-][0-9][0-9][.-][0-9][0-9]')
     if pos ne -1 then begin
        date = strmid(pwd, pos+1, 10)
   endif
  endif
  if n_elements(date) ne 0 then begin
     isodate = (strsplit(date_conv(strreplace(date, '.', '-', n = 2), 'F'), 'T', /extract))[0]
  endif else begin
     print, 'red_plot_r0 : No date given.'
     help, date, dnames[0]
     stop
  endelse

  
  ;; We should now have the date in ISO (YYYY-MM-DD) form:
  print, 'red_plot_r0 : The date is '+isodate
  yr = long((strsplit(isodate, '-', /extract))[0])
  mo = long((strsplit(isodate, '-', /extract))[1])
  dy = long((strsplit(isodate, '-', /extract))[2])


  if Ndirs eq 0 then begin
     print
     print, 'red_plot_r0 : Cannot find any time-stamped directories here.'
  endif else begin
     ;; At this point we should have a list of time-stamped directories
     ;; in dnames.
     print, dnames
  endelse
  
  ;; Download the r0 log file if necessary
  red_logdata, isodate, r0time, r0 = r0data


  ;; Find the timestamp directories
  timeregex = '[0-2][0-9]:[0-5][0-9]:[0-6][0-9]'
  timedirs = red_find_matching_dirs(timeregex, rootdir = datedir, count = Ntd)


  ;; Analyze the timestamp directories
  print, 'red_plot_r0 : Find start and stop times for each timestamp directory'
  tmin_data = 24.
  tmax_data = 0.
  tstarts = dblarr(Ntd)
  tstops = dblarr(Ntd)
  for i = 0, Ntd-1 do begin

     print, timedirs[i]

     tdir = strreplace(timedirs[i],'/./','/')
     spawn, 'cd '+tdir+' ; ls', cdirs
     
     if n_elements(cdirs) gt 0 then begin
        spawn, 'cd '+tdir+'/'+cdirs[0]+' ; ls -rt ', fnames

        Nf = n_elements(fnames)

        if Nf gt 0 then begin
           head = fzhead(tdir+'/'+cdirs[0]+'/'+fnames[0])
           istart = strpos(head, 'Ts=')
           istop = strpos(head, 'Te=')
           len = istop - istart
           tstarts[i] = total(double(strsplit((strsplit(strmid(head, istart, len) $
                                                        , ' ', /extr))[1],':',/extr)) * [3600.,60., 1.])

           head = fzhead(tdir+'/'+cdirs[0]+'/'+fnames[Nf-1])
           istart = strpos(head, 'Ts=')
           istop = strpos(head, 'Te=')
           len = istop - istart
           tstops[i] = total(double(strsplit((strsplit(strmid(head, istart, len) $
                                                       , ' ', /extr))[1],':',/extr)) * [3600.,60., 1.])
           if ~strmatch(timedirs[i],'*dark*',/fold) $
              and ~strmatch(timedirs[i],'*flat*', /fold) $
              and ~strmatch(timedirs[i],'*pinh*', /fold) $
              and ~strmatch(timedirs[i],'*polcal*', /fold) then begin
              
              ;; This is a data directory!
              
              tmin_data = min([tmin_data, tstarts[i]/3600.])
              tmax_data = max([tmax_data, tstops[i]/3600.])
              print, 'red_plot_r0 : New tmin_data,tmax_data:', tmin_data, tmax_data

           endif                ; data type
           
        endif                   ; Nf
     endif                      ; cdirs
  endfor                        ; i



  if keyword_set(plot24) or keyword_set(plot8) $
     or ((keyword_set(markdata) or keyword_set(onlydata)) and n_elements(dir) gt 0) then begin

     ;; Plot r0 for the entire (or part of) day. If we have
     ;; time-stamped directories, possibly mark in the plot where they
     ;; start and end.
     
     ;; Set up the plotting window for the day plot
     cgwindow

     if keyword_set(markdata) or keyword_set(onlydata) then begin
        ;; Find the data directory for this day

        if n_elements(datedir) eq 0 then begin
           dateregex = string(yr, format = '(i04)')+'[.-]'+string(mo, format = '(i02)')+'[.-]'+string(dy, format = '(i02)')
           datedir = file_search(dir+'/'+dateregex, count = Ndd)
           datedir = red_find_matching_dirs(dateregex, rootdir = dir, count = Ndd, maxdepth = 3)

        endif else Ndd = 1

        if Ndd gt 0 then begin

           if Ntd gt 0 then begin

              if keyword_set(onlydata) then begin
                 ;; Set up the plot.
                 cgplot, /add, [tmin_data-0.2, tmax_data+0.2], [0, r0max] $
                         , /nodata, /ystyle, /xstyle $
                         , xtitle = 't (UT)', ytitle = 'r$\sub0$ / 1 m', title = isodate 
              endif else begin
                 ;; Set up the plot.
                 cgplot, /add, [tmin, tmax], [0, r0max] $
                         , /nodata, /xstyle, /ystyle $
                         , xtitle = 't (UT)', ytitle = 'r$\sub0$ / 1 m', title = isodate 
              endelse


              for i = 0, Ntd-1 do begin

                 ;; Now do the marking of directories
                 
                 if tstarts[i] ne 0.0d then begin

                    tinterval = [tstarts[i], tstops[i]]/3600.
                    xpoly = [tstarts[i], tstops[i], tstops[i], tstarts[i], tstarts[i]]/3600.
                    ypoly = [0, 0, 1, 1, 0]*0.01

                    if strmatch(timedirs[i],'*dark*',/fold) then begin
                       ;; Darks
                       cgwindow, 'cgColorFill', /loadcmd, xpoly, ypoly*2+0.01, color = color_dark
                    endif else if strmatch(timedirs[i],'*flat*', /fold) then begin
                       ;; Flats
                       cgwindow, 'cgColorFill', /loadcmd, xpoly, ypoly*2+0.01, color = color_flat
                    endif else if strmatch(timedirs[i],'*pinh*', /fold) then begin
                       ;; Pinholes
                       cgwindow, 'cgColorFill', /loadcmd, xpoly, ypoly*2+0.01, color = color_pinhole
                    endif else if strmatch(timedirs[i],'*polcal*', /fold) then begin
                       ;; Polcal
                       cgwindow, 'cgColorFill', /loadcmd, xpoly, ypoly*2+0.01, color = color_polcal
                    endif else begin
                       if strmatch(timedirs[i],'*blue*', /fold) then begin
                          ;; Blue data
                          cgwindow, 'cgColorFill', /loadcmd, xpoly, ypoly+0.010, color = color_blue
                       endif else begin
                          ;; Red data
                          cgwindow, 'cgColorFill', /loadcmd, xpoly, ypoly+0.020, color = color_red
                       endelse  ; red/blue
                    endelse     ; data type
                    
                 endif          ; starts[i]
              endfor            ; i
           endif                ; Ntd
        endif                   ; Ndd
        
     endif else begin           ; markdata
        ;; Set up the plot.
        cgplot, /add, [tmin, tmax], [0, r0max], /nodata, /xstyle, /ystyle $
                , xtitle = 't (UT)', ytitle = 'r$\sub0$ / 1 m', title = isodate 
     endelse

     tmin = !x.crange[0]
     tmax = !x.crange[1]

     ;; Flatline r0 values outside these limits:
     r0_top = r0max * 0.99
     r0_bot = 0.040
     
     ;; Smooth over a dt_mean*60. wide window because r0time has a 1
     ;; sec step.

     ;; Do the r0 plotting
     if keyword_set(plot24) then begin
        indx = where(r0data[0, *] gt r0_top, Ntop)
        cgplot, /add, /over, r0time[0, *]/3600., r0data[0, *] >r0_bot<r0_top, color = color_24, psym = 16, symsize = 0.1
        cgplot, /add, /over, r0time[0, indx]/3600., replicate(r0_top, Ntop), color = color_24, psym = 5, symsize = 0.15
        cgplot, /add, /over, r0time[0, *]/3600., smooth(r0data[0, *], dt_mean*60.) >r0_bot, color = color_24
     endif
     if keyword_set(plot8) then begin
        indx = where(r0data[1, *] gt r0_top, Ntop)
        cgplot, /add, /over, r0time[0, *]/3600., r0data[1, *] >r0_bot<r0_top, color = color_8, psym = 16, symsize = 0.1
        cgplot, /add, /over, r0time[0, indx]/3600., replicate(r0_top, Ntop), color = color_8, psym = 5, symsize = 0.15
        cgplot, /add, /over, r0time[0, *]/3600., smooth(r0data[1, *], dt_mean*60.) >r0_bot, color = color_8
     endif

     cgcontrol, output = 'r0_'+isodate+'.pdf'

  endif 

end
