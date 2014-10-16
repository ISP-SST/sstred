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
;       Lower time limit for plot (in hours after midnight). If given
;       as a variable, on output it contains the real limit, e.g, as
;       detemined by keyword onlydata. 
;
;    tmax : in, out, optional, type=scalar
;
;       Upper time limit for plot (in hours after midnight). If given
;       as a variable, on output it contains the real limit, e.g, as
;       detemined by keyword onlydata.
;
;    onlydata : in, optional, type=boolean
;
;       Limit time interval to when data were collected.
;
;    dt_mean : in, optional, type=scalar, default=3
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
;     2014-08-18 : MGL. Use red_strreplace(), not strreplace(). 
;  
;     2014-10-10 : TL. Fixed an issue with calculating the isodate
;                  string.  
; 
;     2014-10-10 : MGL. 8x8 r0 data does not exist before 2013-10-28.
;                  Set plot24 to true to replace plot8 for earlier
;                  dates. Make sure to keep within the plot range when
;                  marking directories. All cg command to just load
;                  themselves, not add, and then execute only when
;                  writing to file. Mark TRIPPEL data directories with
;                  their own color (pink).
;
;     2014-10-10 : MGL. Bugfix: do not assume the r0time array is a 2D
;                  array. 
;
;     2014-10-16 : MGL. Store time intervals of the different
;                  directories to disk to speed up next plot. Encode
;                  plot interval in plot file name so we can do
;                  multiple plots without overwriting. Make tick mark
;                  labels as properly formatted times. Accept string
;                  values for tmin and tmax.
;
;
;
;-
pro red_plot_r0, dir = dir, today = today, date = date $
                 , lapalma = lapalma $
                 , local = local, r0max = r0max, extension = extension $
                 , tmin = tmin, tmax = tmax $
                 , plotvertical = plotvertical $
                 , plot24 = plot24, plot8 = plot8 $
                 , markdata = markdata, onlydata = onlydata $
                 , dt_mean = dt_mean


  if n_elements(r0max) eq 0 then r0max = 0.40             
  if n_elements(extension) eq 0 then extension = 'jpg'  
  if n_elements(plotvertical) eq 0 then plotvertical = 0
  if n_elements(dt_mean) eq 0 then dt_mean = 3.
  if n_elements(tmax) eq 0 then tmax = 19
  if n_elements(tmin) eq 0 then tmin =  7
 
  ;; Colors to use for plotting
  color_8       = 'black'
  color_24      = 'goldenrod'
  color_dark    = 'black'
  color_flat    = 'gray'
  color_polcal  = 'magenta'
  color_pinhole = 'green'

  color_blue    = 'blue'
  color_red     = 'red'
  color_spec    = 'pink'

  ;; Are tmin or tmax strings?
  if n_elements(tmin) eq 1 then if size(tmin, /tname) eq 'STRING' then tmin = red_time2double(tmin)/3600.
  if n_elements(tmax) eq 1 then if size(tmax, /tname) eq 'STRING' then tmax = red_time2double(tmax)/3600.

  ;; Did we specify a directory somehow?
  if keyword_set(local) then dir = './' ; "local" keyword highest prio
  if n_elements(dir) eq 0 then begin    ; Given "dir" second prio
     if keyword_set(lapalma) then begin ; La Palma data tree
        dir = '/data/'                  ;
     endif else begin                   ; Root
        dir = '/'                       ;
     endelse                            ;
  endif

  timeregex = '[0-2][0-9]:[0-5][0-9]:[0-6][0-9]'
;  dateregex = '[0-9][0-9][0-9][0-9][.-][0-9][0-9][.-][0-9][0-9]'
  dateregex = '20[0-2][0-9][.-][01][0-9][.-][0-3][0-9]'

  ;; Find the directory/directories to plot.
  
  ;; Does the directory name contain a time-stamp?
;  regex = '/[0-9][0-9]:[0-9][0-9]:[0-9][0-9]'
  pos = stregex(dir,'/'+timeregex)      
  if pos ne -1 then begin

     ;; Yes, this is a time-stamped directory.

     dnames = [dir]

  endif else begin

     ;; No, this is not a time-stamped directory.
     ;; Does it contain time-stamped directories?

     dnames = file_search(dir+'/'+timeregex, count = Ndirs)
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
        pos = stregex(dir,'/'+timeregex)      
        if pos ne -1 then begin

           ;; Yes, this is a time-stamped directory.

           dnames = [dir]

        endif else begin

           ;; No, this is still not a time-stamped directory.
           ;; But does it contain time-stamped directories this time?

           dnames = file_search(dir+'/'+timeregex, count = Ndirs)
       endelse
        
     endif

  endelse



  ;; Do we have a date?
  if keyword_set(today) then date = systime(/UTC,/JULIAN)
  if n_elements(date) eq 0 then begin
     if n_elements(dir) gt 0 then begin
        print, 'red_plot_r0 : Try to get date from the selected directory.'
        pos = stregex(dir,'/'+dateregex)
        if pos ne -1 then begin
           date = strmid(dir, pos+1, 10)
           datedir = strmid(dir, 0, pos+12)
        endif
     endif
  endif
  if n_elements(date) eq 0 then begin
     print, 'red_plot_r0 : Try to extract a date from the current directory.'
     pwd = getenv('PWD')
     pos = stregex(pwd,'/'+dateregex)
     if pos ne -1 then begin
        date = strmid(pwd, pos+1, 10)
     endif
  endif
  if n_elements(date) ne 0 then begin
     dashdate=red_strreplace(date, '.', '-', n = 2)
     isodate = (strsplit(date_conv(dashdate[0], 'F'), 'T', /extract))[0]
     ;; isodate = (strsplit(date_conv(red_strreplace(date, '.', '-', n = 2), 'F'), 'T', /extract))[0]
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


  thisdateregex = string(yr, format = '(i04)') + '[.-]' + string(mo, format = '(i02)') $
              + '[.-]' + string(dy, format = '(i02)')
  
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
  Nr0 = (size(r0data, /dim))[0]
  if keyword_set(plot8) and Nr0 eq 1 then begin
     print, 'red_plot_r0 : 8x8 r0 data does not exist before 2013-10-28.'
     print, '              Using /plot24 instead.'
     plot8  = 0
     plot24 = 1
  endif


  analysis_dir = 'dir-analysis/'+isodate+'/'

  if file_test(analysis_dir+'timedirs.txt') then begin

     openr, flun, /get_lun, analysis_dir+'timedirs.txt'
     tmp = ''
     readf, flun, tmp
     Ntd = long((strsplit(tmp,' ',/extr))[1])
     timedirs = strarr(Ntd)
     readf, flun, timedirs
     free_lun, flun

  endif else begin

     ;; Find the timestamp directories
     timedirs = red_find_matching_dirs(timeregex, rootdir = datedir, count = Ntd)
     if Ntd gt 0 then timedirs = red_strreplace(timedirs,datedir,'')

     file_mkdir, analysis_dir

     openw, flun, /get_lun, analysis_dir+'timedirs.txt'
     printf, flun, datedir, ' '+strtrim(Ntd, 2)
     printf, flun, red_strreplace(timedirs,datedir,''), format='(a0)'
     free_lun, flun

  endelse

  print, 'red_plot_r0 : Analyze the timestamp directories'
  for i = 0, Ntd-1 do begin
     
     print, timedirs[i]
     
     intfile = analysis_dir+red_strreplace(timedirs[i],datedir,'')+'/interval.txt'

     if ~file_test(intfile) then begin

        tdir = datedir+red_strreplace(timedirs[i],'/./','/')
;        print, 'cd '+tdir+' ; ls'
        spawn, 'cd '+tdir+' ; ls', cdirs
        
        if n_elements(cdirs) gt 0 then begin
;           print, 'cd '+tdir+'/'+cdirs[0]+' ; ls -rt '
           spawn, 'cd '+tdir+'/'+cdirs[0]+' ; ls -rt ', fnames
           
           Nf = n_elements(fnames)
           
           if Nf gt 0 then begin

              head = fzhead(tdir+'/'+cdirs[0]+'/'+fnames[0])
              istart = strpos(head, 'Ts=')
              istop = strpos(head, 'Te=')
              len = istop - istart
              tstart = total(double(strsplit((strsplit(strmid(head, istart, len) $
                                                       , ' ', /extr))[1],':',/extr)) * [3600.,60., 1.])
              
              head = fzhead(tdir+'/'+cdirs[0]+'/'+fnames[Nf-1])
              istart = strpos(head, 'Ts=')
              istop = strpos(head, 'Te=')
              len = istop - istart
              tstop = total(double(strsplit((strsplit(strmid(head, istart, len) $
                                                      , ' ', /extr))[1],':',/extr)) * [3600.,60., 1.])

              tinterval = [tstart, tstop]/3600.

              file_mkdir, analysis_dir+red_strreplace(timedirs[i],datedir,'')

              openw, flun, /get_lun, intfile
              printf, flun, tinterval, format = '(f7.4)'
              free_lun, flun

           endif                ; Nf
        endif                   ; cdirs
     endif                      ; intfile
  endfor                        ; i

  print, 'red_plot_r0 : Find start and stop times for each timestamp directory'
  tmin_data = 24.
  tmax_data = 0.
  tmin_all = 24.
  tmax_all = 0.
  for i = 0, Ntd-1 do begin

     intfile = analysis_dir+red_strreplace(timedirs[i],datedir,'')+'/interval.txt'
        
     if file_test(intfile) then begin
        openr, flun, /get_lun, intfile
        tinterval = fltarr(2)
        readf, flun, tinterval
        free_lun, flun
        
        if ~strmatch(timedirs[i],'*dark*',/fold) $
           and ~strmatch(timedirs[i],'*flat*', /fold) $
           and ~strmatch(timedirs[i],'*pinh*', /fold) $
           and ~strmatch(timedirs[i],'*polcal*', /fold) then begin
           
           ;; This is a data directory!
           
           tmin_data = min([tmin_data, tinterval[0]])
           tmax_data = max([tmax_data, tinterval[1]])
        endif                   ; data?

        tmin_all = min([tmin_all, tinterval[0]])
        tmax_all = max([tmax_all, tinterval[1]])

     endif                      ; intfile
  
  endfor                        ; i

;  if n_elements(tmin) eq 0 then tmin = tmin_all - 0.1
;  if n_elements(tmax) eq 0 then tmax = tmax_all + 0.1
 
  result = label_date(date_format = '%H:%I')

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

           datedir = file_search(dir+'/'+thisdateregex, count = Ndd)
           datedir = red_find_matching_dirs(thisdateregex, rootdir = dir, count = Ndd, maxdepth = 3)

        endif else Ndd = 1

        if Ndd gt 0 then begin

           if Ntd gt 0 then begin

              if keyword_set(onlydata) then begin
                 ;; Set up the plot.
                 cgplot, /add, [tmin_data-0.2, tmax_data+0.2]*3600., [0, r0max] $
                         , /nodata, /ystyle $
                         , xtitle = 't (UT)', ytitle = 'r$\sub0$ / 1 m', title = isodate $
                         , xstyle = 5, xrange = [tmin_data-0.2, tmax_data+0.2]*3600.
              endif else begin
                 ;; Set up the plot.
                 cgplot, /add, [tmin, tmax]*3600., [0, r0max] $
                         , /nodata, /ystyle $
                         , xtitle = 't (UT)', ytitle = 'r$\sub0$ / 1 m', title = isodate $
                         , xstyle = 5, xrange = [tmin, tmax]*3600.
              endelse
     
              for i = 0, Ntd-1 do begin

                 ;; Now do the marking of directories
                 
                 intfile = analysis_dir+red_strreplace(timedirs[i],datedir,'')+'/interval.txt'
                 
                 if file_test(intfile) then begin
                    openr, flun, /get_lun, intfile
                    tinterval = fltarr(2)
                    readf, flun, tinterval
                    free_lun, flun
                    
                    ;; Proceed only if this directory is at least
                    ;; partly within the plot range. 
                    if (tinterval[0]*3600 ge !x.crange[0] and tinterval[0]*3600 le !x.crange[1]) or $
                       (tinterval[1]*3600 ge !x.crange[0] and tinterval[1]*3600 le !x.crange[1]) then begin
                       
                       ;; Take plot range into account
                       tinterval = [tinterval[0] >!x.crange[0]/3600, tinterval[1] <!x.crange[1]/3600]
                       
                       xpoly = tinterval[[0, 1, 1, 0, 0]]
                       ypoly = [0, 0, 1, 1, 0]*0.01

                       if strmatch(timedirs[i],'*dark*',/fold) then begin
                          ;; Darks
                          cgwindow, 'cgColorFill', /loadcmd, xpoly*3600, (ypoly*2+0.01), color = color_dark
                       endif else if strmatch(timedirs[i],'*flat*', /fold) then begin
                          ;; Flats
                          cgwindow, 'cgColorFill', /loadcmd, xpoly*3600, (ypoly*2+0.01), color = color_flat
                       endif else if strmatch(timedirs[i],'*pinh*', /fold) or $
                          strmatch(timedirs[i], '*[Gg]rid*', /fold) then begin
                          ;; Pinholes or TRIPPEL grid
                          cgwindow, 'cgColorFill', /loadcmd, xpoly*3600, (ypoly*2+0.01), color = color_pinhole
                       endif else if strmatch(timedirs[i],'*polcal*', /fold) then begin
                          ;; Polcal
                          cgwindow, 'cgColorFill', /loadcmd, xpoly*3600, (ypoly*2+0.01), color = color_polcal
                       endif else begin
                          
                          ;; Is this TRIPPEL data?
                          spnames = file_search(timedirs[i]+'Spec*', count = Nspec)
                          if Nspec gt 0 then begin
                             
                             ;; TRIPPEL data
                             cgwindow, 'cgColorFill', /loadcmd, xpoly*3600, (ypoly*2+0.01), color = color_spec
                             
                          endif else if strmatch(timedirs[i],'*blue*', /fold) then begin
                             ;; Blue data
                             cgwindow, 'cgColorFill', /loadcmd, xpoly*3600, (ypoly+0.010), color = color_blue
                          endif else begin
                             ;; Red data
                             cgwindow, 'cgColorFill', /loadcmd, xpoly*3600, (ypoly+0.020), color = color_red
                          endelse ; red/blue
                       endelse    ; data type
                    endif         ; within range
                 endif            ; intfile
              endfor              ; i
           endif                  ; Ntd
        endif                     ; Ndd
        
     endif else begin           ; markdata
        ;; Set up the plot.
        cgplot, /add, [tmin, tmax]*3600, [0, r0max], /nodata, /ystyle $
                , xtitle = 't (UT)', ytitle = 'r$\sub0$ / 1 m', title = isodate 
     endelse

     cgwindow, /add, 'mgl_timeaxis', color = 'opposite' 

     ;; Return time limits of plot
     tmin = !x.crange[0]/3600.
     tmax = !x.crange[1]/3600.

     ;; Flatline r0 values outside these limits:
     r0_top = r0max * 0.99
     r0_bot = 0.040
     
     ;; Smooth over a dt_mean*60. wide window because r0time has a 1
     ;; sec step.

     ;; Do the r0 plotting
     if keyword_set(plot24) then begin
        indx = where(r0data[0, *] gt r0_top, Ntop)
        cgwindow, 'cgplot', /loadcmd, /over, r0time, r0data[0, *] >r0_bot<r0_top $
                  , color = color_24, psym = 16, symsize = 0.1
        if Ntop gt 0 then cgwindow, 'cgplot', /loadcmd, /over, r0time[indx], replicate(r0_top, Ntop) $
                                    , color = color_24, psym = 5, symsize = 0.15
        cgwindow, 'cgplot', /loadcmd, /over, r0time, smooth(r0data[0, *], dt_mean*60.) >r0_bot $
                  , color = color_24
     endif
     if keyword_set(plot8) then begin
        indx = where(r0data[1, *] gt r0_top, Ntop)
        cgwindow, 'cgplot', /loadcmd, /over, r0time, r0data[1, *] >r0_bot<r0_top $
                  , color = color_8, psym = 16, symsize = 0.1
        if Ntop gt 0 then cgwindow, 'cgplot', /loadcmd, /over, r0time[indx], replicate(r0_top, Ntop) $
                                    , color = color_8, psym = 5, symsize = 0.15
        cgwindow, 'cgplot', /loadcmd, /over, r0time, smooth(r0data[1, *], dt_mean*60.) >r0_bot $
                  , color = color_8
     endif

     pname = 'r0_' + isodate + '_' $
             + string(long(tmin),format='(i02)')+':'+string(60*(tmin-long(tmin)),format='(i02)') $
             + '-' $
             + string(long(tmax),format='(i02)')+':'+string(60*(tmax-long(tmax)),format='(i02)') $
             + '.pdf'

     cgcontrol, output = pname

  endif                         ; plot24 or plot8

end
