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
;    instruments : in, optional, type=strarr, default=['chromis']
; 
;      A list of instruments to plot scans and statistics for.  
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
;       Set this to plot the large-FOV r0 value.
;
;    plot8  : in, optional, type=boolean
;
;       Set this to plot the small-FOV r0 value.
;
;    scan24  : in, optional, type=boolean
;
;       Set this to plot the large-FOV r0 for CRISP/CHROMIS scans.
;
;    scan8  : in, optional, type=boolean
;
;       Set this to plot the small-FOV r0 for CRISP/CHROMIS scans.
;
;    plotstats : in, optional, type=boolean 
;
;       Set this to plot r0 statistics vs scan number.
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
;       Time window for running mean in minutes. Do not plot running
;       mean if dt_mean is set to zero.
;
;    psym_plot : in, optional, type=integer, default=16
;
;       Plot symbol to use for plot8 and plot24 plots.
;
;    psym_scan : in, optional, type=integer, default=16
;
;       Plot symbol to use for scan8 and scan24 plots.
;
;    symsize_plot : in, optional, type=float, default=0.25
;
;       Plot symbol size to use for plot8 and plot24 plots.
;
;    symsize_scan : in, optional, type=float, default=1.0
;
;       Plot symbol size to use for scan8 and scan24 plots.
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
;     2015-06-28 : MGL. New keyword "noplot". Fixed bug in calculating
;                  when data directory started. Remove old code that
;                  dealt with time-stamped directories, keyword dir
;                  can now be just the date directory. New keywords
;                  "scan24" and "scan8", added code to plot r0 for
;                  scans.  
;
;     2015-08-10 : MGL. Write all plots to subdir dir-analysis. Now
;                  allowed to give dir without trailing slash. 
;
;     2015-08-12 : MGL. Stop using red_timeaxis.pro to get away from
;                  its many dependencies. Import Pit's gen_timeaxis to
;                  the crispred namespace and use it instead. 
;
;     2015-08-17 : MGL. Can now get directory from config.txt.
;
;     2016-05-19 : MGL. Update for CHROMIS. Read file headers with
;                  red_readhead, this gives us FITS keywords. Works
;                  now with plot8 and plot24. New keywords crispscans
;                  and chromisscans to plot CRISP and CHROMIS scans
;                  separately with scan8 and scan24.
;
;     2016-08-15 : MGL. New keywords allowing change of plot symbols.
;                  Also, can now set dt_mean=0 to remove the smoothed
;                  curves.
;
;     2016-09-07 : MGL. Parts moved to analyze-directories method.
;                  Make plotting work with output from that method.
;
;     2016-09-08 : MGL. New keyword plotstats, implemented plotting r0
;                  statistics vs scan number. Remove unused keywords
;                  crispscans and chromisscans.
;
;     2016-09-12 : MGL. Plotstats plots: add a time axis on top, shade
;                  between min and max. 
;
;     2016-09-27 : MGL. New keyword instruments.
;
;-
pro red_plot_r0, dir = dir $
                 , today = today $
                 , date = date $
                 , lapalma = lapalma $
                 , local = local $
                 , r0max = r0max $
                 , extension = extension $
                 , tmin = tmin $
                 , tmax = tmax $
                 , plotvertical = plotvertical $
                 , plot24 = plot24 $
                 , plot8 = plot8 $
                 , psym_plot = psym_plot $
                 , symsize_plot = symsize_plot $
                 , scan24 = scan24 $
                 , scan8 = scan8 $
                 , psym_scan = psym_scan $
                 , symsize_scan = symsize_scan $
                 , markdata = markdata $
                 , onlydata = onlydata $
                 , dt_mean = dt_mean $
                 , noplot = noplot $
                 , plotstats = plotstats $
                 , instruments = instruments
  


  if n_elements(instruments) eq 0 then instruments = ['chromis']; ['chromis', 'crisp']
  if n_elements(r0max) eq 0 then r0max = 0.40             
  if n_elements(extension) eq 0 then extension = 'jpg'  
  if n_elements(plotvertical) eq 0 then plotvertical = 0
  if n_elements(dt_mean) eq 0 then dt_mean = 3.
  if n_elements(tmax) eq 0 then tmax = 19
  if n_elements(tmin) eq 0 then tmin =  7
 
  if n_elements(psym_plot) eq 0 then psym_plot = 16
  if n_elements(psym_scan) eq 0 then psym_scan = 16
  if n_elements(symsize_plot) eq 0 then symsize_plot = 0.25
  if n_elements(symsize_scan) eq 0 then symsize_scan = .5

  Ninstruments = n_elements(instruments)

  ;; Colors to use for plotting
  color_8       = 'black'
  color_24      = 'goldenrod'
  color_dark    = 'black'
  color_flat    = 'gray'
  color_polcal  = 'magenta'
  color_pinhole = 'green'

  ;; Colors used to mark which instruments were collecting data
  color_blue    = 'blue'
  color_red     = 'red'
  color_spec    = 'pink'
  
  ;; Colors used for scan backgrounds.
  precolors = ['sky blue', 'wheat', 'spring green', 'plum', 'medium gray']

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
        if file_test('config.txt') then begin
           ;; Get root_dir from config file
           spawn, 'grep root_dir config.txt', rd
           if size(rd,/n_dim) gt 0 then begin
              rd = rd[0]
              dir = (strsplit(rd,' =',/extract))[1] 
           endif
        endif
     endelse                            ;
  endif

  print, 'Dir: ', dir

  timeregex = '[0-2][0-9]:[0-5][0-9]:[0-6][0-9]'
  dateregex = '20[0-2][0-9][.-][01][0-9][.-][0-3][0-9]'

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

     print, 'red_plot_r0 : Please run dir-analysis first.'
     retall
     
  endelse


  print, 'red_plot_r0 : Find start and stop times for all data'
  tmin_data = 24.
  tmax_data = 0.
  tmin_all = 24.
  tmax_all = 0.
  for i = 0, Ntd-1 do begin

     print,strtrim(i,2)+'/'+strtrim(Ntd,2)

     intfile = analysis_dir + red_strreplace(timedirs[i],datedir,'') $
               + '/interval.txt'
        
     if file_test(intfile) then BEGIN

        openr, flun, /get_lun, intfile
        tinterval = fltarr(2)
        readf, flun, tinterval
        free_lun, flun
        
        if ~strmatch(timedirs[i],'*dark*',/fold) $
           and ~strmatch(timedirs[i],'*flat*', /fold) $
           and ~strmatch(timedirs[i],'*pinh*', /fold) $
           and ~strmatch(timedirs[i],'*pfscan*', /fold) $
           and ~strmatch(timedirs[i],'*polcal*', /fold) then begin
           
           ;; This is a data directory!
           print,intfile
           
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

  
  ;; Return now if all we wanted was finding directories etc, but no
  ;; actual plots.
  if keyword_set(noplot) then return

  ;; Templates for reading information about the scans:
  r0scantemplate = { version : 1.0 $
                     , datastart : 0L $
                     , delimiter : 32B $
                     , missingvalue : !Values.F_NaN $
                     , commentsymbol : '#' $
                     , fieldcount : 3L $
                     , fieldtypes : [4L, 4L, 4L] $
                     , fieldnames : ['time', 'r0_24x24', 'r0_8x8'] $
                     , fieldlocations : [7L, 21L,  33L] $
                     , fieldgroups : lindgen(3) $
                   }

  ;; Need new template with text filter name
  scantemplate = { version : 1.0 $
                   , datastart : 0L $
                   , delimiter : 32B $
                   , missingvalue : !Values.F_NaN $
                   , commentsymbol : '#' $
                   , fieldcount : 12L $
                   , fieldtypes : [3L, 7L, replicate(4L, 10)] $
                   , fieldnames : ['scanno', 'prefilter', 'tstart', 'tstop' $
                                   , 'r0_24x24_min', 'r0_24x24_mean' $
                                   , 'r0_24x24_median', 'r0_24x24_max' $
                                   , 'r0_8x8_min', 'r0_8x8_mean' $
                                   , 'r0_8x8_median', 'r0_8x8_max'] $
                   , fieldlocations : [0L, 2L, 13L, 22L, 30L, 39L, 48L, 57L $
                                       , 68L, 77L, 86L, 95L] $
                   , fieldgroups : lindgen(12) $
                 }




  if keyword_set(plot24) $
     or keyword_set(plot8) $
     or ((keyword_set(markdata) or keyword_set(onlydata)) and n_elements(dir) gt 0) $
  then begin

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
                 L = red_gen_timeaxis([tmin_data-0.2, tmax_data+0.2]*3600.)
                 cgplot, /add, [tmin_data-0.2, tmax_data+0.2]*3600., [0, r0max] $
                         , /nodata, /ystyle $
                         , xtitle = 't [UT]', ytitle = 'r$\sub0$ / 1 m', title = isodate $
                         , xrange = [tmin_data-0.2, tmax_data+0.2]*3600. $
                         , XTICKV=L.tickv, XTICKS=L.ticks, XMIN=L.minor, XTICKNAM=L.name
              endif else begin
                 ;; Set up the plot.
                 L = red_gen_timeaxis([tmin, tmax]*3600.)
                 cgplot, /add, [tmin, tmax]*3600., [0, r0max] $
                         , /nodata, /ystyle $
                         , xtitle = 't [UT]', ytitle = 'r$\sub0$ / 1 m', title = isodate $
                         , xrange = [tmin, tmax]*3600. $
                         , XTICKV=L.tickv, XTICKS=L.ticks, XMIN=L.minor, XTICKNAM=L.name
              endelse

              for i = 0, Ntd-1 do begin

                 print,strtrim(i,2)+'/'+strtrim(Ntd,2)

                 ;; Now do the marking of directories
                 
                 intfile = analysis_dir + red_strreplace(timedirs[i],datedir,'') $
                           + '/interval.txt'
                 
                 if file_test(intfile) then begin
                    openr, flun, /get_lun, intfile
                    tinterval = fltarr(2)
                    readf, flun, tinterval
                    free_lun, flun
                    
                    ;; Proceed only if this directory is at least
                    ;; partially within the plot range.
                    if (tinterval[0]*3600 ge !x.crange[0] and tinterval[0]*3600 le !x.crange[1]) or $
                       (tinterval[1]*3600 ge !x.crange[0] and tinterval[1]*3600 le !x.crange[1]) then begin
                       
                       ;; Take plot range into account
                       tinterval = [tinterval[0] >!x.crange[0]/3600, tinterval[1] <!x.crange[1]/3600]
                       
                       xpoly = tinterval[[0, 1, 1, 0, 0]]
                       ypoly = [0, 0, 1, 1, 0]*0.01

                       if strmatch(timedirs[i],'*dark*',/fold) then begin
                          ;; Darks
                          cgwindow, 'cgColorFill', /loadcmd, xpoly*3600, (ypoly*2+0.01) $
                                    , color = color_dark
                       endif else if strmatch(timedirs[i],'*flat*', /fold) then begin
                          ;; Flats
                          cgwindow, 'cgColorFill', /loadcmd, xpoly*3600, (ypoly*2+0.01) $
                                    , color = color_flat
                       endif else if strmatch(timedirs[i],'*pinh*', /fold) or $
                          strmatch(timedirs[i], '*[Gg]rid*', /fold) then begin
                          ;; Pinholes or TRIPPEL grid
                          cgwindow, 'cgColorFill', /loadcmd, xpoly*3600, (ypoly*2+0.01) $
                                    , color = color_pinhole
                       endif else if strmatch(timedirs[i],'*polcal*', /fold) then begin
                          ;; Polcal
                          cgwindow, 'cgColorFill', /loadcmd, xpoly*3600, (ypoly*2+0.01) $
                                    , color = color_polcal
                       endif else begin
                          
                          ;; This is not calibration data. So which
                          ;; instrument is it from? Could be more than
                          ;; one! 

                          subdirs = file_search(datedir+timedirs[i]+'*')

                          ;; Blue data
                          if total(strmatch(file_basename(subdirs),'Chromis-?')) gt 0 $
                             or strmatch(timedirs[i],'*blue*', /fold) then begin
                             cgwindow, 'cgColorFill', /loadcmd, xpoly*3600, (ypoly+0.010) $
                                       , color = color_blue
                          endif ; Blue/CHROMIS

                          ;; Red data
                          if total(strmatch(file_basename(subdirs),'Crisp-?')) gt 0 then begin
                             cgwindow, 'cgColorFill', /loadcmd, xpoly*3600, (ypoly+0.020) $
                                       , color = color_red
                          endif ; Red/CRISP

                          ;; TRIPPEL data?
                          if total(strmatch(file_basename(subdirs),'Spec*')) gt 0 then begin
                             cgwindow, 'cgColorFill', /loadcmd, xpoly*3600, (ypoly*2+0.01) $
                                       , color = color_spec
                          endif ; TRIPPEL
                          
                       endelse    ; data type
                    endif         ; within range
                 endif            ; intfile
              endfor              ; i
           endif                  ; Ntd
        endif                     ; Ndd
        
     endif else begin           ; markdata or onlydata
        ;; Set up the plot.
        L = red_gen_timeaxis([tmin, tmax]*3600.)
        cgplot, /add, [tmin, tmax]*3600., [0, r0max], /nodata, /ystyle $
                , xtitle = 't [UT]', ytitle = 'r$\sub0$ / 1 m' $
                , title = isodate $
                , xrange = [tmin, tmax]*3600. $
                , XTICKV=L.tickv, XTICKS=L.ticks, XMIN=L.minor, XTICKNAM=L.name

     endelse

     ;; Return time limits of plot
     tmin = !x.crange[0]/3600.
     tmax = !x.crange[1]/3600.

     ;; Flatline r0 values outside these limits:
     r0_top = r0max * 0.99
     r0_bot = 0.040
     
     ;; Smooth over a dt_mean*60 wide window because r0time has a 1
     ;; sec step.

     ;; Do the r0 plotting
     if keyword_set(plot24) then begin
        indx = where(r0data[0, *] gt r0_top, Ntop)
        cgwindow, 'cgplot', /loadcmd, /over, r0time, r0data[0, *] >r0_bot<r0_top $
                  , color = color_24, psym = psym_plot, symsize = symsize_plot
        if Ntop gt 0 then cgwindow, 'cgplot', /loadcmd, /over $
                                    , r0time[indx], replicate(r0_top, Ntop) $
                                    , color = color_24, psym = 5, symsize = 0.15
        if dt_mean ne 0 then cgwindow, 'cgplot', /loadcmd, /over $
                                       , r0time, smooth(r0data[0, *], dt_mean*60.) >r0_bot $
                                       , color = color_24
     endif

     if keyword_set(plot8) then begin
        indx = where(r0data[1, *] gt r0_top, Ntop)
        cgwindow, 'cgplot', /loadcmd, /over, r0time, r0data[1, *] >r0_bot<r0_top $
                  , color = color_8, psym = psym_plot, symsize = symsize_plot
        if Ntop gt 0 then cgwindow, 'cgplot', /loadcmd, /over $
                                    , r0time[indx], replicate(r0_top, Ntop) $
                                    , color = color_8, psym = 5, symsize = 0.15
        if dt_mean ne 0 then cgwindow, 'cgplot', /loadcmd, /over $
                                       , r0time, smooth(r0data[1, *], dt_mean*60.) >r0_bot $
                                       , color = color_8
     endif

     pname = 'dir-analysis/' + 'r0_' + isodate + '_' $
             + string(long(tmin),format='(i02)')+':'+string(60*(tmin-long(tmin)),format='(i02)') $
             + '-' $
             + string(long(tmax),format='(i02)')+':'+string(60*(tmax-long(tmax)),format='(i02)') $
             + '.pdf'

     cgcontrol, output = pname
     
  endif                         ; plot24 or plot8

  if keyword_set(scan24) or keyword_set(scan8) then begin
     
     ;; Plot r0 for all scans.
     
     cgwindow

     for iinstrument = 0, Ninstruments-1 do begin
        
        print,'red_plot_r0 : Plot r0 '+instruments[iinstrument]+' scans.'

        scanfiles = file_search(analysis_dir+'/*/*/'+instruments[iinstrument]+'-scans.txt', count = Nscanfiles, /fold)
        
        for iscanfile = 0, Nscanfiles-1 do begin

           print,strtrim(iscanfile,2)+'/'+strtrim(Nscanfiles,2)
           
           print, scanfiles[iscanfile]
           
           instrument = file_basename(scanfiles[iscanfile],'-scans.txt')
           timedir = file_dirname(scanfiles[iscanfile])
           
           if ~query_ascii(scanfiles[iscanfile]) then continue ; Skip non-valid files
           scaninfo = read_ascii(scanfiles[iscanfile], template = scantemplate)

           sindx = sort(scaninfo.tstart)
           Ns = n_elements(sindx) ; Total number of scans of all types

           upref = scaninfo.prefilter(uniq(scaninfo.prefilter, sort(scaninfo.prefilter)))
           Npref = n_elements(upref)
           print, upref
           
           Nscans = max(scaninfo.scanno)+1
           
           interval_lengths = fltarr(Npref)
           for ipref = 0L, Npref-1 do begin
              indx = where(scaninfo.prefilter eq upref[ipref])
              interval_lengths[ipref] = median((scaninfo.tstop-scaninfo.tstart)[indx])*3600
           endfor               ; ipref

           ;; We base number of plots on the shortest scans, which
           ;; should be no less than 1/15 of the plot time interval.
           tfrac = 20
           Nplots = ceil(3600.*(max(scaninfo.tstop)-min(scaninfo.tstart))/(tfrac*min(interval_lengths)))
           Tplot = (max(scaninfo.tstop)-min(scaninfo.tstart))/Nplots
           Tstart = min(scaninfo.tstart) - Tplot
           Tstop = Tstart + Tplot + max(scaninfo.tstop-scaninfo.tstart)

           iplot = 0

           for is = 0, Ns-1 do begin
              
              print, is, Ns, format = '(i0, "/", i0)'
              
              r0scanfile = timedir $
                           + '/r0data_scan' + strtrim(scaninfo.scanno[sindx[is]], 2) + '_pre' $
                           + strtrim(scaninfo.prefilter[sindx[is]], 2) $
                           + '.txt'

              if file_test(r0scanfile) then begin
                 
                 if ~query_ascii(scanfiles[iscanfile]) then continue ; Skip non-valid files
                 r0info = read_ascii(r0scanfile, template = r0scantemplate)
                 
                 sstart = min(r0info.time) ; Start and stop of this particular scan
                 sstop  = max(r0info.time)
                 nofile = 0
                 
              endif else nofile = 1
              
              
              if sstop gt Tstop*3600 or is eq 0 or nofile then begin

                 ;; Time to start a new plot.
                 if is ne 0 then begin
                    iplot += 1
                    pname = 'dir-analysis/' + 'r0plot_' + isodate + '_' $
                            + (strsplit(timedir,'/',/extr,count=Nsplit))[Nsplit-1] $
                            + '_' + string(iplot, format = '(i04)') + '.' + extension
                    cgcontrol, output = pname
                    
                 endif
                 
                 cgcontrol, /delete, /all
                 Tstart += Tplot
                 Tstop += Tplot
                 if keyword_set(scan24) then title = 'r0 24x24 : ' else title = 'r0 8x8 : '
                 title += instrument + ' ' + isodate + ' ' + file_basename(timedir)
                 L = red_gen_timeaxis([Tstart, Tstop]*3600)
                 cgwindow, /add, 'cgplot', /nodata, [0], [0] $
                           , ticklen = -!p.ticklen $
                           , xrange = [Tstart, Tstop]*3600 $
                           , yrange = [0, r0max] $
                           , xtitle = 't [UT]', ytitle = 'r$\sub0$ / 1 m' $
                           , title = title $
                           , XTICKV=L.tickv, XTICKS=L.ticks, XMIN=L.minor, XTICKNAM=L.name
              endif
              
              cgwindow, /load, 'cgColorFill' $
                        , [sstart, sstop, sstop, sstart, sstart] >!x.crange[0] <!x.crange[1] $
                        , !y.crange([0, 0, 1, 1, 0]) $
                        , color = precolors[where(upref eq scaninfo.prefilter[sindx[is]])]
              
              if keyword_set(plotvertical) or Npref eq 1 then begin
                 cgwindow, /load, 'cgplot', /over, [1, 1]*sstart >!x.crange[0] <!x.crange[1], !y.crange([0, 1]) 
                 cgwindow, /load, 'cgplot', /over, [1, 1]*sstop  >!x.crange[0] <!x.crange[1], !y.crange([0, 1]) 
              endif
              
              tpos = (sstart+sstop)/2.
              if tpos gt !x.crange[0] and tpos lt !x.crange[1] then begin
                 cgwindow, /load, 'cgtext', tpos, !y.crange[1]*.95 $
                           , strtrim(string(scaninfo.scanno[sindx[is]]), 2) $
                           , align = 0.5 
                 if scaninfo.scanno[sindx[is]] eq 0 then $
                    cgwindow, /load, 'cgtext', tpos, !y.crange[1]*.03 $
                              , strtrim(scaninfo.prefilter[sindx[is]], 2) $
                              , align = 0, charsize = .75, orientation = 90 
              endif
              
              if keyword_set(scan24) then begin
                 cgwindow, /load, 'cgplot', /over, r0info.time, r0info.r0_24x24 $
                           , psym = psym_scan, symsize = symsize_scan
                 cgwindow, /load, 'cgplot', /over, color = 'blue' $
                           , [sstart, sstop], [1, 1]*scaninfo.r0_24x24_mean[sindx[is]] 
                 cgwindow, /load, 'cgplot', /over, color = 'red' $
                           , [sstart, sstop], [1, 1]*scaninfo.r0_24x24_median[sindx[is]] 
              endif
              
              if keyword_set(scan8) then begin
                 cgwindow, /load, 'cgplot', /over, r0info.time, r0info.r0_8x8 $
                           , psym = psym_scan, symsize = symsize_scan
                 cgwindow, /load, 'cgplot', /over, color = 'blue' $
                           , [sstart, sstop], [1, 1]*scaninfo.r0_8x8_mean[sindx[is]] 
                 
                 cgwindow, /load, 'cgplot', /over, color = 'red' $
                           , [sstart, sstop], [1, 1]*scaninfo.r0_8x8_median[sindx[is]] 
              endif
              
              
           endfor               ; is
           
        endfor                  ; iscanfile
     endfor                     ; iinstrument

  endif                         ; plotscans
  
  
  if keyword_set(plotstats) then begin
     ;; Now plot for each timedir only the r0 statistics vs scan
     ;; number. 

     cgwindow

     for iinstrument = 0, Ninstruments-1 do begin
        
        print,'red_plot_r0 : Plot r0 statistics vs scan number for '+instruments[iinstrument]+'.'

        scanfiles = file_search(analysis_dir+'/*/*/'+instruments[iinstrument]+'-scans.txt', count = Nscanfiles, /fold)

        for iscanfile = 0, Nscanfiles-1 do begin

           print,strtrim(iscanfile,2)+'/'+strtrim(Nscanfiles,2)
           
           print, scanfiles[iscanfile]
           
           instrument = file_basename(scanfiles[iscanfile],'-scans.txt')
           timedir = file_dirname(scanfiles[iscanfile])
           print, timedir

           if ~query_ascii(scanfiles[iscanfile]) then continue ; Skip non-valid files
           scaninfo = read_ascii(scanfiles[iscanfile], template = scantemplate)

           upref = scaninfo.prefilter(uniq(scaninfo.prefilter, sort(scaninfo.prefilter)))
           Npref = n_elements(upref)
           print, upref
           
           for ipref = 0, Npref-1 do begin

              pindx = where(scaninfo.prefilter eq upref[ipref], Ns)

              cgcontrol, /delete, /all
              title = instrument + ' ' + upref[ipref] + ' ' + isodate + ' ' + file_basename(timedir)
              
              cgwindow, /add, 'cgplot', /nodata, [0], [0] $
                        , xrange = [min(scaninfo.scanno), max(scaninfo.scanno)] $
                        , yrange = [0, r0max] $
                        , xtitle = 'scan #', ytitle = 'r$\sub0$ / 1 m' $
                                ;, title = title $
                        , xstyle = 8 

              ;; Add a time axis on top
              tmin = min(scaninfo.tstart)
              tmax = max(scaninfo.tstop) 
              L = red_gen_timeaxis([tmin, tmax]*3600.)

              cgwindow, /add, 'cgaxis', xaxis = 1 $ ; , xrange = [min(scaninfo.scanno), max(scaninfo.scanno)] $
                                ;, xtitle = 'scan #' $
                        , /xstyle $
                        , xrange = [tmin, tmax]*3600. $
                        , XTICKV=L.tickv, XTICKS=L.ticks, XMIN=L.minor, XTICKNAM=L.name

              cgtext, mean(!x.crange), !y.crange[1]*.03, title, /data, align = 0.5,/add

              statscolors = ['green', 'red', 'blue', 'cyan']
              statslines = [0, 1]

              xpoly = scaninfo.scanno[[pindx, reverse(pindx), 0]]
              ypoly = [scaninfo.r0_24x24_min[pindx], scaninfo.r0_24x24_max[reverse(pindx)], scaninfo.r0_24x24_min[0]]
              cgwindow, /add, 'cgcolorfill', xpoly, ypoly, color = 'Light Gray'

              xpoly = scaninfo.scanno[[pindx, reverse(pindx), 0]]
              ypoly = [scaninfo.r0_8x8_min[pindx], scaninfo.r0_8x8_max[reverse(pindx)], scaninfo.r0_8x8_min[0]]
              cgwindow, /add, 'cgcolorfill', xpoly, ypoly, color = 'Pale Goldenrod'

;           cgwindow, /add, /over, 'cgplot', scaninfo.scanno[pindx], scaninfo.r0_8x8_min[pindx] $
;                     , linestyle = statslines[0], color = statscolors[3]
;           cgwindow, /add, /over, 'cgplot', scaninfo.scanno[pindx], scaninfo.r0_8x8_max[pindx] $
;                     , linestyle = statslines[0], color = statscolors[0]
              cgwindow, /add, /over, 'cgplot', scaninfo.scanno[pindx], scaninfo.r0_8x8_mean[pindx] $
                        , linestyle = statslines[0], color = statscolors[1]
              cgwindow, /add, /over, 'cgplot', scaninfo.scanno[pindx], scaninfo.r0_8x8_median[pindx] $
                        , linestyle = statslines[0], color = statscolors[2]
              
;           cgwindow, /add, /over, 'cgplot', scaninfo.scanno[pindx], scaninfo.r0_24x24_min[pindx] $
;                     , linestyle = statslines[1], color = statscolors[3]
;           cgwindow, /add, /over, 'cgplot', scaninfo.scanno[pindx], scaninfo.r0_24x24_max[pindx] $
;                     , linestyle = statslines[1], color = statscolors[0]
              cgwindow, /add, /over, 'cgplot', scaninfo.scanno[pindx], scaninfo.r0_24x24_mean[pindx] $
                        , linestyle = statslines[1], color = statscolors[1]
              cgwindow, /add, /over, 'cgplot', scaninfo.scanno[pindx], scaninfo.r0_24x24_median[pindx] $
                        , linestyle = statslines[1], color = statscolors[2]

;           cglegend, /add, align = 0, location = [.15, .85], linestyle = 0 $
;                     , title = ['max', 'mean', 'median', 'min'], colors = statscolors
              cglegend, /add, align = 0, location = [.15, .85], linestyle = 0 $
                        , title = ['mean', 'median'], colors = statscolors[[1, 2]]

              pname = 'dir-analysis/' + 'r0stats_'+ upref[ipref] + '_' + isodate + '_' $
                      + (strsplit(timedir,'/',/extr,count=Nsplit))[Nsplit-1] $
                      + '.pdf'
              cgcontrol, output = pname

           endfor               ; ipref
        endfor                  ; iscanfile
     endfor                     ; iinstrument

  endif

end
