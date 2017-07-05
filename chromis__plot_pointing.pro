pro chromis::plot_pointing, file = file

  if n_elements(file) gt 0 then cgPS_Open, file, /decomposed

  psym_calib = 16
  psym_data = 16

  symsize_calib = 0.3
  symsize_data = 1
  
  ;; Get pointing data
  red_logdata, self.isodate, time_pig, pig = metadata_pig, rsun = rsun

  red_plot_diskcoordinates, metadata_pig, rsun=rsun, rplot=1.01, color = 'gray' $
                            , title = 'SST/CHROMIS '+self.isodate + ' pointing'

  dirs = file_search('dir-analysis/'+self.isodate+'/CHROMIS-*', count = Ndirs)

  for idir = 0, Ndirs-1 do begin

    subdirs = file_search(dirs[idir]+'/??:??:??', count = Nsubdirs)
    
    if Nsubdirs gt 0 then begin

      if file_basename(dirs[idir]) eq 'CHROMIS-data' then begin

        ;; Find clumping of pointings
        for isubdir = 0, Nsubdirs-1 do begin
          openr, lun, /get_lun, subdirs[isubdir]+'/interval.txt'
          readf, lun, tbeg
          readf, lun, tend
          free_lun, lun
          indx = where(time_pig ge tbeg*3600 and time_pig le tend*3600, count)
          h2 = hist_2d(metadata_pig[0, indx]/Rsun,metadata_pig[1, indx]/Rsun $
                       ,bin1=.01,bin2=.01,min1=-1.1,min2=-1.1,max1=1.1,max2=1.1)
          if isubdir eq 0 then h2tot = h2 else h2tot += h2
        endfor                  ; isubdir
        dims = size(h2tot, /dim)
        bim = bytarr(dims+2)
        bim[1, 1] = label_region(h2tot gt 0) ; Big image
        lim = bim[1:dims[0], 1:dims[1]]      ; Label image
        Nclumps = max(lim)

        if Nclumps gt 8 then begin
          print, 'Too many disjoint pointing regions'
          stop
          retall
        endif
        
        rgb = distinct_colors(N_COLORS = Nclumps)
        r = bytarr(dims)+255b
        g = bytarr(dims)+255b
        b = bytarr(dims)+255b
        for iclump = 0, Nclumps-1 do begin
          indx = where(lim eq iclump+1, count)
          if count gt 0 then begin
            r[indx] = rgb[0, iclump]
            g[indx] = rgb[1, iclump]
            b[indx] = rgb[2, iclump]
          endif
        endfor                  ; iclump
        rgbim = [[[r]],[[g]],[[b]]]
        cgimage, rgbim, /keep, /axes, xrange = [-1.1, 1.1], yrange = [-1.1, 1.1]


stop
        
        hours = strmid(file_basename(subdirs),0,2)
        uhours = hours(uniq(hours))
        
        ucolors = distinct_colors(N_COLORS = n_elements(uhours), /num)
        colors = lonarr(Nsubdirs)
        for isubdir = 0, Nsubdirs-1 do colors[isubdir] = ucolors[where(uhours eq hours[isubdir])]
        
        psym = psym_data
        symsize = symsize_data

        upsym = replicate(psym, n_elements(uhours))
        usymsize = replicate(symsize, n_elements(uhours))

      endif else begin
        colors = replicate(0, Nsubdirs)
        psym = psym_calib
        symsize = symsize_calib
      endelse
      
      for isubdir = 0, Nsubdirs-1 do begin

        openr, lun, /get_lun, subdirs[isubdir]+'/interval.txt'
        readf, lun, tbeg
        readf, lun, tend
        free_lun, lun

        indx = where(time_pig ge tbeg*3600 and time_pig le tend*3600, count)

        if count gt 0 then $
           red_plot_diskcoordinates, /over, metadata_pig, indx = indx, rsun=rsun, rplot=1.01 $
                                     , color = colors[isubdir], psym = psym, symsize = symsize

      endfor                    ; isubdir
    endif
    
  endfor                        ; idir

  titles =  ['Data '+string(uhours, format = '(i02)')+'-'+string(uhours+1, format = '(i02)') $
             , 'Calibrations']
  
  cglegend, titles = titles, colors = [ucolors, 0] $
            , align = 5, location = [0, -1], /data $
            , psym = [upsym, psym_calib] $
            , symsize = symsize_data, length = 0.0

  ;; cglegend does not allow symsize to be an array :(
  
  if n_elements(file) gt 0 then cgPS_Close;, /PDF, /Delete_PS

end
