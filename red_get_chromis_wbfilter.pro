function red_get_chromis_wbfilter, nbinfo, date, workdir = workdir

  if n_elements(workdir) eq 0 then workdir = './'
  
  linedefdir = workdir + '/' + 'downloads/'
  ldfiles = file_search(linedefdir + '/linedef.*', count = nLD)
  ;; Find latest file that is before date, for now use first file
  ild = 0
  
  linedef = red_readlinedef(ldfiles[ild])

  ;; Get the nb wheel number from nbinfo
  case 1 of
    
    isa(nbinfo, /number) : begin
      ;; A number, could be the wheel number or the wavelength in Ångström
      if nbinfo lt 100 then begin
        ;; Could imagine >10 filters in a wheel, not >100.
        nbwheel = long(nbinfo)
      endif else begin
        ;; Filter wavelength
        nbfilt = string(nbinfo, format = '(i04)')
      endelse
    end

    isa(nbinfo, /string) : begin
      ;; Could be the four-digit filter tag, the wheel number
      if strmatch(nbinfo,'[0-9][0-9][0-9][0-9]') then begin
        ;; Filter tag
        nbfilt = nbinfo
      endif else begin
        ;; Wheel number
        nbwheel = long(red_strreplace(nbinfo, 'wheel', ''))        
      endelse
    end

  endcase

  if n_elements(nbfilt) gt 0 then begin
    ;;nbwheel = linedef.wb_map[where(linedef.nb_wl eq nbfilt, Nbmatch)]
    nbwheel = where(linedef.nb_wl eq nbfilt, Nbmatch)
    if Nbmatch eq 0 then stop
  endif

  wbwheel = linedef.wb_map[nbwheel]
  
  wbfilt = strtrim(linedef.wb_wl[wbwheel], 2)
  
  return, wbfilt
  
end

print, red_get_chromis_wbfilter('wheel00005', workdir = '/scratch/mats/2019-11-11/CHROMIS/')
print, red_get_chromis_wbfilter('wheel00005', workdir = '/scratch/mats/2020-08-18/CHROMIS/')
print, red_get_chromis_wbfilter('wheel00004', workdir = '/scratch/mats/2020-08-18/CHROMIS/')
print, red_get_chromis_wbfilter('3969', workdir = '/scratch/mats/2020-08-18/CHROMIS/')
print, red_get_chromis_wbfilter('3934', workdir = '/scratch/mats/2020-08-18/CHROMIS/')
print, red_get_chromis_wbfilter('4862', workdir = '/scratch/mats/2020-08-18/CHROMIS/')

end
