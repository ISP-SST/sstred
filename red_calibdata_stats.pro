; docformat = 'rst'

;+
; Calculate statistics on automatically collected calibration data.
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
;   dates : in, optional, type=strarr, default="All dates"
;
;      Dates for which to calculate statistics.
; 
;   force : in, optional, type=boolean
;   
;      Remove any existing data and calculate new stats.
; 
; 
; :History:
; 
;   2017-10-12 : MGL. First version.
; 
; 
; 
; 
;-
pro red_calibdata_stats, dir = dir, force = force, dates = dates

  if n_elements(dir) eq 0 then dir = '/freija-scratch/asukh/calibration-data/'
  if ~strmatch( dir, '*/' ) then dir += '/'

  instruments = ['CHROMIS', 'CRISP'] 
  data_types = ['darks', 'flats', 'pinhs', 'polcal_sums']

  Ninstr = n_elements(instruments)
  Ntypes = n_elements(data_types)

  if n_elements(date) eq 0 then begin
    dates = file_basename(file_search(dir+'/????.??.??', count = Ndates))
  endif else begin
    Ndates = n_elements(dates)
  endelse
  
  for idate = 0, Ndates-1 do begin

    ;; Remove old statistics if force
    if keyword_set(force) then file_delete, dates[idate], /recursive, /allow_nonexistent

    ;; Avoid recalculating existing statistics
    if ~file_test(dates[idate]) then begin

      print, dates[idate]
      
      ;; Loop over instruments and data types
      for iinstr = 0, Ninstr-1 do begin
        for itype = 0, Ntypes-1 do begin

          timestamps = file_basename(file_search(dir + dates[idate]+'/' $
                                                 + instruments[iinstr]+'-calibrations/' $
                                                 + data_types[itype]+'/' $
                                                 + '??:??:??' $
                                                 , count = Ntimes ))

          if Ntimes gt 0 then begin

            tdir_in = dir + dates[idate]+'/' $
                      + instruments[iinstr]+'-calibrations/' $
                      + data_types[itype]+'/'

            tdir_out = dates[idate]+'/' $
                       + instruments[iinstr]+'/' $
                       + data_types[itype]+'/'

            file_mkdir, tdir_out

            print, tdir_in

            for itime = 0, Ntimes-1 do begin

              case data_types[itype] of
                'darks'       : fregex = 'cam*dark.fits'
                'flats'       : fregex = 'cam*_summed.flat.fits'
                'pinhs'       : fregex = 'cam*pinh'
                'polcal_sums' : fregex = 'Crisp-'+['R', 'T']+'/cam*.fits'
              endcase

              fnames = file_search(tdir_in $
                                   + timestamps[itime]+'/' $
                                   + fregex $
                                   , count = Nfiles )

              if Nfiles gt 0 then begin

                ;; Open a file in which to write statistics for files
                ;; with a certain timestamp, one line per file.
                openw, lun, /get_lun, tdir_out+'stats_'+timestamps[itime]+'.txt'

                printf, lun, '# filename min max median moments (mean, variance, skewness, kurtosis) '
                
                for ifile = 0, Nfiles-1 do begin

                  print, fnames[ifile]
                  statsline = fnames[ifile] + ' '
                  data = red_readdata(fnames[ifile], /silent)

                  if size(data,/n_dim) eq 0 && data_types[itype] eq 'flats' then begin

                    ;; Kludge to reconstruct flats, only ordinary
                    ;; flats were stored for many data sets and not
                    ;; the summed, not-dark-corrected versions that we
                    ;; want to use.
                    
                    h = headfits(fnames[ifile])
                    Nsumexp = fxpar(h, 'NSUMEXP')
                    detector = strtrim(fxpar(h, 'DETECTOR'), 2)

                    ;; Read the ordinary flat and undo dark
                    ;; correction with the actual dark that was
                    ;; used. 
                    fn = red_strreplace(fnames[ifile], '_summed', '')
                    data = red_readdata(fn, /silent)
                    dn = file_dirname(tdir_in) + '/darks/'+detector+'.dark.fits'
                    dark = red_readdata(dn, /silent)

                    data = long((data + dark)*Nsumexp)
                    
                  endif

                  if size(data,/n_dim) eq 0 then begin
                    ;; This was not a flat that could be reconstructed!
                    statsline += 'No data'
                  endif else begin
                    
                    statsline += string(min(data)) + ' '
                    statsline += string(max(data)) + ' '
                    statsline += string(median(data)) + ' '

                    m = moment(data)
                    statsline += string(m[0]) + ' '
                    statsline += string(m[1]) + ' '
                    statsline += string(m[2]) + ' '
                    statsline += string(m[3]) + ' '
                  endelse
                  
                  printf, lun, statsline
                  
                endfor          ; ifile

                ;; Close the file
                free_lun, lun
                
              endif
              
            endfor              ; itime
            
          endif
          
        endfor                  ; itype
      endfor                    ; iinstr
    endif
stop
  endfor                        ; idate
  
end
