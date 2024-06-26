; docformat = 'rst'

;+
; Sets up a work directory for the TRIPPEL instrument.
;
; :Categories:
;
;   SST pipeline
;
; :Author:
;
;   Mats Löfdahl, Institute for Solar Physics
;
; :Params:
;
;   work_dir : in, type = string
;
;     A name of the directory to set up the directory tree at and to
;     where the config file and the doit.pro script should be put.
;
;   root_dir : in, type = string
;
;     A name of the directory where the raw data is stored for a date
;     given in isodate.
;
;   cfgfile : in, type = string
;
;     A name of the config file, the default value is 'config.txt'.
;
;   scriptfile : in, type = string
;
;     A name of the set-up script, the default value is 'doit.pro'.
;
;   isodate : in, type = string
;
;     A date given in ISO format (e.g., 2017-08-23) to search raw data
;     for in root_dir.
;
; :Keywords:
;
;   calibrations_only : in, optional, type = boolean,
;
;      Set up to process calibration data only. Currently, a dummy
;      keyword that does nothing.
; 
;   no_observer_metadata : in, optional, type=boolean
;
;      Do not offer to add OBSERVER metadata.
;
;   old_dir : in, optional, type = string
;
;      Copy files from this directory, in particular summed
;      calibration data. Useful if you summed the calibration data in
;      La Palma.
;
; :History:
;
;   2017-08-21 : AVS. A calibrations_only keyword is added. An rst
;                format header is added as well.
;
;   2019-07-01 : MGL. New keyword old_dir.
; 
;   2021-01-25 : MGL. New keyword no_observer_metadata.
;
;-
pro red_setupworkdir_trippel, work_dir, root_dir, cfgfile, scriptfile, isodate $
                              , calibrations_only = calibrations_only $
                              , no_observer_metadata = no_observer_metadata $
                              , old_dir = old_dir 

  inam = red_subprogram(/low, calling = inam1)

  file_mkdir, work_dir

  red_metadata_store, fname = work_dir + '/info/metadata.fits' $
                      , [{keyword:'INSTRUME', value:'TRIPPEL' $
                          , comment:'Name of instrument'} $
                         , {keyword:'TELCONFG', value:'Schupmann, spectrograph table', $
                            comment:'Telescope configuration'}]

  ;; Offer to add the OBSERVER metadata keyword
  if ~keyword_set(no_observer_metadata) then begin
    observer = ''
    read, 'Add names for OBSERVER metadata keyword for the TRIPPEL workdir or hit return:', observer
    ;; Write it to the metadata file
    if observer ne '' then begin
      print
      print, inam+' : Adding to CHROMIS metadata, OBSERVER = '+observer
      red_metadata_store, fname = work_dir + '/info/metadata.fits' $
                          , [{keyword:'OBSERVER', value:observer $
                              , comment:'Observer name(s)'}]
    endif
  endif
  
  ;; Open two files for writing. Use logical unit Clun for a Config
  ;; file and Slun for a Script file.

  openw, Clun, work_dir + cfgfile,    /get_lun
  openw, Slun, work_dir + scriptfile, /get_lun

  ;; Specify the date in the config file, ISO format.
  print, 'Date'
  printf, Clun, '#'
  printf, Clun, '# --- Settings'
  printf, Clun, '#'
  printf, Clun,'isodate = '+isodate

  printf, Slun, 'a = trippelred("'+cfgfile+'")'
  printf, Slun, 'root_dir = "' + root_dir + '"'

  ;; Download SST log files and optionally some other data from the web.
  print, 'Log files'
  printf, Clun, '#'
  printf, Clun, '# --- Download SST log files'
  printf, Clun, '#'
  printf, Slun
  printf, Slun, 'a -> download ; add ", /all" to get also HMI images and AR maps.'

  print, 'Cameras'
  printf, Clun, '#'
  printf, Clun, '# --- Cameras'
  printf, Clun, '#'
  printf, Clun, 'camera = SpecC'
  printf, Clun, '#'
  printf, Clun, 'root_dir = ' + root_dir
  printf, Clun, '#'

  print, 'Output'
  printf, Clun, '#'
  printf, Clun, '# --- Output'
  printf, Clun, '#'
  printf, Clun, 'out_dir = ' + work_dir

  print, 'Darks'
  printf, Clun, '#'
  printf, Clun, '# --- Darks'
  printf, Clun, '#'

  darksubdirs = red_find_instrumentdirs(root_dir, 'spec', '*dark*' $
                                        , count = Nsubdirs)
  if Nsubdirs gt 0 then begin
    darkdirs = file_dirname(darksubdirs)
    darkdirs = darkdirs[uniq(darkdirs, sort(darkdirs))]
    printf, Slun
    for idir = 0, n_elements(darkdirs)-1 do begin
      printf, Clun, 'dark_dir = '+red_strreplace(darkdirs[idir] $
                                                 , root_dir, '')
      printf, Slun, 'a -> sumdark, /check, dirs=root_dir+"' $
              + red_strreplace(darkdirs[idir], root_dir, '') + '"'
    endfor                      ; idir
  endif                         ; Nsubdirs

  print, 'Flats'
  printf, Clun, '#'
  printf, Clun, '# --- Flats'
  printf, Clun, '#'

  flatsubdirs = red_find_instrumentdirs(root_dir, 'spec', '*flat*' $
                                        , count = Nsubdirs)
  if Nsubdirs gt 0 then begin
    ;; There are spectral flats!

    ;; Directories with camera dirs below:
    flatdirs = file_dirname(flatsubdirs)
    flatdirs = flatdirs[uniq(flatdirs, sort(flatdirs))]
    Nflatdirs = n_elements(flatdirs)

    ;; Loop over the flatdirs, write each to the config file and
    ;; to the script file
    printf, Slun
    for idir = 0, Nflatdirs-1 do begin

      ;; Config file

      printf, Clun, 'flat_dir = '+red_strreplace(flatdirs[idir] $
                                                 , root_dir, '')

      ;; Script file

      ;; Look for wavelengths in those flatsubdirs that match
      ;; flatdirs[idir]! Also collect prefilters.
      indx = where(strmatch(flatsubdirs, flatdirs[idir]+'*'))
      fnames = file_search(flatsubdirs[indx]+'/cam*', count = Nfiles)

      if Nfiles gt 0 then begin

        camdirs = strjoin(file_basename(flatsubdirs[indx]), ' ')

        red_extractstates, fnames, /basename, pref = wls
        wls = wls[uniq(wls, sort(wls))]
        wls = wls[WHERE(wls ne '')]
        wavelengths = strjoin(wls, ' ')
        ;; Print to script file
        printf, Slun, 'a -> sumflat, /check, dirs=root_dir+"' $
                + red_strreplace(flatdirs[idir], root_dir, '') $
                + '"  ; ' + camdirs+' ('+wavelengths+')'

        red_append, prefilters, wls

      endif                     ; Nfiles
    endfor                      ; idir
  endif                         ; Nsubdirs

  if n_elements(prefilters) gt 0 then $
     prefilters = prefilters[uniq(prefilters, sort(prefilters))]
  Nprefilters = n_elements(prefilters)

  if Nprefilters eq 0 then begin
    ;; This can happen if flats were already summed in La Palma. Look
    ;; for prefilters in the summed flats directory instead.
    spawn, 'ls flats/cam*.flat | cut -d. -f2|sort|uniq', prefilters
    Nprefilters = n_elements(prefilters)
  endif

  printf, Slun, "a -> makegains, /preserve"

  print, 'Pinholes'
  printf, Clun, '#'
  printf, Clun, '# --- Pinholes'
  printf, Clun, '#'
  pinhdirs = file_search(root_dir+'/*grid*/*', count = Ndirs, /fold)
  for i = 0, Ndirs-1 do begin
    pinhsubdirs = file_search(pinhdirs[i]+'/spec*', count = Nsubdirs, /fold)
    if Nsubdirs gt 0 then begin
      printf, Slun
      printf, Clun, 'pinh_dir = '+red_strreplace(pinhdirs[i], root_dir, '')
      printf, Slun, 'a -> setpinhdir, root_dir+"' $
              + red_strreplace(pinhdirs[i], root_dir, '')+'"'
;        printf, Slun, 'a -> sumpinh_new'
      for ipref = 0, Nprefilters-1 do begin
        printf, Slun, "a -> sumpinh, /pinhole_align, pref='"+prefilters[ipref]+"'"
      endfor                    ; ipref
    endif else begin
      pinhsubdirs = file_search(pinhdirs[i]+'/*', count = Nsubdirs)
      printf, Slun
      for j = 0, Nsubdirs-1 do begin
        pinhsubsubdirs = file_search(pinhsubdirs[j]+'/spec*' $
                                     , count = Nsubsubdirs, /fold)
        if Nsubsubdirs gt 0 then begin
          printf, Clun, 'pinh_dir = ' $
                  + red_strreplace(pinhsubdirs[j], root_dir, '')
          printf, Slun, 'a -> setpinhdir, root_dir+"' $
                  + red_strreplace(pinhsubdirs[j], root_dir, '')+'"'
;              printf, Slun, 'a -> sumpinh_new'
          for ipref = 0, Nprefilters-1 do begin
            printf, Slun, "a -> sumpinh, /pinhole_align, pref='" $
                    + prefilters[ipref]+"'"
          endfor                ; ipref
        endif
      endfor                    ; j
    endelse
  endfor                        ; i


  print, 'Science'
  printf, Clun, '#'
  printf, Clun, '# --- Science data'
  printf, Clun, '# '

  ;;  sciencedirs = file_search(root_dir+'/sci*/*', count = Ndirs, /fold)
  red_append, nonsciencedirs, pinhdirs
  red_append, nonsciencedirs, darkdirs
  red_append, nonsciencedirs, flatdirs
  sciencedirs = file_search(root_dir+'/*/*', count = Ndirs)

  for i = 0, Ndirs-1 do begin

    if total(sciencedirs[i] eq nonsciencedirs) eq 0 then begin
      sciencesubdirs = file_search(sciencedirs[i]+'/spec*' $
                                   , count = Nsubdirs, /fold)
      if Nsubdirs gt 0 then begin
        red_append, dirarr, red_strreplace(sciencedirs[i], root_dir, '')
      endif else begin
        sciencesubdirs = file_search(sciencedirs[i]+'/*', count = Nsubdirs)
        for j = 0, Nsubdirs-1 do begin
          sciencesubsubdirs = file_search(sciencesubdirs[j]+'/spec*' $
                                          , count = Nsubsubdirs, /fold)
          if Nsubsubdirs gt 0 then begin
            red_append, dirarr, red_strreplace(sciencesubdirs[j], root_dir, '')
          endif
        endfor                  ; j
      endelse
    endif
  endfor
  if n_elements(dirarr) gt 0 then printf, Clun, "data_dir = ['"+strjoin(dirarr, "','")+"']"



  free_lun, Clun
  free_lun, Slun


end
