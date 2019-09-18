; docformat = 'rst'

;+
; Make links to raw data in the form that is required for momfbd
; processing with aid of crisp__extractstates_db 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;    Oleksii Andriienko, Institute for Solar Physics
; 
; :Params:
; 
; 
; :Keywords:
;
;    dirs: in, optional, type = strarr
;
;       Names of directories (or just timestamps) with raw data
; 
;    link_dir : in, optional, type = string
;
;       A name of directory for links
;   
;    uscan : in, optional, type = string
;   
;       Only process scan 'scan'
;   
;    pref : in, optional, type=string
;
;       Only process prefilter 'pref'
; 
; :History:
; 
;   2019-07-10 : OA. Derived from crisp__link_data
;
;   2019-07-25 :  OA. Added check for sst_db.
;  
;-
pro crisp::link_data_db, dirs = dirs $
                      , pref = pref $
                      , uscan = uscan $
                      , link_dir = link_dir 

  if n_elements(link_dir) eq 0 then link_dir = 'data'
  if n_elements(uscan) eq 0 then uscan = ''

  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)

  dataset = 0B
  Ndirs=n_elements(dirs)
  if Ndirs eq 0 then begin
    ;; No dirs given in keyword, use default
    dirs = *self.data_dirs
  endif else begin
    for idir = 0, Ndirs-1 do begin
      if ~file_test(dirs[idir], /dir) then begin
        ;; The keyword doesn't point to an existing directory.
        ;; Try to interpret as a selection from the default dirs.
        dirs[idir] = file_dirname((*self.data_dirs)[0])+'/'+dirs[idir]
      endif
    endfor                                 ;idir
  endelse

  Ndirs = n_elements(dirs)
  if( Ndirs eq 0) then begin
    print, inam+' : ERROR : no directories defined'
    return
  endif 

  linkerdir = self.out_dir + '/' + 'link_scripts' + '/'
  file_mkdir, linkerdir

  cams = *self.cameras
  Ncams = n_elements(cams)

  ;; Create file list
  for idir = 0L, Ndirs - 1 do begin
    data_dir = dirs[idir]
    print, inam + ' : Folder -> ' + data_dir
    folder_tag = strsplit(data_dir,'/',/extract)
    nn = n_elements(folder_tag) - 1
    folder_tag = folder_tag[nn]

    if self.db_present then begin
      dat = stregex(data_dir, '20[0-2][0-9][.-][01][0-9][.-][0-3][0-9]', /extract)
      Ndots   = n_elements(strsplit(dat, '.', /extract))
      Ndashes = n_elements(strsplit(dat, '-', /extract))
      if Ndots gt Ndashes then dat = red_strreplace(dat, '.', '-', n = 2)    
      timestamp = stregex(data_dir, '[0-2][0-9]:[0-5][0-9]:[0-6][0-9]', /extract)
      dataset = dat + ' '  + timestamp
      self->extractstates_db,dataset,sts
    endif

    for icam = 0L, Ncams-1 do begin

      cam = cams[icam]
      detector = self->getdetector( cam )
      iswb = strmatch(cam,'*-[DW]')

      inn = where(sts.camera eq cam)
      if n_elements(inn) eq 1 then if inn eq -1 then continue ; skipping the camera
      states = sts(inn)

      files = states.filename      

      if self.db_present then begin
        inn = where(sts.camera eq cam)
        if n_elements(inn) eq 1 then if inn eq -1 then continue ; skipping the camera
        states = sts(inn)
        files = states.filename 
      endif else begin
        files = file_search(data_dir + '/' + cam + '/*cam*', count = Nfiles)
        files = files(where(strpos(files, '.lcd.') LT 0, nf))
        if(files[0] eq '') then begin
          print, inam+' : ERROR : '+cam+': no files found in: '+$
                 data_dir +' : skipping camera!'
          continue
          ;; Only one prefilter?
          if keyword_set(pref) then begin
            idx = where(states.prefilter eq pref, Np)
            if Np eq 0 then begin
              print, inam+' : ERROR : '+cam+': no files matching prefilter '+pref
              continue
            endif
            files = files[idx]
            states = states[idx]
          endif 
        endif
        self->extractstates_nondb, files, states
        ;; Only one prefilter?
        if keyword_set(pref) then begin
          idx = where(states.prefilter eq pref, Np)
          if Np eq 0 then begin
            print, inam+' : ERROR : '+cam+': no files matching prefilter '+pref
            continue
          endif
          files = files[idx]
          states = states[idx]
        endif 
        ;; Check for complete scans only
        if ~keyword_set(all_data) then begin
          scans = states[uniq(states.scannumber, sort(states.scannumber))].scannumber
          Nscans = n_elements(scans)
          f_scan = lonarr(Nscans)
          for iscan = 0L, Nscans-1 do $
             f_scan[iscan] = n_elements(where(states.scannumber eq scans[iscan]))
          mask = replicate(1B, Nfiles)
          for iscan = 1L, Nscans-1 do begin
            if f_scan[iscan]-f_scan[0] lt 0 then begin
              print, inam + ' : WARNING : ' + cam + ': Incomplete scan nr ' + strtrim(scans[iscan], 2)
              print, inam + '             only ' + strtrim(f_scan[iscan], 2) + ' of ' $
                     + strtrim(f_scan(0), 2) + ' files.  Skipping it'
              mask[where(states.scannumber EQ scans[iscan])] = 0
            endif
          endfor
          idx = where(mask)
          files = (temporary(files))[idx]
          Nfiles = n_elements(files)
          states = states[idx]
        endif
      endelse 
      
      ;; Create linker script
      Nfiles = n_elements(files)

      linkername = linkerdir + cam + '_science_linker_'+folder_tag+'.sh'
      openw, lun, linkername, /get_lun
      printf, lun, '#!/bin/bash'      
      
      ;; Create folders
      outdir = self.out_dir + '/' + link_dir + '/' + folder_tag+ '/' + cam + '/'
      file_mkdir, outdir
      if iswb then begin
        outdir1 = self.out_dir + '/' + link_dir + '/' + folder_tag+ '/' + cam + '_nostate/'
        file_mkdir, outdir1
      endif
      
      for ifile = 0L, Nfiles - 1 do begin
        if uscan ne '' then if states.scannumber[ifile] NE uscan then continue

        namout = outdir + detector $
                 + '_' + string(states[ifile].scannumber, format = '(i05)') $
                 + '_' + strtrim(states[ifile].fullstate, 2) $
                 + '_' + string(states[ifile].framenumber, format = '(i07)')
        
        printf, lun, 'ln -sf '+ files[ifile] + ' ' + namout
        
        if iswb then begin
          namout = outdir1 + detector $
                   + '_' + string(states[ifile].scannumber, format = '(i05)') $
                   + '_' + strtrim(states[ifile].prefilter, 2) $
                   + '_' + string(states[ifile].framenumber, format = '(i07)')
          
          printf, lun, 'ln -sf '+ files[ifile] + ' ' + namout
        endif

        red_progressbar, ifile, Nfiles, inam+' : creating linker for '+cam
        
      endfor                    ; ifile
      
      free_lun, lun
      
      ;; Link data 
      print, inam + ' : executing '+  linkername
      spawn, 'chmod a+x ' + linkername
      spawn, '/bin/bash ' + linkername
            
    endfor                      ; icam
  endfor                        ; idir

end
