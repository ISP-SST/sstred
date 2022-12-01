; docformat = 'rst'

;+
; Make links to raw data in the form that is required for momfbd
; processing.
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
; 
; 
; :Params:
; 
; 
; :Keywords:
; 
;    link_dir : 
;   
;   
;    uscan : 
;   
;   
;    all_data : in, optional, type=boolean
;
;       Link not only complete sequences, but everything found 
;   
;    pref : in, optional, type=string
;
;       Process only prefilter 'pref'
; 
; :History:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2013-08-27 : MGL. Added support for logging. Let the subprogram
;                find out its own name.
; 
;   2013-09-11 : MGL. Fixed default values for keywords.
; 
;   2013-12-16 : PS. Keyword ALL_DATA, by default only link complete
;                scans keyword PREF
;
;   2014-01-02 : PS. Change subdirectory names to also use cam, not
;                camtag remove nremove and no_remove keywords (all
;                done in prepmomfbd)
;
;   2016-05-30 : MGL. Rename from red::link_data, make it a CHROMIS
;                method. Rewrite for CHROMIS.
;
;   2016-05-31 : MGL. Added dirs keyword. Don't zero the scannumber. 
;
;   2016-06-02 : MGL. Remove some keywords to extractstates.
;
;   2016-08-23 : THI. Rename camtag to detector and channel to camera,
;                so the names match those of the corresponding SolarNet
;                keywords.
;
;   2016-09-14 : MGL. Modify regex for finding files. Bugfix. Make
;                link script executable.
; 
;   2019-07-10 : OA. Rewritten to use database
;
;   2019-07-24 : OA. Added check for sst_db.
;
;   2022-08-21 : MGL. Version for CRISP2 based on the CHROMIS version.
;                Make nonstate links to WB data files, links to
;                directories for NB data.
;
;   2022-11-06 : MGL. This version now used for CRISP2 and CHROMIS
;                (2022-11-03 or later). Special case is mosaic data.
;   
;-
pro red::link_data, all_data = all_data $
                    , dirs = dirs $
                    , pref = pref $
                    , uscan = uscan $
                    , link_dir = link_dir

  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)

  if ((typename(self)).tolower()) eq 'chromis' $
     && self.isodate lt '2022-11-03' then begin
    ;; Use the old CHROMIS method for CHROMIS data with undecoded
    ;; wheel and hrz file names.
    self -> link_data_wheelhrz, all_data = all_data $
                                , dirs = dirs $
                                , pref = pref $
                                , uscan = uscan $
                                , link_dir = link_dir
    return
  end
  
  
  if n_elements(link_dir) eq 0 then link_dir = 'data'
  if n_elements(uscan) eq 0 then uscan = ''

  Ndirs = n_elements(dirs)
  if Ndirs eq 0 then begin
    ;; No dirs given in keyword, use default
    dirs = *self.data_dirs
    Ndirs = n_elements(dirs) 
  endif else begin
    for idir = 0, Ndirs-1 do begin
      if ~file_test(dirs[idir], /dir) then begin
        ;; This dir is not an existing directory.
        ;; Try to interpret as a selection from the default dirs.
        dirs[idir] = file_dirname((*self.data_dirs)[0])+'/'+dirs[idir]
      endif
    endfor                      ; idir
  endelse
  
  if Ndirs eq 0 then begin
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
      self->extractstates_db,'',datasets=dataset,sts
    endif
      

    for icam = 0L, Ncams-1 do begin
      
      cam = cams[icam]
      detector = self->getdetector( cam )
      iswb = strmatch(cam,'*-[DW]')
      
      ;; Make links to cam directories, the files in these have the NB
      ;; state info in them.
      file_mkdir,  self.out_dir + '/' + link_dir + '/' + folder_tag+ '/'
      file_link, /allow, data_dir + '/' + cam, self.out_dir + '/' + link_dir + '/' + folder_tag+ '/' 
      
      if iswb then begin

        ;; Wideband data, also make nostate links to individual files.
        
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
        outdir1 = self.out_dir + '/' + link_dir + '/' + folder_tag+ '/' + cam + '_nostate/'
        file_mkdir, outdir1
        
        for ifile = 0L, Nfiles - 1 do begin
          if uscan ne '' then if states.scannumber[ifile] NE uscan then continue
          
          namout = outdir1 + detector $
                   + '_' + string(states[ifile].scannumber, format = '(i05)') $
                   + '_' + strtrim(states[ifile].prefilter, 2) $
                   + '_' + string(states[ifile].framenumber, format = '(i07)') $
                   + '.fits'
          
          printf, lun, 'ln -sf '+ files[ifile] + ' ' + namout
          
          red_progressbar, ifile, Nfiles, inam+' : creating linker for '+cam
          
        endfor                  ; ifile
        
        free_lun, lun
        
        ;; Link data
        print, inam+' : executing '+  linkername
        spawn, 'chmod a+x ' + linkername
        spawn, '/bin/bash ' + linkername
      endif
      
    endfor                      ; icam
  endfor                        ; idir
  
end
