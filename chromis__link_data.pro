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
;-
pro chromis::link_data, all_data = all_data $
                        , dirs = dirs $
                        , nremove = nremove $
                        , pref = pref $
                        , uscan = uscan $
                        , link_dir = link_dir 

  if n_elements(link_dir) eq 0 then link_dir = 'data'
  if n_elements(uscan) eq 0 then uscan = ''

  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)

  if n_elements(dirs) eq 0 then begin
    ;; No dirs given in keyword, use default
    dirs = *self.data_dirs
  endif else begin
    if ~file_test(dirs, /dir) then begin
      ;; The keyword doesn't point to an existing directory.
      ;; Try to interpret as a selection from the default dirs.
      dirs = file_dirname((*self.data_dirs)[0])+'/'+dirs
    endif
  endelse
  
  Ndirs = n_elements(dirs)
  if( Ndirs eq 0) then begin
    print, inam+' : ERROR : no directories defined'
    return
  endif ;else begin
;    if Ndirs gt 1 then dirstr = '['+ strjoin(dirs,';') + ']' $
;    else dirstr = dirs[0]
;  endelse

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

    for icam = 0L, Ncams-1 do begin
      
      cam = cams[icam]
      detector = self->getdetector( cam )
      iswb = strmatch(cam,'*-[DW]')

;      case cam of
;        'Chromis-N' : iswb = 0B
;        'Chromis-W' : iswb = 1B
;        'Chromis-D' : iswb = 1B
;        else: stop
;      endcase

      files = file_search(data_dir + '/' + cam + '/*cam*', count = Nfiles)
      
      if(files[0] eq '') then begin
        print, inam+' : ERROR : '+cam+': no files found in: ' $
               + data_dir + ' : skipping camera!'
        continue
      endif

      ;; Sort files by image number
;      files = red_sortfiles(files)
      
      ;; Get states
      self->extractstates, files, states

      ;; Only one prefilter?
      if keyword_set(pref) then begin
        idx = where(states.prefilter eq pref, Np)
        if Np eq 0 then begin
          print, inam+' : ERROR : '+cam+': no files matching prefilter '+pref
          continue
        endif
        Nfiles = Np
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
        endfor                  ; iscan
        idx = where(mask)
        files = (temporary(files))[idx]
        Nfiles = n_elements(files)
        states = states[idx]
      endif

      ;; Flag nremove
      
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
;           if(stat.star[ifile]) then continue
        if uscan ne '' then if states.scannumber[ifile] NE uscan then continue
        
        namout = outdir + detector $
                 + '_' + string(states[ifile].scannumber, format = '(i05)') $
                 + '_' + strtrim(states[ifile].fullstate, 2) $
                 + '_' + string(states[ifile].framenumber, format = '(i07)') $
                 + '.fits'
        
        printf, lun, 'ln -sf '+ files[ifile] + ' ' + namout
        
        if iswb then begin
          namout = outdir1 + detector $
                   + '_' + string(states[ifile].scannumber, format = '(i05)') $
                   + '_' + strtrim(states[ifile].prefilter, 2) $
                   + '_' + string(states[ifile].framenumber, format = '(i07)') $
                   + '.fits'
          
          printf, lun, 'ln -sf '+ files[ifile] + ' ' + namout
        endif

        red_progressbar, ifile, Nfiles, inam+' : creating linker for '+cam
        
      endfor                    ; ifile
      free_lun, lun
      
      ;; Link data
      print, inam+' : executing '+  linkername
      spawn, 'chmod a+x ' + linkername
      spawn, '/bin/bash ' + linkername
      
    endfor                      ; icam
  endfor                        ; idir
  
end
