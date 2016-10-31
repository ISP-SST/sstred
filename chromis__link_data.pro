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
;    link_dir  : 
;   
;   
;    uscan  : 
;   
;   
;    all_data    : Not only link complete sequences, but everything found 
;   
;    pref        : Only process prefilter 'pref'
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
pro chromis::link_data, link_dir = link_dir $
                        , uscan = uscan $
                        , all_data = all_data $
                        , pref = pref $
                        , dirs = dirs $
                        , nremove = nremove

  if n_elements(link_dir) eq 0 then link_dir = 'data'
  if n_elements(uscan) eq 0 then uscan = ''

  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  ;; Logging
  help, /obj, self, output = selfinfo 
  red_writelog, selfinfo = selfinfo

  if n_elements(dirs) gt 0 then begin
     dirs = [dirs] 
  endif else begin
     if ~ptr_valid(self.data_dirs) then begin
        print, inam+' : ERROR : undefined data_dir'
        return
     endif
     dirs = *self.data_dirs
  endelse

  Ndirs = n_elements(dirs)
  if( Ndirs eq 0) then begin
     print, inam+' : ERROR : no directories defined'
     return
  endif else begin
     if Ndirs gt 1 then dirstr = '['+ strjoin(dirs,';') + ']' $
     else dirstr = dirs[0]
  endelse


  linkerdir = self.out_dir + '/' + 'link_scripts' + '/'
  file_mkdir, linkerdir
  
  ;; Create file list
  for idir = 0L, Ndirs - 1 do begin
     print, inam + ' : Folder -> ' + dirs[idir]
     data_dir = dirs[idir]
     folder_tag = strsplit(data_dir,'/',/extract)
     nn = n_elements(folder_tag) - 1
     folder_tag = folder_tag[nn]

     cams = *self.cameras
     Ncams = n_elements(cams)

     for icam = 0L, Ncams-1 do begin
        
        cam = cams[icam]
        detector = self->getdetector( cam )

        case cam of
           'Chromis-N' : wb = 0B
           'Chromis-W' : wb = 1B
           'Chromis-D' : wb = 1B
           else: stop
        endcase

        files = file_search(data_dir + '/' + cam + '/*cam*', count = Nfiles)
        
        if(files[0] eq '') then begin
           print, inam+' : ERROR : '+cam+': no files found in: '+$
                  data_dir +' : skipping camera!'
           continue
        endif

        ;; Sort files by image number
        files = red_sortfiles(files)
        
        ;; Get states
        self->extractstates, files, states

        ;; Only one prefilter?
        IF keyword_set(pref) THEN BEGIN
           idx = where(states.prefilter EQ pref, np)
           IF np EQ 0 THEN BEGIN
              print, inam+' : ERROR : '+cam+': no files matching prefilter '+pref
              CONTINUE
           ENDIF
           files = files[idx]
           Nfiles = np
           states = states[idx]
        ENDIF

        ;;; check for complete scans only
        IF ~keyword_set(all_data) THEN BEGIN
           scans = states[uniq(states.scannumber, sort(states.scannumber))].scannumber

           Nscans = n_elements(scans)
           f_scan = lonarr(Nscans)
           FOR iscan = 0L, Nscans-1 DO $
              f_scan[iscan] = n_elements(where(states.scannumber EQ scans[iscan]))
           mask = replicate(1b, Nfiles)
           FOR iscan = 1L, Nscans-1 DO BEGIN
              IF f_scan[iscan]-f_scan[0] LT 0 THEN BEGIN
                 print, inam + ' : WARNING : ' + cam + ': Incomplete scan nr ' + strtrim(scans[iscan], 2)
                 print, inam + '             only ' + strtrim(f_scan[iscan], 2) + ' of ' $
                        + strtrim(f_scan(0), 2) + ' files.  Skipping it'
                 mask[where(states.scannumber EQ scans[iscan])] = 0
              ENDIF
           ENDFOR
           idx = where(mask)
           files = (temporary(files))[idx]
           Nfiles = n_elements(files)
           states = states[idx]
        ENDIF

        ;; Flag nremove
            
        ;; Create linker script
        Nfiles = n_elements(files)

        linkername = linkerdir + cam + '_science_linker_'+folder_tag+'.sh'
        openw, lun, linkername, /get_lun
        printf, lun, '#!/bin/bash'
         
        outdir = self.out_dir + '/' + link_dir + '/' + folder_tag+ '/' + cam + '/'
        outdir1 = self.out_dir + '/' + link_dir + '/' + folder_tag+ '/' + cam + '_nostate/'

        ;; Create folders
        file_mkdir, outdir
        if wb then file_mkdir, outdir1

        for ifile = 0L, Nfiles - 1 do begin
;           if(stat.star[ifile]) then continue
           if uscan ne '' then if states.scannumber[ifile] NE uscan then continue
                                
           namout = outdir + detector $
                    + '_' + string(states[ifile].scannumber, format = '(i05)') $
                    + '_' + strtrim(states[ifile].fullstate, 2) $
                    + '_' + string(states[ifile].framenumber, format = '(i07)') $
                    + '.fits'
         
           printf, lun, 'ln -sf '+ files[ifile] + ' ' + namout
         
           if wb then begin
              namout = outdir1 + detector $
                       + '_' + string(states[ifile].scannumber, format = '(i05)') $
                       + '_' + strtrim(states[ifile].prefilter, 2) $
                       + '_' + string(states[ifile].framenumber, format = '(i07)') $
                       + '.fits'
              
              printf, lun, 'ln -sf '+ files[ifile] + ' ' + namout
           endif

           red_progressbar, ifile, Nfiles, inam+' : creating linker for '+cam, clock = clock
           
        endfor                  ; ifile
        free_lun, lun
         
        ;; Link data
        print, inam+' : executing '+  linkername
        spawn, 'chmod a+x ' + linkername
        spawn, '/bin/bash ' + linkername
        
     endfor                     ; icam
  endfor                        ; idir
  
end
