; docformat = 'rst'

;+
; 
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
;-
pro chromis::link_data, link_dir = link_dir $
                        , uscan = uscan $
                        , all_data = all_data $
                        , pref = pref $
                        , nremove = nremove

  if n_elements(link_dir) eq 0 then link_dir = 'data'
  if n_elements(uscan) eq 0 then uscan = ''

  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  ;; Logging
  help, /obj, self, output = selfinfo 
  red_writelog, selfinfo = selfinfo

  if ~ptr_valid(self.data_dirs) then begin
     print, inam+' : ERROR : undefined data_dir'
     return
  endif

  Ndirs = n_elements(*self.data_dirs)

  linkerdir = self.out_dir + '/' + 'link_scripts' + '/'
  file_mkdir, linkerdir
  
  ;; Create file list
  for idir = 0L, Ndirs - 1 do begin
     print, inam + ' : Folder -> ' + (*self.data_dirs)[idir]
     data_dir = (*self.data_dirs)[idir]
     folder_tag = strsplit(data_dir,'/',/extract)
     nn = n_elements(folder_tag) - 1
     folder_tag = folder_tag[nn]

     cams = *self.cam_channels
     Ncams = n_elements(cams)

     for icam = 0L, Ncams-1 do begin
        
        cam = cams[icam]
        camtag = self->getcamtag( cam )

        case cam of
           'Chromis-N' : wb = 0B
           'Chromis-W' : wb = 1B
           'Chromis-D' : wb = 1B
           else: stop
        endcase

        files = file_search(data_dir + '/' + cam + '/cam*', count = Nfiles)
        
        if(files[0] eq '') then begin
           print, inam+' : ERROR : '+cam+': no files found in: '+$
                  data_dir +' : skipping camera!'
           continue
        endif

        ;; Sort files by image number
        files = red_sortfiles(files)
        
        ;; Get states
        self->extractstates, files, states, /basename, /cam, /prefilter, /fullstate

        ;; Only one prefilter?
        IF keyword_set(pref) THEN BEGIN
           idx = where(states.prefilter EQ pref, np)
           IF np EQ 0 THEN BEGIN
              print, inam+' : ERROR : '+cam+': no files matching prefilter '+pref
              CONTINUE
           ENDIF
           files = files(idx)
           Nfiles = np
           states = states[idx]
        ENDIF


        ;;; check for complete scans only
        IF ~keyword_set(all_data) THEN BEGIN
            scans = states.scannumber[uniq(states.scannumber, sort(states.scannumber))]

           Nscans = n_elements(scans)
           f_scan = lonarr(Nscans)
           FOR iscan = 0L, Nscans-1 DO $
              f_scan[iscan] = n_elements(where(states.scannumber EQ scans[iscan]))
           idx = replicate(1b, Nfiles)
           FOR iscan = 1L, Nscans-1 DO BEGIN
              IF f_scan[iscan]-f_scan[0] LT 0 THEN BEGIN
                 print, inam+' : WARNING : '+cam+': Incomplete scan nr '+scans[iscan]
                 print, inam+'             only ' + strtrim(f_scan[iscan], 2) + ' of ' $
                        + strtrim(f_scan(0), 2) + ' files.  Skipping it'
                 idx[where(states.scannumber EQ scans[iscan])] = 0
              ENDIF
           ENDFOR
           files = (temporary(files))[where(idx)]
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

        red_progressbar, 0, Nfiles, message = inam+' : creating linker for '+cam

        for ifile = 0L, Nfiles - 1 do begin
;           if(stat.star[ifile]) then continue
           if uscan ne '' then if states.scannumber[ifile] NE uscan then continue
                                
           namout = outdir + camtag $
                    + '_' + string(states[ifile].scannumber, format = '(i05)') $
                    + '_' + states[ifile].fullstate $
                    + '_' + string(states[ifile].framenumber, format = '(i07)') $
                    + '.fits'
         
           printf, lun, 'ln -sf '+ files[ifile] + ' ' + namout
         
           if wb then begin
              namout = outdir1 + camtag $
                       + '_' + string(states[ifile].scannumber, format = '(i05)') $
                       + '_' + states[ifile].prefilter $
                       + '_' + string(states[ifile].framenumber, format = '(i07)') $
                       + '.fits'
              
              printf, lun, 'ln -sf '+ files[ifile] + ' ' + namout
           endif

           red_progressbar, ifile, Nfiles, message = inam+' : creating linker for '+cam
           
        endfor                  ; ifile
        free_lun, lun
        red_progressbar, message = inam+' : creating linker for '+cam, /finished
         
        ;; Link data
        print, inam+' : executing '+  linkername
;        spawn, '/bin/bash ' + linkername
        
     endfor                     ; icam
  endfor                        ; idir
  
end