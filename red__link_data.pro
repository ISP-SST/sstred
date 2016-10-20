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
; 
; :returns:
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
;                scans keyword PREF.
;
;   2014-01-02 : PS. Change subdirectory names to also use cam, not
;                camtag remove nremove and no_remove keywords (all
;                done in prepmomfbd).
;   
;-
pro red::link_data, link_dir = link_dir, uscan = uscan, ALL_DATA = all_data, PREF = pref, nremove=nremove

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
  
  ;; Create file list
  for ff = 0L, Ndirs - 1 do begin
     print, inam + ' : Folder -> ' + self.data_dirs[ff]
     data_dir = self.data_dirs[ff]
     folder_tag = strsplit(data_dir,'/',/extract)
     nn = n_elements(folder_tag) - 1
     folder_tag = folder_tag[nn]

     for ic = 0L, 2 do begin
        case ic of
           0: begin
              if(self.docamt eq 0B) then begin
                 print, inam+' : Nothing to do for ', self.camt
              endif
              cam = self.camt
              doit = self.docamt
              wb = 0B
           end
           1: begin
              if(self.docamr eq 0B) then begin
                 print, inam+' : Nothing to do for '+ self.camr
              endif
              cam = self.camr
              doit = self.docamr
              wb = 0B
           end
           2: begin
              if(self.docamwb eq 0B) then begin
                 print, inam+' : Nothing to do for '+ self.camwb
              endif
              cam = self.camwb
              doit = self.docamwb
              wb = 1B
           end
        endcase
        
        if(~doit) then continue

        files = file_search(data_dir + '/' + cam + '/camX*')
        files = files(where(strpos(files, '.lcd.') LT 0, nf))
        
        if(files[0] eq '') then begin
           print, inam+' : ERROR : '+cam+': no files found in: '+$
                  data_dir +' : skipping camera!'
           continue
        endif

        ;; Sort files by image number
        files = red_sortfiles(files)
                                
        ;; Get states
        stat = red_getstates(files)
        

        ;; only one prefilter?
        IF keyword_set(pref) THEN BEGIN
            idx = where(stat.pref EQ pref, np)
            IF np EQ 0 THEN BEGIN
                print, inam+' : ERROR : '+cam+': no files matching prefilter '+pref
                CONTINUE
            ENDIF
            files = files(idx)
            nf = np
            stat = red_getstates(files)
        ENDIF


        ;;; check for complete scans only
        IF ~keyword_set(all_data) THEN BEGIN
            scans = stat.scan(uniq(stat.scan))
            n_scans = n_elements(scans)
            f_scan = lonarr(n_scans)
            FOR i = 0, n_scans-1 DO $
              f_scan(i) = n_elements(where(stat.scan EQ scans(i)))
            idx = replicate(1b, nf)
            FOR i = 1, n_scans-1 DO BEGIN
                IF f_scan(i)-f_scan(0) LT 0 THEN BEGIN
                    print, inam+' : WARNING : '+cam+': Incomplete scan nr '+scans(i)
                    print, inam+'             only '+strtrim(f_scan(i), 2)+' of '+strtrim(f_scan(0), 2)+' files.  Skipping it'
                    idx(where(stat.scan EQ scans(i))) = 0
                ENDIF
            ENDFOR
              ;;; re-generate stats, numbers
            files = (temporary(files))(where(idx))
            nf = n_elements(files)
            stat = red_getstates(files)
        ENDIF
        
        
        ;; Flag nremove
        red_flagtuning, stat, nremove
            
        ;; Create linker script
        nt = n_elements(files)
        camtag = (strsplit(file_basename(files[0]), '.', /extract))[0]
        
        linkername = self.out_dir + '/' + cam + '_science_linker_'+folder_tag+'.sh'
        openw, lun, linkername, /get_lun
        printf, lun, '#!/bin/bash'
         
        ;; Print links
        ntot = 100. / (nt - 1.0)
        bb = string(13b)
         
        outdir = self.out_dir + '/' + link_dir + '/' + folder_tag+ '/' + cam + '/'
        outdir1 = self.out_dir + '/' + link_dir + '/' + folder_tag+ '/' + cam + '_nostate/'
         
        for ii = 0L, nt - 1 do begin
           if(stat.star[ii]) then continue
           if uscan ne '' then if stat.scan[ii] NE uscan then continue
                                ;
           namout = outdir + camtag + '.' + stat.scan[ii] +'.'+ stat.state[ii] + '.' +stat.nums[ii]
         
           printf, lun, 'ln -sf '+ files[ii] + ' ' + namout
         
           if(wb) then begin
              namout = outdir1 + camtag + '.' + stat.scan[ii]+ '.' +stat.pref[ii] + '.' +stat.nums[ii]
              printf, lun, 'ln -sf '+ files[ii] + ' ' + namout
           endif
         
           print, bb, inam+' : creating linker for '+cam+$
                  ' -> ', ii * ntot, '%', FORMAT = '(A,A,F5.1,A,$)'
        endfor
        free_lun, lun
        print, ' '
         
        ;; Create folder and link data
        file_mkdir, outdir
        if(wb) then file_mkdir, outdir1
         
        print, inam+' : executing '+  linkername
        spawn, '/bin/bash ' + linkername
         
     endfor
  endfor
         
  return
end
