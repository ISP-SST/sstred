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
;    link_dir : in, optional
;   
;   
;    uscan : in, optional
;   
;   
;    all_data : in, optional, type=boolean
;
;       Link not only complete sequences, but everything found 
;   
;    pref : in, optional, type=string
;
;       Only process prefilter 'pref'
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
;   2016-05-30 : MGL. Rename from red::link_data, make it a CRISP
;                method.
; 
;   2018-04-12 : MGL. Adapt to new codebase.
;  
;-
pro crisp::link_data, all_data = all_data $
                      , dirs = dirs $
                      , nremove=nremove $
                      , pref = pref $
                      , uscan = uscan $
                      , link_dir = link_dir 

  if n_elements(link_dir) eq 0 then link_dir = 'data'
  if n_elements(uscan) eq 0 then uscan = ''

  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)

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
   endfor     ;idir
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
      iswb = strmatch(cams,'*-W') or strmatch(cams,'*-D')

;       case icam of
;           0: begin
;              if(self.docamt eq 0B) then begin
;                 print, inam+' : Nothing to do for ', self.camt
;              endif
;              cam = self.camt
;              doit = self.docamt
;              wb = 0B
;           end
;           1: begin
;              if(self.docamr eq 0B) then begin
;                 print, inam+' : Nothing to do for '+ self.camr
;              endif
;              cam = self.camr
;              doit = self.docamr
;              wb = 0B
;           end
;           2: begin
;              if(self.docamwb eq 0B) then begin
;                 print, inam+' : Nothing to do for '+ self.camwb
;              endif
;              cam = self.camwb
;              doit = self.docamwb
;              wb = 1B
;           end
;        endcase
      
;      if(~doit) then continue

      files = file_search(data_dir + '/' + cam + '/cam*', count = Nfiles)
      files = files(where(strpos(files, '.lcd.') LT 0, nf))
      Nfiles = n_elements(files)
      
      if(files[0] eq '') then begin
        print, inam+' : ERROR : '+cam+': no files found in: '+$
               data_dir +' : skipping camera!'
        continue
      endif

      ;; Sort files by image number
;      files = red_sortfiles(files)
      
      ;; Get states
      self->extractstates, files, states
;     stat = red_getstates(files)
      

      ;; only one prefilter?
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
;        scans = stat.scan(uniq(stat.scan))
        
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
;        Nscans = n_elements(scans)
;        f_scan = lonarr(Nscans)
;        for iscan = 0L, Nscans-1 do $
;           f_scan[iscan] = n_elements(where(states.scannumber eq scans[iscan]))
;        idx = replicate(1b, nf)
;        FOR iscan = 1, Nscans-1 DO BEGIN
;          IF f_scan(i)-f_scan(0) LT 0 THEN BEGIN
;            print, inam+' : WARNING : '+cam+': Incomplete scan nr '+scans(i)
;            print, inam+'             only '+strtrim(f_scan(i), 2)+' of '+strtrim(f_scan(0), 2)+' files.  Skipping it'
;            idx(where(stat.scan EQ scans(i))) = 0
;          ENDIF
;        ENDFOR
;              ;;; re-generate stats, numbers
;        files = (temporary(files))(where(idx))
;        nf = n_elements(files)
;        stat = red_getstates(files)
;      endif
      
      
      ;; Flag nremove
;      red_flagtuning, stat, nremove ; Did we move this functionality to momfbd? red_setupworkdir_crips doesn't call link_data with nremove set anyway!
      
      ;; Create linker script
      Nfiles = n_elements(files)

      linkername = linkerdir + cam + '_science_linker_'+folder_tag+'.sh'
      openw, lun, linkername, /get_lun
      printf, lun, '#!/bin/bash'
      
;      nt = n_elements(files)
;      camtag = (strsplit(file_basename(files[0]), '.', /extract))[0]
      
;      linkername = self.out_dir + '/' + cam + '_science_linker_'+folder_tag+'.sh'
;      openw, lun, linkername, /get_lun
;      printf, lun, '#!/bin/bash'
      
;     ;; Print links
;      Ntot = 100. / (Nt - 1.0)
;      bb = string(13b)
      
      ;; Create folders
      outdir = self.out_dir + '/' + link_dir + '/' + folder_tag+ '/' + cam + '/'
      file_mkdir, outdir
      if iswb[icam] then begin
        outdir1 = self.out_dir + '/' + link_dir + '/' + folder_tag+ '/' + cam + '_nostate/'
        file_mkdir, outdir1
      endif
      
      for ifile = 0L, Nfiles - 1 do begin
;           if(stat.star[ifile]) then continue
        if uscan ne '' then if states.scannumber[ifile] NE uscan then continue

        namout = outdir + detector $
                 + '_' + string(states[ifile].scannumber, format = '(i05)') $
                 + '_' + strtrim(states[ifile].fullstate, 2) $
                 + '_' + string(states[ifile].framenumber, format = '(i07)')
        
        printf, lun, 'ln -sf '+ files[ifile] + ' ' + namout
        
        if iswb[icam] then begin
          namout = outdir1 + detector $
                   + '_' + string(states[ifile].scannumber, format = '(i05)') $
                   + '_' + strtrim(states[ifile].prefilter, 2) $
                   + '_' + string(states[ifile].framenumber, format = '(i07)')
          
          printf, lun, 'ln -sf '+ files[ifile] + ' ' + namout
        endif

        red_progressbar, ifile, Nfiles, inam+' : creating linker for '+cam
        
      endfor                    ; ifile
      
;      for ii = 0L, nt - 1 do begin
;        if(stat.star[ii]) then continue
;        if uscan ne '' then if stat.scan[ii] NE uscan then continue
;                                ;
;        namout = outdir + camtag $
;                  + '.' + stat.scan[ii] $
;                  + '.' + stat.state[ii] $
;                  + '.' + stat.nums[ii]
;        
;        printf, lun, 'ln -sf '+ files[ii] + ' ' + namout
;        
;        if(wb) then begin
;          namout = outdir1 + camtag + '.' + stat.scan[ii]+ '.' +stat.pref[ii] + '.' +stat.nums[ii]
;          printf, lun, 'ln -sf '+ files[ii] + ' ' + namout
;        endif
;        
;        print, bb, inam+' : creating linker for '+cam+$
;               ' -> ', ii * ntot, '%', FORMAT = '(A,A,F5.1,A,$)'
;      endfor
      free_lun, lun
      
      ;; Link data 
      print, inam + ' : executing '+  linkername
      spawn, 'chmod a+x ' + linkername
      spawn, '/bin/bash ' + linkername
      
      
    endfor                      ; icam
  endfor                        ; idir
;      ;; Create folder and link data
;      file_mkdir, outdir
;      if(wb) then file_mkdir, outdir1
;      
;      print, inam+' : executing '+  linkername
;      spawn, '/bin/bash ' + linkername
;      
;    endfor
;  endfor
;  
;  return
end
