
; docformat = 'rst'

;+
; Split data cube CHROMIS files into single-frame files.
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, ISP
; 
; :Params:
; 
; 
; :Keywords:
; 
;    split_dir  : 
;   
;   
;    uscan  : 
;   
;   
;    all_data    : Not only split complete sequences, but everything found 
;   
;    pref        : Only process prefilter 'pref'
; 
; :History:
; 
;   2016-05-30 : MGL. New method, based on chromis::link_data.
; 
;   2016-05-31 : MGL. Added dirs keyword.
;   
;-
pro chromis::split_data, split_dir = split_dir $
                         , uscan = uscan $
                         , all_data = all_data $
                         , pref = pref $
                         , dirs = dirs $
                         , nremove = nremove

  if n_elements(split_dir) eq 0 then split_dir = 'data'
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
     dirs = *self.dark_dir
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
        
        ;; Create linker script for WB data
        Nfiles = n_elements(files)

        linkername = linkerdir + cam + '_science_linker_' + folder_tag + '.sh'
        openw, lun, linkername, /get_lun
        printf, lun, '#!/bin/bash'
        
        outdir  = self.out_dir + '/' + split_dir + '/' + folder_tag+ '/' + cam + '/'
        outdir1 = self.out_dir + '/' + split_dir + '/' + folder_tag+ '/' + cam + '_nostate/'

        ;; Create folders
        file_mkdir, outdir
        if wb then file_mkdir, outdir1

        red_progressbar, 0, Nfiles, message = inam+' : splitting files for '+cam

        for ifile = 0L, Nfiles - 1 do begin
;           if(stat.star[ifile]) then continue
           if uscan ne '' then if states.scannumber[ifile] NE uscan then continue

           ;; Read the data cube
           cube = red_readdata(files[ifile], header = head)

           dims = size(cube, /dim)
           Nframes = dims[2]

           ;; Get the frame number for the first frame in the cube,
           ;; remove from header
           frame1   = fxpar(head, 'FRAME1')
           xposure  = fxpar(head, 'XPOSURE')
           date_beg = fxpar(head, 'DATE-BEG')
           cadence  = fxpar(head, 'CADENCE')
           
           date = (strsplit(date_beg, 'T', /extract))[0]
           time_beg = (strsplit(date_beg, 'T', /extract))[1]

           ;; Delete and modify keywords so the header can be used for
           ;; the individual-frame files.
           sxaddpar, head, 'DATE', red_timestamp(/utc, /iso) $
                     , 'Creation date of FITS header', before = 'TIMESYS'
           sxdelpar, head, 'FRAME1'
           sxdelpar, head, 'CADENCE'
           sxdelpar, head, 'NAXIS3'
           sxaddpar, head, 'NAXIS', 2, /savecomment

           for iframe = 0L, Nframes-1 do begin

              frameno = frame1 + iframe
              
              ;; DATE-BEG, DATE-AVE, DATE-END
              sxaddpar, head, 'DATE-BEG', /savecomment, after = 'TIMESYS' $
                        , date + 'T' + red_time2double(time_beg+iframe*cadence, /dir)
              sxaddpar, head, 'DATE-END', /savecomment, after = 'TIMESYS' $
                        , date + 'T' + red_time2double(time_beg+iframe*cadence+xposure, /dir)
              sxaddpar, head, 'DATE-AVE', /savecomment, after = 'TIMESYS' $
                        , date + 'T' + red_time2double(time_beg+iframe*cadence+xposure/2, /dir)
              
              
              ;; Write one frame
              
              namout = outdir + camtag $
                       + '_' + string(states[ifile].scannumber, format = '(i05)') $
                       + '_' + states[ifile].fullstate $
                       + '_' + string(frameno, format = '(i07)') $
                       + '.fits'
              
              red_writedata, namout, cube[*, *, iframe], header = head, filetype = 'FITS'
              
              if wb then begin

                 ;; Write link command 
                 namout = outdir1 + camtag $
                          + '_' + string(states[ifile].scannumber, format = '(i05)') $
                          + '_' + states[ifile].prefilter $
                          + '_' + string(frameno, format = '(i07)') $
                          + '.fits'
                 
                 printf, lun, 'ln -sf '+ files[ifile] + ' ' + namout

              endif
              
           endfor               ; iframe

           red_progressbar, ifile, Nfiles, message = inam+' : splitting files for '+cam
           
        endfor                  ; ifile

        red_progressbar, message = inam+' : creating linker for '+cam, /finished

        free_lun, lun
        
        ;; Link data
        print, inam+' : executing '+  linkername
;        spawn, '/bin/bash ' + linkername
        
     endfor                     ; icam
  endfor                        ; idir
  
end
