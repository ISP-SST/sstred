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
; 
; :Keywords:
; 
;    split_dir : in, optional, type=string, default='data'
;
;        Name of directory under which the split data will be stored. 
;   
;    cameras : in, optional, type=strarr
;
;        The cameras (like Chromis-N) for which you want to split
;        data.
;    
;    filenums : in, optional, type=array
;   
;        Specify which files to split. 
;    
;    scannums : in, optional, type=array
;   
;        Specify which scans to split. 
;   
;    uscan : in, optional
;   
;   
;    all_data : Not only split complete sequences, but everything found 
;   
;    pref : Only process prefilter 'pref'
; 
; :History:
; 
;   2016-05-30 : MGL. New method, based on chromis::link_data.
;   
;   2016-05-31 : JLF. Start using red_keytab to keep track of SOLARNET
;                keywords. 
; 
;   2016-05-31 : MGL. Added dirs keyword. Link WB data on the
;                fly. Don't zero the scannumber. Bugfix.
;
;   2016-06-02 : MGL. Remove some keywords to extractstates.
;
;   2016-06-03 : MGL. Read data silently. Bugfix. Save frames as INTs. 
;
;   2016-08-23 : THI. Rename camtag to detector and channel to camera,
;                so the names match those of the corresponding SolarNet
;                keywords. 
;
;   2016-08-31 : MGL. Find files with new style names. Use header
;                keyword FRAMENUM. 
;
;   2016-09-01 : MGL. Remove unused keyword nremove. New keywords
;                cameras and filenums.
;
;   2016-09-05 : MGL. New keyword scannums.
;
;   2016-09-08 : MGL. Speed up.
;
;-
pro chromis::split_data, split_dir = split_dir $
                         , uscan = uscan $
                         , all_data = all_data $
                         , pref = pref $
                         , dirs = dirs $
                         , cameras = cameras $
                         , filenums = filenums $
                         , scannums = scannums

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
     dirs = *self.data_dirs
  endelse

  Ndirs = n_elements(dirs)
  if Ndirs eq 0 then begin
     print, inam+' : ERROR : no directories defined'
     return
  endif else begin
     if Ndirs gt 1 then dirstr = '['+ strjoin(dirs,';') + ']' $
     else dirstr = dirs[0]
  endelse

  if n_elements(cameras) eq 0 then cams = *self.cameras else cams = cameras
  Ncams = n_elements(cams)


  ;; Create file list
  for idir = 0L, Ndirs - 1 do begin
     print, inam + ' : Folder -> ' + dirs[idir]
     data_dir = dirs[idir]
     folder_tag = strsplit(data_dir,'/',/extract)
     nn = n_elements(folder_tag) - 1
     folder_tag = folder_tag[nn]

     for icam = 0L, Ncams-1 do begin
        
        cam = cams[icam]
        detector = self->getdetector( cam )

        print, inam + ' : Camera -> ' + cam


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

        ;; Get states
        self->extractstates, files, states

        if n_elements(filenums) ne 0 then begin
           ;; Filter the files list so we only split the wanted files. 
           match2, filenums, states.framenumber, suba, subb
           indx = where(suba ne -1) ; Should this be where(subb ne -1)?
           files = files[indx]
           states = states[indx]
           Nfiles = n_elements(indx)
           stop
        endif

        if n_elements(scannums) ne 0 then begin
           ;; Filter the files list so we only split the files that
           ;; belong to the specified scans. 
           match2, scannums, states.scannumber, suba, subb
           indx = where(subb ne -1)
           files = files[indx]
           states = states[indx]
           Nfiles = n_elements(indx)
        endif

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
                 print, inam+' : WARNING : '+cam+': Incomplete scan nr '+scans[iscan]
                 print, inam+'             only ' + strtrim(f_scan[iscan], 2) + ' of ' $
                        + strtrim(f_scan(0), 2) + ' files.  Skipping it'
                 mask[where(states.scannumber EQ scans[iscan])] = 0
              ENDIF
           ENDFOR
           idx = where(mask)
           files = (temporary(files))[idx]
           Nfiles = n_elements(files)
           states = states[idx]
        ENDIF

        
        Nfiles = n_elements(files)

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
           cube = red_readdata(files[ifile], header = head, /silent)
           
           dims = size(cube, /dim)
           Nframes = dims[2]

           ;; Get the frame number for the first frame in the cube,
           ;; remove from header

           xposure  = fxpar(head, 'XPOSURE')
           ;;frame1   = fxpar(head, red_keytab('frame'))
           frame1   = fxpar(head, 'FRAMENUM')
           date_beg = fxpar(head, 'DATE-BEG')
           cadence  = fxpar(head, 'CADENCE')
           
           date = (strsplit(date_beg, 'T', /extract))[0]
           time_beg = (strsplit(date_beg, 'T', /extract))[1]

           ;; Delete and modify keywords so the header can be used for
           ;; the individual-frame files.
           sxaddpar, head, 'DATE', red_timestamp(/utc, /iso) $
                     , 'Creation date of FITS header', before = 'TIMESYS'
           sxdelpar, head, 'CADENCE'
           sxdelpar, head, red_keytab('frame')
;           sxdelpar, head, 'NAXIS3'
;           sxaddpar, head, 'NAXIS', 2, /savecomment

           check_fits, fix(cube[*, *, 0]), head, /UPDATE, /SILENT

           for iframe = 0L, Nframes-1 do begin

              frameno = frame1 + iframe
              sxaddpar, head, red_keytab('framenumber'), frameno, before = 'COMMENT'

              ;; DATE-BEG, DATE-AVE, DATE-END
              sxaddpar, head, 'DATE-BEG', /savecomment, after = 'TIMESYS' $
                        , date + 'T' + red_time2double(time_beg+iframe*cadence, /dir)
              sxaddpar, head, 'DATE-END', /savecomment, after = 'TIMESYS' $
                        , date + 'T' + red_time2double(time_beg+iframe*cadence+xposure, /dir)
              sxaddpar, head, 'DATE-AVE', /savecomment, after = 'TIMESYS' $
                        , date + 'T' + red_time2double(time_beg+iframe*cadence+xposure/2, /dir)
              
              
              ;; Write one frame
              if wb then statestring = strtrim(states[ifile].fpi_state, 2) $
              else statestring = strtrim(states[ifile].fullstate, 2)
              
              namout = outdir + detector $
                       + '_' + string(states[ifile].scannumber, format = '(i05)') $
                       + '_' + statestring $
                       + '_' + string(frameno, format = '(i07)') $
                       + '.fits'
              
              red_writedata, namout, fix(cube[*, *, iframe]), header = head $
                             , filetype = 'FITS', /overwrite
              
              if wb then begin

                 ;; Link name 
                 namout1 = outdir1 + detector $
                          + '_' + string(states[ifile].scannumber, format = '(i05)') $
                          + '_' + strtrim(states[ifile].prefilter, 2) $
                          + '_' + string(frameno, format = '(i07)') $
                          + '.fits'
 
                 ;; Do the linking
                 file_link, namout, namout1, /allow_same

              endif
              
           endfor               ; iframe

           red_progressbar, ifile, Nfiles, message = inam+' : splitting files for '+cam
           
        endfor                  ; ifile

        red_progressbar, message = inam+' : splitting files for '+cam, /finished
        
     endfor                     ; icam
  endfor                        ; idir
  
end
