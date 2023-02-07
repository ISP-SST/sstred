; docformat = 'rst'

;+
; Make an adjusted version of a fitscube that is suitable for
; exporting to an SVO, and copy it to its archival directory.
;
; Copy also the spectral cube, if there is one.
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, Institute for Solar Physics
; 
; 
; :Params:
; 
;    filename : in, type=string
; 
;      The path of the fitscube file. 
; 
; 
; :Keywords:
;
;    iframe : in, optional, type=integer, default="based on ituning, istokes, and iscan"
;
;       The frame index in the data cube seen as a 3D cube, to be used
;       for generating a thumbnail. Takes precedence over iscan,
;       istokes, and ituning.
;
;    iscan : in, optional, type=integer, default="Selected for contrast."
;
;       The scan index, used to calculate iframe, to be used for
;       generating a thumbnail.
;   
;    istokes  : in, optional, type=integer, default=0
;
;       The stokes index, used to calculate iframe, to be used for
;       generating a thumbnail and video.
;
;    ituning  : in, optional, type=integer, default="Selected for brightness"
;
;       The tuning index, used to calculate iframe, to be used for
;       generating a thumbnail and video.
; 
;    smooth_width: in, optional, type=integer
;  
;      If present, blur all images in the exported cube with
;      smooth(frame, blur). This is to facilitate making sample cubes
;      without actually releasing data. If unity, interpreted as a
;      boolean and set to a default width of 11 pixels.
;
;    keywords: in, optional, type=struct
;
;      Keywords and values to be added to the header. Keywords are
;      upcased and checked against a list of SOLARNET-approved
;      keywords.
;
;    no_spectral_file: in, optional, type=boolean
;
;      Do not copy the sp version of the input file, if there is one.
;
;    no_thumbnail: in, optional, type=boolean
;
;      Do not make a thumbnail image.
;
;    no_video: in, optional, type=boolean
;
;      Do not make a video.
;
;    outdir: in, out, optional, type=string, default='/storage_new/science_data/YYYY-MM-DD/'
;   
;      The directory in which to store the new file/to which the file
;      was exported.
;
;    outfile: out, optional, type=string
;   
;      The name of the file to which the data were exported.
;
;    overwrite: in, optional, type=boolean
;
;      Set this to overwrite an existing file.
;
;    release_date : in, optional, type=string
;
;      The value of the RELEASE keyword, a date after which the data
;      are not proprietary.
;
;    release_comment : in, optional, type=string
;
;      The value of the RELEASEC keyword, a comment with details about
;      the proprietary state of the data, and whom to contact.
;
;    svo_api_key : in, optional, type=string
;
;      Ingest the metadata of the cube into the prototype SVO using
;      this API_KEY.
;
;    svo_username : in, optional, type=string
;
;      Ingest the metadata of the cube into the prototype SVO using
;      this USERNAME.
;
;    thumb_shrink_fac : in, optional, type=integer, default=4
;
;      If making a thumbnail image, shrink it by this factor.
;
;    video_shrink_fac : in, optional, type=integer, default=4
;
;      If making a video, shrink the FOV by this factor.
;
;    trust_datasum : in, optional, typy=boolean
;
;      Set this keyword to avoid checking datasum.
;
;
; :History:
; 
;   2019-06-12 : MGL. First version.
; 
;   2019-06-17 : MGL. New keyword help.
; 
;   2019-06-18 : MGL. New keywords smooth_width and outfile. Add
;                prstep header info.
; 
;   2019-07-04 : MGL. Optionally (opt-out) add checksum and datasum.
; 
;   2019-10-18 : MGL. New keywords svo_api_key and svo_username.
; 
;   2020-11-02 : MGL. Remove stuff that is now in fitscube_finalize. 
; 
;   2020-12-04 : MGL. Generate an OID. 
; 
;   2020-12-07 : MGL. Generate thumbnail image and video.
; 
;   2020-12-10 : MGL. New keywords thumb_shrink_fac and
;                video_shrink_fac. 
;
;   2021-01-11 : OA. Submit information about exported cubes to the
;                database.
;
;   2022-03-24 : OA. Added trust_datasum keyword.
; 
;-
pro red::fitscube_export, filename $
                          , help = help $
                          , keywords = keywords $
                          , no_spectral_file = no_spectral_file $
                          , no_wb_file = no_wb_file $
                          , no_thumbnail = no_thumbnail $
                          , no_video = no_video $
                          , outdir = outdir $
                          , outfile = outfile $
                          , overwrite = overwrite $
                          , smooth_width = smooth_width $
                          , svo_api_key = svo_api_key $
                          , svo_username = svo_username $
                          , thumb_shrink_fac = thumb_shrink_fac $
                          , video_shrink_fac = video_shrink_fac $
                          , allowed_users = allowed_users $
                          , swedish_data = swedish_data $
                          , trust_datasum = trust_datasum 
  
  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)
  ;; Make prpara
  red_make_prpara, prpara, keywords     
  red_make_prpara, prpara, no_spectral_file    
  red_make_prpara, prpara, smooth_width

  ;; Get the current date
  time = Systime(UTC=Keyword_Set(utc))
  day = String(StrMid(time, 8, 2), Format='(I2.2)') ; Required because UNIX and Windows differ in time format.
  month = StrUpCase(Strmid(time, 4, 3))
  year = Strmid(time, 20, 4)
  months = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']
  m = (Where(months EQ month)) + 1
  current_date =  year + '-' + string(m, FORMAT='(I2.2)') + '-' + day

  hdr = red_readhead(filename)
    
  release_date = fxpar(hdr,'RELEASE')
  if current_date lt release_date then begin
     if ~keyword_set(allowed_users) then begin
        allowed_users = ''
        reply = ''
        print, 'This data should be protected from unauthorized downloading.'
        print, "You have to enter users' e-mails."
        while reply ne 'x' do begin
           read,"Enter a user's e-mail (put 'x' to escape) >  ", reply
           allowed_users += reply+';'
        endwhile
        allowed_users = strmid(allowed_users, 0, strlen(allowed_users)-2)
      endif
     print, "Don't forget to ask the administrator to add user to the database and generate password."
     if allowed_users eq '' then begin
        print, 'You have to provide at least one username to make the cube available for downloading.'
        return
     endif
  endif

  if n_elements(smooth_width) gt 0 then begin
    ;; We need to implement reading/writing frames from/to spectral
    ;; fitscubes for this to work easily for spectral cubes.
    no_spectral_file = 1
    ;; If smooth_width is set to the inconsequential value of unity,
    ;; interpret it as a boolean and set to a default width of 11.
    if smooth_width eq 1 then smooth_width = 11
  endif
  
  if n_elements(thumb_shrink_fac) eq 0 then thumb_shrink_fac = 4
  if n_elements(video_shrink_fac) eq 0 then video_shrink_fac = 4

  ;; This will be used for a new DATE header and also to indicate the
  ;; version of the cube as part of the file name.
  new_DATE = red_timestamp(/utc,/iso)

  indir  = file_dirname(filename)+'/'
  infile = file_basename(filename)
  
  date_beg_split = strsplit(fxpar(hdr,'DATE-BEG'),'T',/extract)
  date_beg = date_beg_split[0]
  time_beg = (strsplit(date_beg_split[1],'.',/extract))[0]
  time_end = (strsplit(fxpar(hdr,'DATE-END'),'T',/extract))[1]
  time_end = strmid(time_end,0,8)

  ;; We could do the outdir default per site, just like we do for the
  ;; raw data directories.
  if n_elements(outdir) eq 0 then outdir = '/storage_new/science_data/' $
                                           + date_beg + '/' $
                                           + strtrim(fxpar(hdr,'INSTRUME'),2) + '/'

  ;; Any old versions of this file?
  spl=strsplit(infile,'_',/extract)
  srch='*'+strjoin(spl[1:5],'_')+'*'
  oldfiles = file_search(outdir + '/' + srch, count = Noldfiles)
  if Noldfiles gt 0 then begin
    s = ''
    print, inam + ' : There are '+strtrim(Noldfiles, 2)+' older versions of exported files in '+outdir+':'
    print, file_basename(oldfiles), format = '(a0)'
    read, '    Delete them [N]? ', s
    if s eq '' then s = 'N'
    if strupcase(strmid(s, 0, 1)) eq 'Y' then begin
      file_delete, oldfiles
      ;; Should perhaps also add code that removes the deleted files'
      ;; metadata from the SVO?
    endif
  endif

  if n_elements(smooth_width) gt 0 then begin
    outfile = red_strreplace(infile, '_im.fits', '_blurred_export'+new_DATE+'_im.fits')
  endif else begin
    outfile = red_strreplace(infile, '_im.fits', '_export'+new_DATE+'_im.fits')
  endelse

  if file_same(filename, outdir+outfile) then begin
    print, inam + ' : This operation would overwrite the input file.'
    return
  endif
  
  if file_test(outdir+'/'+outfile) and ~keyword_set(overwrite) then begin
    print, inam + ' : The output file already exists: '
    spawn, 'ls -l '+outdir+'/'+outfile
    return
  endif

  ;; Delete header keywords that should not be there
  red_fitsdelkeyword, hdr, 'STATE'

  ;; Add keywords after this one:
  anchor = 'TIMESYS'

  ;; Update the DATE keyword
  red_fitsaddkeyword, anchor = anchor, hdr, 'DATE', new_DATE
  red_fitsaddkeyword, anchor = anchor, hdr, 'FILENAME', outfile

  file_mkdir, outdir

  if ~keyword_set(no_wb_file) then begin

     indx = where(strmatch(hdr, '*Align reference:*'), Nwhere)     
     if Nwhere eq 0 then begin
        print, 'There is no information about align reference in the header.'
        print, 'Check the header (may be you should run red_fitscube_correct_prstep)'
        print, 'or rerun fitscube_export with /no_wb_file option.'
        return
     endif
     key = strtrim((strsplit(hdr[indx[-1]], '=', /extract))[0], 2)
     
     wcfile = red_strreplace(fxpar(hdr, key), 'Align reference: ', '')
     if file_dirname(wcfile) eq '.' then begin
        steps = fxpar(hdr,'PRSTEP*')
        if where(strmatch(steps,'*DATA-CURATION*')) ne -1 then $
          wcfile = 'cubes_converted/' + wcfile $
        else $
          wcfile = 'cubes_wb/' + wcfile
     endif
     if ~file_test(wcfile) then begin
        print, 'WB file ',wcfile,' does not exist.'
        print, 'Make the cube or rerun fitscube_export with /no_wb_file option.'
        return
     endif
     
     wboutfile = red_strreplace(file_basename(wcfile), '_im.fits', '_export'+new_DATE+'_im.fits')     
     red_fitsaddkeyword, hdr, key $
                         , 'Align reference: '+wboutfile $
                         , 'WB cube file name'
     wbhdr = headfits(wcfile)
     red_fitsdelkeyword, wbhdr, 'STATE'
     red_fitsaddkeyword, anchor = 'TIMESYS', wbhdr, 'DATE', new_DATE
     red_fitsaddkeyword, anchor = 'TIMESYS', wbhdr, 'FILENAME', wboutfile
     print, inam + ' : Copying the wideband cube...'
     file_copy, wcfile, outdir+wboutfile
     red_fitscube_newheader, outdir+wboutfile, wbhdr
     print, inam + ' : Wrote '+outdir+wboutfile
     red_fitscube_checksums, outdir+wboutfile, trust_datasum = trust_datasum 
  endif

  ;; Copy the file and write the new header
  print, inam + ' : Copying the fitscube...'
  file_copy, filename, outdir+outfile
  red_fitscube_newheader, outdir+outfile, hdr
  print, inam + ' : Wrote '+outdir+outfile
  ;; Checksums need to be checked. And updated for the main hdu.
  red_fitscube_checksums, outdir+outfile, trust_datasum = trust_datasum 
  ;; The fitscube files should not be changed after this point!

    ;; Any spectral file to copy?
  if ~keyword_set(no_spectral_file) then begin
    spfile = red_strreplace(infile, '_im.fits', '_sp.fits')
    if ~file_test(indir + spfile) then begin
      print, 'File ',spfile,' does not exist.'
      ans=''
      read, 'Would you like to continue without spectral cube (Y/n) : ',ans
      if strupcase(ans) eq 'N' then return
      print
    endif
    print, inam + ' : Copying the spectral cube...'
    spoutfile = red_strreplace(outfile, '_im.fits', '_sp.fits')
    file_copy, indir+spfile, outdir+spoutfile
    ;; Change some keywords in the spectral file
    sphdr = headfits(outdir+spoutfile)
    ;; Delete header keywords that should not be there
    red_fitsdelkeyword, sphdr, 'STATE'
    ;; Add keywords after this one:
    spanchor = 'TIMESYS'
    red_fitsaddkeyword, anchor = spanchor, sphdr, 'DATE', new_DATE
    red_fitsaddkeyword, anchor = spanchor, sphdr, 'FILENAME', spoutfile
    red_fitscube_newheader, outdir+spoutfile, sphdr
    print, inam + ' : Wrote '+outdir+spoutfile
  endif
  
  ;; Add FITS keywords with info available in the file but not as
  ;; headers. 

  
  if n_elements(smooth_width) gt 0 then begin
    ;; Blur the frames in the fitscube(s)
    naxis = fxpar(hdr, 'NAXIS*')
    Nframes = round(product(naxis[2:*]))

    ;; The statistics will not be right after this!
    
    for iframe = 0, Nframes-1 do begin

      red_progressbar, iframe, Nframes, /predict, 'Blurring frames'

      red_fitscube_getframe, outdir+outfile, frame, iframe = iframe
      red_fitscube_addframe, outdir+outfile, smooth(frame, smooth_width), iframe = iframe

    endfor                      ; iframe
    
  endif

  dimensions = long(fxpar(hdr, 'NAXIS*'))
  Nx      = dimensions[0]
  Ny      = dimensions[1]
  Ntuning = dimensions[2]
  Nstokes = dimensions[3]
  Nscans  = dimensions[4]

  if ~keyword_set(no_thumbnail) then begin
    ;; Get thumbnail from keyword or generate one from the cube?

    if n_elements(iframe) ne 0 then begin

      ;; Get the selected frame
      red_fitscube_getframe, filename, thumb $
                             , iframe = iframe
    endif else begin

      ;; No iframe selected. Use sharpest I image from brightest
      ;; tuning.

      if n_elements(istokes) eq 0 then istokes = 0 ; Stokes I
      
      if n_elements(ituning) eq 0 then begin
        if Ntuning eq 1 then ituning = 0 else begin
          ;; Select the brightest tuning, should be (near-)continuum
          tmp = red_fitsgetkeyword(filename,'DATAMEDN',var=datamedn)
          medn = median(datamedn.values[0, *, istokes, *], dim = 4)
          tmp = max(reform(medn), ituning)
        endelse
      endif
      
      if n_elements(iscan) eq 0 then begin
        if Nscans eq 1 then iscan = 0 else begin
          ;; Select highest RMS contrast, should be the best scan
          tmp = red_fitsgetkeyword(filename,'DATANRMS',var=datanrms)
          tmp = max(reform(datanrms.values[0, ituning, istokes, *]), iscan)
        endelse
      endif

      ;; Get the selected frame
      red_fitscube_getframe, filename, thumb $
                             , iframe = iframe $
                             , iscan = iscan $
                             , istokes = istokes $
                             , ituning = ituning                                
    endelse

    Nxx = Nx - (Nx mod thumb_shrink_fac)
    Nyy = Ny - (Ny mod thumb_shrink_fac)

    thumb = red_centerpic(thumb, xs = Nxx, ys = Nyy)
    thumb = rebin(thumb $
                  , Nxx/thumb_shrink_fac $
                  , Nyy/thumb_shrink_fac $
                  , /samp)
    
    ;; Scale the frame.
    thumb = bytscl(red_histo_opt(thumb))

    ;; Shrink it to thumb size?
    
    ;; Save it as a png image.
    write_png, outdir + red_strreplace(outfile, '.fits', '.png'), thumb
    
  endif

  
  if ~keyword_set(no_video) then begin

    ;; istokes and ituning may be given as keywords or
    ;; automatically selected above. But they might not, so repeat
    ;; those lines here.

    if n_elements(istokes) eq 0 then istokes = 0 ; Stokes I
    
    if n_elements(ituning) eq 0 then begin
      if Ntuning eq 1 then ituning = 0 else begin
        ;; Select the brightest tuning, should be (near-)continuum
        tmp = red_fitsgetkeyword(filename,'DATAMEDN',var=datamedn)
        medn = median(datamedn.values[0, *, istokes, *], dim = 4)
        tmp = max(reform(medn), ituning)
      endelse
    endif

    video_extension = 'mov'
    self -> fitscube_video, filename $
                            , istokes = istokes $
                            , ituning = ituning $
                            , outdir = outdir $
                            , outfile = red_strreplace(outfile, '.fits' $
                                                       , '.'+video_extension) $
                            , shrink_fac = video_shrink_fac
    
  endif
  
  

  
  if n_elements(svo_username) gt 0 and n_elements(svo_api_key) gt 0 then begin

    ;; Ingest metadata into the SVO database.
    
    ;; Look for the python script
    cmd = file_search(strsplit(!path,':',/extract)+'/submit_record.py', count=cnt)
    if cnt eq 0 then begin

      print, inam + " : Didn't find submit_record.py in !path"
      stop
      
    endif else begin

      ;; Keywords to the SVO ingestion script
      file_url = 'https://dubshen.astro.su.se/sst_archive/download/' + outfile
      instrument = strtrim(fxpar(hdr,'INSTRUME'),2)
      
      ;; Generate OID = unique observation ID of the metadata?
      ;; Is date + timestamp enough or should this be a hash-like
      ;; string that changes with versioning etc?
      ;;
      ;; Let's do DATE-BEG_PREFILTER_SCANNUMBERS:
      scn = red_fitsgetkeyword(filename, 'SCANNUM', comment = comment, variable_values = scannum)
      filter = strtrim(fxpar(hdr, 'FILTER1'),2)
      scans = rdx_ints2str(scannum.values)
      oid = date_beg+'T'+time_beg + '_' $
            + filter + '_' $
            + scans

      ;; 'file-path' in SVO is a fake for SST data (not uploaded
      ;; there) but it must be unique.
      ;; Let's make it same as 'oid' but without semicolons
      ;; which are not alloweded.
      tt = strjoin(strsplit(time_beg,':',/extract),'')
      sc = red_strreplace(scans,':','_',n=100)
      file_path = date_beg+'T'+ tt + '_' $
            + filter + '_' + sc
      
      ;; Build the command string to submit metadata to SVO
      cmd = cmd[0]      
      cmd += ' ' + outdir+outfile ; The FITS file to submit to the SVO
      cmd += ' --dataset ' + strupcase(instrument) ; The dataset ID in the SOLARNET Data Archive/SVO
      cmd += ' --file-url ' + file_url                          ; The URL of the file
      cmd += ' --file-path ' + file_path                        ; Fake path of the file
      if n_elements(thumbnail_url) gt 0 then $                  ;
         cmd += ' --thumbnail-url ' + thumbnail_url             ; The URL of the thumbnail
      cmd += ' --username ' + svo_username                      ; The SVO username of the user owning the data
      cmd += ' --api-key ' + svo_api_key                        ; The SVO API key of the user owning the data
      cmd += ' --oid ' + oid                                    ; The unique observation ID of the metadata
      
      ;; Spawn running the script
      spawn, cmd, status
      print, status

      ;; Submit metadata to sst_archive
      ;; Unfortunately we can do it only from dubshen
      spawn, 'whoami',result
      if result eq 'olad6860' then begin
        
        cmd = '/usr/bin/ssh olad6860@dubshen " sudo /root/bin/submit_cube ' + $
              date_beg + '/' + instrument + '/' + outfile 
        if keyword_set(allowed_users) then begin
          cmd += ' ' + allowed_users
          if keyword_set(swedish_data) then cmd += ' --swedish-data'
        endif
        cmd += ' > /dev/null 2>&1" &'
        spawn, cmd
        wait,0.5

        ;; Check submitting status and display it
        spawn, 'tail -1 /home/olexa/submit/status.log', status
        if status ne 'Submitting has been started.' then begin
          print,'We have troubles.'
          return
        endif else begin
          bb = string(13B)      ; CR w/o LF
          outline = string(replicate(32B, 70))
          tw = 0 & ttw = 0
          while status ne 'Submitting has been finished.' do begin
            wait,1
            ttw++ & tw++
            bar = string(replicate(61B, tw))   ; Replicated '='
            bar += string(replicate(45B, 41-tw)) ; Replicated '-'
            strput, outline, string('Submitting: [' + bar + '] ',  ttw, ' sec', format = '(A,I3,A,$)')
            print, bb, outline, FORMAT = '(A,A,$)'
            if tw eq 40 then tw = 0
            print
            spawn, 'tail -1 /home/olexa/submit/status.log', status
          endwhile
        endelse

        ;; Read the log and print it.
        print
        spawn, 'tail -12 /home/olexa/submit/cube_submit.log', result
        print, result
        
      endif else begin
        print
        print,'!!!!!!!!!!!!!!!!!!!!!!!!'
        print, "Unfortunately only admin can submit metadata to the SST archive."
        print, "Please don't forget to ask the administrator to do it."
        print,'!!!!!!!!!!!!!!!!!!!!!!!!'
      endelse      

      ;; We still need to fill information in data_cubes table to use
      ;; idl procedure to download cubes.
      ;; (Current sst_archive doesn't support downloads with wb and sp cubes.)

      ;; Put datacube information to the database
      release_comment = fxpar(hdr,'RELEASEC')
      if fxpar(hdr,'NAXIS4') eq 4 then pol = '1' else pol = '0'
      red_mysql_check, handle
      
      if current_date lt release_date then begin
        filename = outdir+outfile
        query = "INSERT INTO data_cubes (date, time_beg, time_end, wvlnth, pol, scans, filename, release_date,  instrument," + $
                'release_comment, allowed_users) VALUES ("' + date_beg +'", "'+ time_beg+ '", "'+ time_end+'", '+ filter+', '+pol+', "'+ $
                scans+'", "'+outfile +'", "'+ release_date +'", "'+ instrument + '", "' + release_comment + '", "' + allowed_users +'") '+ $
                "ON DUPLICATE KEY UPDATE filename=VALUES(filename), release_date=VALUES(release_date), release_comment=VALUES(release_comment), allowed_users=VALUES(allowed_users);"        
      endif else begin
        query = "INSERT INTO data_cubes (date, time_beg, time_end, wvlnth, pol, scans, filename, release_date,  instrument," + $
                'release_comment) VALUES ("' + date_beg +'", "'+ time_beg +'", "'+ time_end+'", '+ filter+', '+ pol+', "'+ scans+'", "'+ $
                outfile +'", "'+ release_date +'", "'+ instrument + '", "' + release_comment +'") '+ $
                "ON DUPLICATE KEY UPDATE filename=VALUES(filename), release_date=VALUES(release_date), release_comment=VALUES(release_comment);"
      endelse
      
      red_mysql_cmd, handle, query, ans, nl
      free_lun,handle
      
    endelse
  endif else begin

    print, inam + " : Keywords svo_username or svo_api_key not provided."
    print, inam + " : Will not ingest metadata into the SVO database."

  endelse

  
end

;; Code for testing. Test only with smaller cubes because we have to
;; read the entire cube before calling fits_test_checksum below.

cd, '/scratch/mats/2016.09.19/CHROMIS-jan19/'
a = chromisred(/dev)

odir = './export/'
infile = 'cubes_nb/nb_3950_2016-09-19T09:28:36_scans=0-3_corrected_im.fits'

;a -> fitscube_finalize, infile
a -> fitscube_export, infile, outdir = odir


stop

if 0 then begin
  cd, '/scratch/mats/2016.09.19/CRISP-aftersummer/'
  filename = 'cubes_nb/nb_6302_2016-09-19T09:30:20_scans=12-16_stokes_corrected_im.fits'
  a = crispred(/dev)
endif else begin
  cd, '/scratch/mats/2016.09.19/CHROMIS-jan19/'
  filename = 'cubes_nb/nb_3950_2016-09-19T09:28:36_scans=69-75_corrected_im.fits'
  a = chromisred(/dev)
endelse



a -> fitscube_export, filename, outdir = outdir, outfile = outfile $
                      , /smooth $
                      , keywords = { observer:'Some person' $
                                     , requestr:'Boss person' $
                                   } $
                      , release_comment = 'These are test data, never to be released' $
                      , release_date = '2999-09-19'


search_keywords = [ 'STARTOBS' $ 
                    , 'DATE' $
                    , 'DATE-BEG' $
                    , 'DATE-AVG' $
                    , 'DATE-END' $
                    , 'RELEASE' $
                    , 'RELEASEC' $
                    , 'INSTRUME' $
                    , 'TELESCOP' $
                    , 'WCSNAMEA' $
                    , 'CRVAL1A' $
                    , 'CRVAL2A' $
                    , 'NAXIS1' $
                    , 'NAXIS2' $
                    , 'NAXIS3' $
                    , 'NAXIS4' $
                    , 'NAXIS5' $
                    , 'CTYPE1' $
                    , 'CTYPE2' $
                    , 'CTYPE3' $
                    , 'CTYPE4' $
                    , 'CTYPE5' $
                    , 'WAVEMIN' $
                    , 'WAVEMAX' $
                    , 'WAVEUNIT' $
                  ]

h = headfits(outdir+'/'+outfile)
for i = 0, n_elements(search_keywords)-1 do begin
  keyval = fxpar(h, search_keywords[i], count = n)
  if n gt 0 then begin
    print,search_keywords[i], ' : ', keyval
  endif else begin
    print,search_keywords[i], ' : -- '
  endelse
endfor

print, 'DATASUM:  ',fxpar(h, 'DATASUM')
print, 'CHECKSUM: ',fxpar(h, 'CHECKSUM')

d = red_readdata(outdir+'/'+outfile)
test = fits_test_checksum(h, d, ERRMSG = errmsg)
print, 'Test checksums: ', test
if test eq -1 then print, errmsg
end
