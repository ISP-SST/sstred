; docformat = 'rst'

;+
; Copy summed darks, flats, pinholes, polcal from a work directory
; created with the old CRISPRED pipeline.
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
; :Returns:
; 
; 
; :Params:
; 
;   olddir : in, type=string
; 
;     The path to the old work directory from which the summed
;     calibration data should be copied.
; 
; 
; :Keywords:
; 
;   all : in, optional, type=boolean
;
;      Set this to copy all kinds of summed calibration data.
;   
;   darks : in, optional, type=boolean
;
;      Set this to copy summed darks data.
;   
;   flats : in, optional, type=boolean
;
;      Set this to copy summed flats data.
;   
;   pinholes : in, optional, type=boolean
;
;      Set this to copy summed pinhole array data.
; 
;   polcal : in, optional, type=boolean
;
;      Set this to copy summed polcal data.
;   
; 
; 
; :History:
; 
;   2019-05-10 : MGL. First version.
; 
;   2019-08-09 : MGL. Allow variations in the flats directory names.
; 
;-
pro crisp::copy_oldsums, olddir $
                         , all = all $
                         , darks = darks $a
                         , flats = flats $a
                         , overwrite = overwrite $
                         , pinholes = pinholes $
                         , polcal = polcal 
                         
  inam = red_subprogram(/low, calling = inam1)

  print
  print, inam + ' : Attempting to copy summed calibration data from'
  print, olddir

  red_make_prpara, prpara, olddir

  if ~file_test(olddir) then begin
    print
    print, inam + ' : The specified directory does not exist: '
    print, olddir
    return
  end

  if keyword_set(all) then begin
    darks = 1
    flats = 1
    pinholes = 1
    polcal = 1
  endif
  
  ;; Darks
  if keyword_set(darks) then begin

    if file_test(olddir+'/darks') then begin
      
      print
      print, 'Copy dark sums '

      oldfiles = file_search(olddir+'/darks/cam*.dark', count = Nfiles)
      if Nfiles gt 0 then begin

        newdir = self.out_dir+'/darks/'
        file_mkdir, newdir
        
        self -> extractstates, oldfiles, oldstates
        newnames = self -> filenames('dark', oldstates)

        for ifile = 0, Nfiles-1 do begin

          if keyword_set(overwrite) or ~file_test(newnames[ifile]) then begin
            
            data = float(red_readdata(oldfiles[ifile], head = hdr))

            check_fits, data, hdr, /update, /silent
            red_metadata_restore, anchor = 'DATE', hdr

            ;; Remove some keywords
            red_fitsdelkeyword, hdr, 'AO_NMODE'
            red_fitsdelkeyword, hdr, 'OBJECT'
            
            ;; Add/change some keywords
            red_fitsaddkeyword, hdr, 'FILENAME', file_basename(newnames[ifile]), anchor = 'SOLARNET'
            red_fitsaddkeyword, hdr, 'DATE-OBS', self.isodate, anchor = 'DATE' ; The timestamp is not available
            red_fitsaddkeyword, hdr, 'CAMERA', oldstates[ifile].camera, anchor = 'DETECTOR'

            ;; Fake a step that says that the data were summed with
            ;; the old pipeline.
            anchor = 'STATE'
            red_fitsaddkeyword, anchor = anchor, hdr, 'PRSTEP1', 'Summing'      , 'Processing step name'                           
            red_fitsaddkeyword, anchor = anchor, hdr, 'PRPROC1', 'red::sumdark' , 'Name of procedure used'                        
            red_fitsaddkeyword, anchor = anchor, hdr, 'PRLIB1' , 'CRISPRED'     , 'Software library containing red::sumdark'
            red_fitsaddkeyword, anchor = anchor, hdr, 'PRLIB1A', 'IDLAstro'     , 'Additional software library'                    
            red_fitsaddkeyword, anchor = anchor, hdr, 'PRLIB1B', 'Coyote'       , 'Additional software library'                    
            red_fitsaddkeyword, anchor = anchor, hdr, 'PRLIB1C', 'mpfit'        , 'Additional software library'                    
            red_fitsaddkeyword, anchor = anchor, hdr, 'PRLIB1D', 'reduxdlm'     , 'Additional software library'
            ;; This last one needs to be there for headerinfo_addstep
            ;; to know where to put its keywords.
            red_fitsaddkeyword, anchor = anchor, hdr, 'PRBRA1' , ''             , 'Unknown version control branch and versions'

            ;; Then add the copying step
            self -> headerinfo_addstep, hdr, prstep = 'Copying' $
                                        , prproc = inam, prpara = prpara

            ;; Write the data
            red_writedata, newnames[ifile], data, header=hdr, filetype = 'fits', overwrite = overwrite
            
            print, 'Copied ' + file_basename(oldfiles[ifile]) + ' --> ' + file_basename(newnames[ifile])
            
          endif
          
        endfor                  ; ifile
        
      endif else begin
        print
        print, inam + ' : No files in '
        print, olddir
      endelse
    endif else begin
      print
      print, inam + ' : No "darks" subdirectory in '
      print, olddir
    endelse
  endif

  
  ;; Flats
  if keyword_set(flats) then begin
    
    ;; Flats subdirectories are allowed to be the standard name plus
    ;; some extension. Copy all such subdirs.
    fdirs = file_search(olddir+'/flats*', count = Nfdirs)
    for idir = 0, Nfdirs-1 do begin
;   if file_test(olddir+'/flats') then begin

      print
      print, 'Copy flat sums from '+fdirs[idir]

      oldfiles = file_search(fdirs[idir]+'/cam*.flat', count = Nfiles)

      if Nfiles gt 0 then begin

        ;; Don't copy unpol flats, make them with the new pipeline
        indx = where(~strmatch(oldfiles, '*unpol*'), Nfiles)
        if Nfiles gt 0 then oldfiles = oldfiles[indx]
        
      endif
        
      if Nfiles gt 0 then begin

        newdir = self.out_dir+'/'+file_basename(fdirs[idir])+'/'
        file_mkdir, newdir
        
        self -> extractstates, oldfiles, oldstates

        ;; WB flats need some extra care
        indx=where(oldstates.camera eq 'Crisp-W', Nwb)
        for iwb = 0, Nwb-1 do begin
          splitname = strsplit(file_basename(oldfiles[indx[iwb]]), '.', /extract)
          oldstates[indx[iwb]].prefilter = splitname[1]
          ;; For the fullstate, fake the "line" info the make the
          ;; tuning part complete.
          oldstates[indx[iwb]].fullstate = strjoin([oldstates[indx[iwb]].prefilter $
                                                    , oldstates[indx[iwb]].prefilter $
                                                    , '+0000' $
                                                   ], '_')
        endfor                  ; iwb

        newnames = self -> filenames('flat', oldstates)
        ;; In case we are dealing with multiple flats dirs:
        newnames = red_strreplace(newnames,'/flats/','/'+file_basename(fdirs[idir])+'/')

        for ifile = 0, Nfiles-1 do begin

          if keyword_set(overwrite) or ~file_test(newnames[ifile]) then begin
            
            ;; Tweak states some more
            oldstates[ifile].pf_wavelength = float(oldstates[ifile].prefilter)*1e-10
            
            data = float(red_readdata(oldfiles[ifile], head = hdr))
            
            check_fits, data, hdr, /update, /silent

            red_metadata_restore, anchor = 'DATE', hdr

            ;; Remove some keywords
            red_fitsdelkeyword, hdr, 'AO_NMODE'

            ;; Add/change some keywords
            if oldstates[ifile].is_wb then begin
              anchor = 'DETECTOR'
              red_fitsaddkeyword, anchor = anchor, hdr, 'FILTER1', oldstates[ifile].prefilter
              red_fitsaddkeyword, anchor = anchor, hdr, 'STATE', oldstates[ifile].fullstate
              red_fitsaddkeyword, anchor = anchor, hdr, 'WAVEUNIT',  -9, 'WAVELNTH in units 10^WAVEUNIT m = nm'
            endif
            red_fitsaddkeyword, hdr, 'FILENAME', file_basename(newnames[ifile]), anchor = 'SOLARNET'
            red_fitsaddkeyword, hdr, 'DATE-OBS', self.isodate, anchor = 'DATE' ; The timestamp is not available
            red_fitsaddkeyword, hdr, 'CAMERA',   oldstates[ifile].camera, anchor = 'DETECTOR'

            ;; Subset of headerinfo_addstep keywords that say the data
            ;; were summed with the old pipeline.
            anchor = 'STATE'
            red_fitsaddkeyword, anchor = anchor, hdr, 'PRSTEP1', 'Summing'      , 'Processing step name'                           
            red_fitsaddkeyword, anchor = anchor, hdr, 'PRPROC1', 'red::sumflat' , 'Name of procedure used'                        
            red_fitsaddkeyword, anchor = anchor, hdr, 'PRLIB1' , 'CRISPRED'     , 'Software library containing red::sumflat'
            red_fitsaddkeyword, anchor = anchor, hdr, 'PRLIB1A', 'IDLAstro'     , 'Additional software library'                    
            red_fitsaddkeyword, anchor = anchor, hdr, 'PRLIB1B', 'Coyote'       , 'Additional software library'                    
            red_fitsaddkeyword, anchor = anchor, hdr, 'PRLIB1C', 'mpfit'        , 'Additional software library'                    
            red_fitsaddkeyword, anchor = anchor, hdr, 'PRLIB1D', 'reduxdlm'     , 'Additional software library'                    
            ;; This last one needs to be there for the next call of
            ;; headerinfo_addstep to know where to put its keywords:
            red_fitsaddkeyword, anchor = anchor, hdr, 'PRBRA1' , ''             , 'Unknown version control branch and versions'

            ;; Then add this copying step
            self -> headerinfo_addstep, hdr, prstep = 'Copying' $
                                        , prproc = inam, prpara = prpara

            ;; Write the data
            red_writedata, newnames[ifile], data, header=hdr, filetype = 'fits', overwrite = overwrite

            print, 'Copied ' + file_basename(oldfiles[ifile]) + ' --> ' + file_basename(newnames[ifile])
            
          endif
          
        endfor                  ; ifile
        
      endif else begin
        print
        print, inam + ' : No files in '
        print, olddir
      endelse

      
    endfor                      ; idir

    if Nfdirs eq 0 then begin
      print
      print, inam + ' : No "flats" subdirectory in'
      print, olddir
    endif
  endif
    

  ;; Pinholes
  if keyword_set(pinholes) then begin
    if file_test(olddir+'/pinh_align') then begin
      
      print
      print, 'Copy pinhole sums '

      ;; The old pipeline summed LC states separately, while the new
      ;; pipeline sums them together and saves in file names without
      ;; LC state tags. Right now all LC states are copied to the same
      ;; name, making LC3 overwrite LC0-2. This is probably not a
      ;; problem but it would be better to sum the images for
      ;; different LC states (divided by 4).
      
      ;; Search for floating point pinholes
      oldfiles = file_search(olddir+'/pinh_align/cam*.fpinh', count = Nfiles)

      if Nfiles gt 0 then begin ; Fall back on integer pinholes
        oldfiles = file_search(olddir+'/pinh_align/cam*.pinh', count = Nfiles)
      endif
      
      if Nfiles gt 0 then begin

        newdir = self.out_dir+'/pinhs/'
        file_mkdir, newdir
        
        self -> extractstates, oldfiles, oldstates
        
        newnames = self -> filenames('pinh', oldstates)

        for ifile = 0, Nfiles-1 do begin

          if keyword_set(overwrite) or ~file_test(newnames[ifile]) then begin
            
            ;; Tweak states some more
            oldstates[ifile].pf_wavelength = float(oldstates[ifile].prefilter)*1e-10
            
            data = float(red_readdata(oldfiles[ifile], head = hdr))
            
            check_fits, data, hdr, /update, /silent

            red_metadata_restore, anchor = 'DATE', hdr

            ;; Add/change some keywords
            red_fitsaddkeyword, hdr, 'FILENAME', file_basename(newnames[ifile]), anchor = 'SOLARNET'
            red_fitsaddkeyword, hdr, 'DATE-OBS', self.isodate, anchor = 'DATE' ; The timestamp is not available
            red_fitsaddkeyword, hdr, 'CAMERA',   oldstates[ifile].camera, anchor = 'DETECTOR'

            ;; Subset of headerinfo_addstep keywords that say the data
            ;; were summed with the old pipeline.
            anchor = 'STATE'
            red_fitsaddkeyword, anchor = anchor, hdr, 'PRSTEP1', 'Summing'      , 'Processing step name'                           
            red_fitsaddkeyword, anchor = anchor, hdr, 'PRPROC1', 'red::sumflat' , 'Name of procedure used'                        
            red_fitsaddkeyword, anchor = anchor, hdr, 'PRLIB1' , 'CRISPRED'     , 'Software library containing red::sumflat'
            red_fitsaddkeyword, anchor = anchor, hdr, 'PRLIB1A', 'IDLAstro'     , 'Additional software library'                    
            red_fitsaddkeyword, anchor = anchor, hdr, 'PRLIB1B', 'Coyote'       , 'Additional software library'                    
            red_fitsaddkeyword, anchor = anchor, hdr, 'PRLIB1C', 'mpfit'        , 'Additional software library'                    
            red_fitsaddkeyword, anchor = anchor, hdr, 'PRLIB1D', 'reduxdlm'     , 'Additional software library'                    
            ;; This last one needs to be there for the next call of
            ;; headerinfo_addstep to know where to put its keywords:
            red_fitsaddkeyword, anchor = anchor, hdr, 'PRBRA1' , ''             , 'Unknown version control branch and versions'

            ;; Then add this copying step
            self -> headerinfo_addstep, hdr, prstep = 'Copying' $
                                        , prproc = inam, prpara = prpara

            ;; Write the data
            red_writedata, newnames[ifile], data, header=hdr, filetype = 'fits', overwrite = overwrite

            print, 'Copied ' + file_basename(oldfiles[ifile]) + ' --> ' + file_basename(newnames[ifile])
            
          endif
          
        endfor                  ; ifile
      endif 
    endif  else begin
      print
      print, inam + ' : No "pinh_align" subdirectory in'
      print, olddir
    endelse
  endif
  
  
  ;; Polcal
  if keyword_set(polcal) then begin
    if file_test(olddir+'/polcal_sums') then begin

      cams = 'Crisp-'+['T', 'R']
      
      for icam = 0, 1 do begin

        print
        print, 'Copy polcal sums for ' + cams[icam]

        oldfiles = file_search(olddir+'/polcal_sums/'+cams[icam]+'/cam*', count = Nfiles)

        if Nfiles gt 0 then begin
          
          newdir = self.out_dir+'/polcal_sums/' + cams[icam]
          file_mkdir, newdir
          
          self -> extractstates, oldfiles, oldstates, /polcal

          newnames = self -> filenames('pols', oldstates)

          for ifile = 0, Nfiles-1 do begin

            if keyword_set(overwrite) or ~file_test(newnames[ifile]) then begin
              
              ;; Tweak states some more
;              oldstates[ifile].pf_wavelength = float(oldstates[ifile].prefilter)*1e-10
              
              data = float(red_readdata(oldfiles[ifile], head = hdr))
              
              check_fits, data, hdr, /update, /silent

              red_metadata_restore, anchor = 'DATE', hdr

              ;; Remove some keywords
              red_fitsdelkeyword, hdr, 'AO_NMODE'

              ;; Add/change some keywords
              red_fitsaddkeyword, hdr, 'FILENAME', file_basename(newnames[ifile]), anchor = 'SOLARNET'
              red_fitsaddkeyword, hdr, 'DATE-OBS', self.isodate, anchor = 'DATE' ; The timestamp is not available

              anchor = 'CAMERA'
              red_fitsaddkeyword, anchor = anchor, hdr, 'CALIB_QW' $
                                  , oldstates[ifile].qw, '[deg] Quarterwave plate angle'
              red_fitsaddkeyword, anchor = anchor, hdr, 'CALIB_LP' $
                                  , oldstates[ifile].lp, '[deg] Linear polarizer angle'
              red_fitsaddkeyword, anchor = anchor, hdr, 'POL_LC' $
                                  , oldstates[ifile].lc, 'Liquid crystal state number'
              red_fitsaddkeyword, anchor = anchor, hdr, 'STATE' $
                                  , oldstates[ifile].fullstate, 'Polcal state'
              
              ;; Subset of headerinfo_addstep keywords that say the data
              ;; were summed with the old pipeline.
              anchor = 'STATE'
              red_fitsaddkeyword, anchor = anchor, hdr, 'PRSTEP1', 'Summing'      , 'Processing step name'                           
              red_fitsaddkeyword, anchor = anchor, hdr, 'PRPROC1', 'red::sumflat' , 'Name of procedure used'                        
              red_fitsaddkeyword, anchor = anchor, hdr, 'PRLIB1' , 'CRISPRED'     , 'Software library containing red::sumflat'
              red_fitsaddkeyword, anchor = anchor, hdr, 'PRLIB1A', 'IDLAstro'     , 'Additional software library'                    
              red_fitsaddkeyword, anchor = anchor, hdr, 'PRLIB1B', 'Coyote'       , 'Additional software library'                    
              red_fitsaddkeyword, anchor = anchor, hdr, 'PRLIB1C', 'mpfit'        , 'Additional software library'                    
              red_fitsaddkeyword, anchor = anchor, hdr, 'PRLIB1D', 'reduxdlm'     , 'Additional software library'                    
              ;; This last one needs to be there for the next call of
              ;; headerinfo_addstep to know where to put its keywords:
              red_fitsaddkeyword, anchor = anchor, hdr, 'PRBRA1' , ''             , 'Unknown version control branch and versions'

              ;; Then add this copying step
              self -> headerinfo_addstep, hdr, prstep = 'Copying' $
                                          , prproc = inam, prpara = prpara

              ;; Write the data
              red_writedata, newnames[ifile], data, header=hdr, filetype = 'fits', overwrite = overwrite

              print, 'Copied ' + file_basename(oldfiles[ifile]) + ' --> ' + file_basename(newnames[ifile])
              
            endif
            
          endfor                ; ifile
        endif

      endfor                    ; icam
      
    endif else begin
      print, inam + ' : No "polcal_sums" subdirectory in'
      print, olddir
    endelse
  endif
      
      

end

cd, '/scratch/mats/2016.09.19/CRISP-copysums/'
a = crispred(/dev)
a -> copy_oldsums, '/scratch/mats/2016.09.19/CRISP', /over, /all
;a -> copy_oldsums, '/scratch/mats/2016.09.19/CRISP', /over, /pol
;a -> copy_oldsums, '/scratch/mats/2016.09.19/CRISP', /over, /pin
;a -> copy_oldsums, '/scratch/mats/2016.09.19/CRISP', /over, /flat
;a -> copy_oldsums, '/scratch/mats/2016.09.19/CRISP', /over, /dark

end
