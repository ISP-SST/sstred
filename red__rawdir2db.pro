; docformat = 'rst'

;+
; Write metadata for a raw dataset into the database.
;
; The dataset is defined by its timestamp directory.
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
;   timestamp, in, type=string
; 
;     The timestamp-directory to process.
; 
; 
; :Keywords:
; 
;   cams : in, optional, type=strarr, default=all
;   
;     The cameras to include.
; 
; 
; :History:
; 
;   2018-07-02 : MGL. First version. 
; 
;-
pro red::rawdir2db, timestamp, cams = cams

  inam = red_subprogram(/low, calling = inam1)

  if n_elements(cams) eq 0 then begin
    cams = *self.cameras
  endif
  Ncams = n_elements(cams)
  
  for icam = 0, Ncams-1 do begin
    
    dir = self.root_dir + '/' + 'CHROMIS-data/' + timestamp + '/' + cams[icam] + '/'
    
    files = file_search(dir+'*', count = Nfiles)
    
    if Nfiles eq 0 then continue


    ;files = files[0:10]         ; testing...
    Nfiles = n_elements(files)        

    self -> extractstates, files, states



    
    ;; Array with NB prefilters to be indexed with the "wheel" number
    ;; minus 1. Same index to use for lambda_ref, du_ref, convfac.
    nbprefs = ['3925', '3934', '3969', '3978', '3999', '4862']
    Nnbprefs = n_elements(nbprefs)
    lambda_ref = dblarr(Nnbprefs)
    du_ref     = dblarr(Nnbprefs)
    convfac    = dblarr(Nnbprefs)
    for wheel = 1, Nnbprefs do begin
      zfile = self.out_dir + 'info/hrz_zeropoint_' + nbprefs[wheel-1] + '.fz'
      refinfo = f0(zfile)
      lambda_ref[wheel-1] = refinfo[0]
      du_ref[wheel-1]     = refinfo[1]
      convfac[wheel-1]    = refinfo[2]
    endfor                      ; wheel

    
    ;; Make struct for info from header and states, that should (or
    ;; at least might) go into the database.

    strct  = { $
             ;; From the states struct
             CAMERA         : '' $ 
             , CAM_SETTINGS   : '' $
             , DETECTOR       : '' $
             , EXPOSURE       : 0d $
             , FPI_STATE      : '' $
             , FRAMENUMBERS    : lonarr(100) $
             , FULLSTATE      : '' $
             , GAIN           : 0. $
             , NFRAMES        : 0L $
             , PF_WAVELENGTH  : 0d $
             , PREFILTER      : '' $
             , SCANNUMBER     : 0L $
             , SKIP           : 0B $
             , TUNING         : '' $
             , TUN_WAVELENGTH : 0d $
             , FILENAME       : '' $
             , IS_WB          : 0B $
             ;; From the header
             , BITPIX   : 0B $  
             , CADAVG   : 0. $
             , DATE     : '' $
             , DATE_BEGS : strarr(100) $
             , DATE_OBS : '' $
             , DETFIRM  : '' $
             , DETGAIN  : 0. $
             , DETOFFS  : 0L $
             , FILTER1  : '' $
             , FRAMENUM : 0L $
             , INSTRUME : '' $
             , NAXIS    : 0B $
             , NAXIS1   : 0L $          
             , NAXIS2   : 0L $
             , NAXIS3   : 0L $
             , OBJECT   : '' $
             , OBSERVER : '' $
             , ORIGIN   : '' $
             , SCANNUM  : 0L $
             , TELESCOP : '' $
             , WAVEBAND : '' $
             , WAVELNTH : 0. $
             , WAVEMAX  : 0. $
             , WAVEMIN  : 0. $
             , WAVEUNIT : 0  $
             , XPOSURE  : 0. $
             }

    ;; Make array of such structs
    dbinfo = replicate(strct, Nfiles)

    ;; Fill in values in the array of structs
    for ifile = 0, Nfiles-1 do begin

      red_progressbar, ifile, Nfiles, /predict, 'Get DB values for '+cams[icam]
      
      h = red_readhead(files[ifile],date_beg=date_beg,framenumbers=framenumbers)
      s = states[ifile]

      ;; For WB and PD, need to translate FPI_STATE from something
      ;; like "wheel00005_hrz32061" to the proper tuning on the form
      ;; "4000_+1258" (what we get in the NB states).
      if s.is_wb then begin

        fpi_info = stregex(s.fpi_state, 'wheel([0-9]*)_hrz([0-9]*)', /extract, /sub)

        ;; Determine NB prefilter
        wheel = long(fpi_info[1])
        nbpref = nbprefs[wheel]
        
        ;; Tuning in digital units
        du = fpi_info[2]

        ;; Tuning in [m]
        dlambda = convfac[wheel-1] * (du-du_ref[wheel-1]) 
        lambda_ref_string = string(round(lambda_ref[wheel-1]*1d10), format = '(i04)')
        
        tuning_string = strtrim(round(dlambda*1d13), 2)
        if strmid(tuning_string, 0, 1) ne '-' then tuning_string = '+'+tuning_string
        s.fpi_state = lambda_ref_string + '_' + tuning_string

      endif

      ;; From the states struct
      if n_elements(s.FILENAME      ) then dbinfo[ifile].FILENAME       = s.FILENAME        
      if n_elements(s.DETECTOR      ) then dbinfo[ifile].DETECTOR       = s.DETECTOR        
      if n_elements(s.CAMERA        ) then dbinfo[ifile].CAMERA         = s.CAMERA          
      if n_elements(s.FULLSTATE     ) then dbinfo[ifile].FULLSTATE      = s.FULLSTATE       
      if n_elements(s.FPI_STATE     ) then dbinfo[ifile].FPI_STATE      = s.FPI_STATE       
      if n_elements(s.NFRAMES       ) then dbinfo[ifile].NFRAMES        = s.NFRAMES         
      if n_elements(s.SKIP          ) then dbinfo[ifile].SKIP           = s.SKIP            
      if n_elements(s.SCANNUMBER    ) then dbinfo[ifile].SCANNUMBER     = s.SCANNUMBER      
      if n_elements(s.TUNING        ) then dbinfo[ifile].TUNING         = s.TUNING          
      if n_elements(s.PREFILTER     ) then dbinfo[ifile].PREFILTER      = s.PREFILTER       
      if n_elements(s.PF_WAVELENGTH ) then dbinfo[ifile].PF_WAVELENGTH  = s.PF_WAVELENGTH   
      if n_elements(s.TUN_WAVELENGTH) then dbinfo[ifile].TUN_WAVELENGTH = s.TUN_WAVELENGTH  
      if n_elements(s.EXPOSURE      ) then dbinfo[ifile].EXPOSURE       = s.EXPOSURE        
      if n_elements(s.GAIN          ) then dbinfo[ifile].GAIN           = s.GAIN            
      if n_elements(s.CAM_SETTINGS  ) then dbinfo[ifile].CAM_SETTINGS   = s.CAM_SETTINGS    
      if n_elements(s.IS_WB         ) then dbinfo[ifile].IS_WB          = s.IS_WB           
      ;; From the header
      value = fxpar(h, 'BITPIX'  , count=cnt) & if cnt eq 1 then dbinfo[ifile].BITPIX   = value 
      value = fxpar(h, 'NAXIS'   , count=cnt) & if cnt eq 1 then dbinfo[ifile].NAXIS    = value 
      value = fxpar(h, 'NAXIS1'  , count=cnt) & if cnt eq 1 then dbinfo[ifile].NAXIS1   = value           
      value = fxpar(h, 'NAXIS2'  , count=cnt) & if cnt eq 1 then dbinfo[ifile].NAXIS2   = value 
      value = fxpar(h, 'ORIGIN'  , count=cnt) & if cnt eq 1 then dbinfo[ifile].ORIGIN   = value 
      value = fxpar(h, 'TELESCOP', count=cnt) & if cnt eq 1 then dbinfo[ifile].TELESCOP = value 
      value = fxpar(h, 'INSTRUME', count=cnt) & if cnt eq 1 then dbinfo[ifile].INSTRUME = value 
      value = fxpar(h, 'DATE-OBS', count=cnt) & if cnt eq 1 then dbinfo[ifile].DATE_OBS = value 
      value = fxpar(h, 'DATE'    , count=cnt) & if cnt eq 1 then dbinfo[ifile].DATE     = value 
      value = fxpar(h, 'WAVELNTH', count=cnt) & if cnt eq 1 then dbinfo[ifile].WAVELNTH = value 
      value = fxpar(h, 'WAVEMIN' , count=cnt) & if cnt eq 1 then dbinfo[ifile].WAVEMIN  = value 
      value = fxpar(h, 'WAVEMAX' , count=cnt) & if cnt eq 1 then dbinfo[ifile].WAVEMAX  = value 
      value = fxpar(h, 'WAVEBAND', count=cnt) & if cnt eq 1 then dbinfo[ifile].WAVEBAND = value 
      value = fxpar(h, 'WAVEUNIT', count=cnt) & if cnt eq 1 then dbinfo[ifile].WAVEUNIT = value 
      value = fxpar(h, 'OBSERVER', count=cnt) & if cnt eq 1 then dbinfo[ifile].OBSERVER = value 
      value = fxpar(h, 'OBJECT'  , count=cnt) & if cnt eq 1 then dbinfo[ifile].OBJECT   = value 
      value = fxpar(h, 'XPOSURE' , count=cnt) & if cnt eq 1 then dbinfo[ifile].XPOSURE  = value 
      value = fxpar(h, 'CADENCE' , count=cnt) & if cnt eq 1 then dbinfo[ifile].CADAVG   = value 
      value = fxpar(h, 'DETGAIN' , count=cnt) & if cnt eq 1 then dbinfo[ifile].DETGAIN  = value 
      value = fxpar(h, 'DETOFFS' , count=cnt) & if cnt eq 1 then dbinfo[ifile].DETOFFS  = value 
      value = fxpar(h, 'DETFIRM' , count=cnt) & if cnt eq 1 then dbinfo[ifile].DETFIRM  = value 
      value = fxpar(h, 'SCANNUM' , count=cnt) & if cnt eq 1 then dbinfo[ifile].SCANNUM  = value 
      value = fxpar(h, 'FRAMENUM', count=cnt) & if cnt eq 1 then dbinfo[ifile].FRAMENUM = value 
      value = fxpar(h, 'FILTER1' , count=cnt) & if cnt eq 1 then dbinfo[ifile].FILTER1  = value        
      
      ;; From the red_readhead call
      if s.NFRAMES eq 0 or s.NFRAMES gt n_elements( dbinfo[ifile].DATE_BEGS) then begin
        print, inam+' : Please increase the length of the following two arrays to at least '+strtrim(s.NFRAMES)
        retall
      end
      if n_elements(date_beg) ne 0 then dbinfo[ifile].DATE_BEGS = date_beg
      if n_elements(framenumbers) ne 0 then dbinfo[ifile].FRAMENUMBERS = framenumbers

    endfor                      ; ifile
    
    ;; Send the array of structs as input to a command that knows what
    ;; info should go into what table, and writes it there.
    ;;red_rawfile2db, dbinfo
    
  endfor                        ; icam

  stop
end



;; a = chromisred(/dev)

;; a -> rawdir2db, '09:28:36'

;end                             
