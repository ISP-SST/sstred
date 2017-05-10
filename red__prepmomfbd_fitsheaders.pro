; docformat = 'rst'

;+
; Write FITS header files for expected MOMFBD output.
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, ISP
; 
; 
; 
; :Keywords:
;
;    no_pd : in, optional, type=boolean
;   
;       Set this to exclude phase diversity data and processing. 
;   
;    momfbddir : in, optional, type=string, default='momfbd'
;   
;       Top directory of output tree.
;   
; 
;    dirs : in, optional, type=strarr
;   
;       The data "timestamp" directories to process.
;   
;    pref : in, optional, type=string
;
;       Prefilter. 
; 
; :History:
; 
;    2017-03-16 : MGL. First version.
; 
;    2017-03-24 : MGL. Add header information about this processing step.
; 
; 
;-
pro red::prepmomfbd_fitsheaders, dirs = dirs $
                                 , momfbddir = momfbddir $
                                 , pref = pref $
                                 , scanno = scanno $
                                 , no_pd = no_pd 
  
  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

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
  if Ndirs eq 0 then return
 
  if n_elements(momfbddir) eq 0 then begin
    if keyword_set(no_pd) then begin
      momfbddir = 'momfbd_nopd' 
    endif else begin
      momfbddir = 'momfbd' 
    endelse
  endif

  ;; Loop directories
  for idir=0L, Ndirs-1 do begin
    
    dir = dirs[idir]+'/'
    folder_tag = file_basename(dir)
    datestamp = self.isodate+'T'+folder_tag
    
    ;; base output location
    cfg_base_dir = self.out_dir + PATH_SEP() + momfbddir + PATH_SEP() + folder_tag + PATH_SEP()
    
    prefs = file_basename(file_search(cfg_base_dir+'*', count = Nprefs))

    for ipref = 0, Nprefs-1 do begin
      
      cfg_dir = cfg_base_dir + PATH_SEP() + prefs[ipref] + PATH_SEP() + 'cfg/'
      cfg_files = file_search(cfg_dir + 'momfbd_reduc_' + prefs[ipref] $
                              + '_?????.cfg', count = Ncfg)
      
      progress_msg = 'Making fits headers for ' $
                     + red_strreplace(red_strreplace(cfg_dir,'//','/') $
                                      , self.out_dir,'')

      ;; Parse all config files, make fitsheaders for all output files
      for icfg = 0, Ncfg-1 do begin

        red_progressbar, icfg, Ncfg, progress_msg, clock=clock, /predict

    
        if n_elements(scanno) ne 0 then $
           if long(scanno) ne long((strsplit(file_basename(cfg_files[icfg],'.cfg') $
                                             ,'_',/extract))[3]) then $
                                                continue
        

        spawn, 'cat '+cfg_files[icfg], cfg
        Nlines = n_elements(cfg)

        ;; Get prpara (processing parameters) from the global config keywords
        Gstart = (where(strmatch(cfg,'}*'),Nmatch))[Nmatch-1] + 1
        undefine, Gprpara       ; Start fresh
        for iline = Gstart, Nlines-1 do begin
          cfgsplit = strsplit(cfg[iline], '=', /extract)
          
          case cfgsplit[0] of
            ;; Make list of possible global keywords complete!
            'FILE_TYPE' : begin
              file_type = cfgsplit[1]
              red_make_prpara, Gprpara, cfgsplit[0], cfgsplit[1]
            end
            'IMAGE_NUMS'      : begin
              ;; These are actually file numbers, several frames in each file.
              file_nums = red_expandrange(cfgsplit[1]) 
              red_make_prpara, Gprpara, cfgsplit[0], file_nums
            end
            'MODES' : red_make_prpara, Gprpara, cfgsplit[0], red_expandrange(cfgsplit[1])
            'SIM_X' : red_make_prpara, Gprpara, cfgsplit[0], red_expandrange(cfgsplit[1])
            'SIM_Y' : red_make_prpara, Gprpara, cfgsplit[0], red_expandrange(cfgsplit[1])
            ;; Ignored config lines:
            'PROG_DATA_DIR' : 
            else : begin
              case n_elements(cfgsplit) of
                1: red_make_prpara, Gprpara, cfg[iline]
                2: red_make_prpara, Gprpara, cfgsplit[0], cfgsplit[1]
                else: stop
              endcase
            end
          endcase
        endfor                  ; iline
        
        
        ;; Parse objects
        Ostarts = where(strmatch(cfg,'object{*'),Nobj)
        for iobj = 0, Nobj-1 do begin

          prpara = Gprpara      ; Start from the global parameters

          Ostart = Ostarts[iobj]
          Oend   = Ostart + (where(strmatch(cfg[Ostart:*],'}*'),Nmatch))[0]
          Cstarts = Ostart + where(strmatch(cfg[Ostart:Oend],'*channel{*'),Nchan)

          ;; Parse object keywords
          for iline = Ostart, Cstarts[0]-1 do begin
            cfgline = strtrim(cfg[iline], 2)
            if ~strmatch(cfgline,'*{') and ~strmatch(cfgline,'*}') then begin
              cfgsplit = strsplit(cfgline, '=', /extract)
              case cfgsplit[0] of
                'OUTPUT_FILE' : begin
                  output_file = cfgsplit[1]
                  red_make_prpara, prpara, cfgsplit[0], cfgsplit[1]
                end
                else : begin
                  case n_elements(cfgsplit) of
                    1: red_make_prpara, prpara, cfg[iline]
                    2: red_make_prpara, prpara, cfgsplit[0], cfgsplit[1]
                    else: stop
                  endcase
                end
              endcase
            endif
          endfor ; iline

          if n_elements(file_type) eq 0 then file_type = 'ANA' ; Default in momfbd program.
          header_file = cfg_dir + output_file + '.fitsheader'
          case file_type of
            'MOMFBD' : output_file += '.momfbd'
            'ANA'    : output_file += '.f0'
            'FITS'   : output_file += '.fits'
            else: stop
          endcase

          ;; Channels
          if Nchan gt 0 then prmode = 'Phase Diversity' else prmode = ''
          for ichan = 0, Nchan-1 do begin
            
            Cstart = Cstarts[ichan]
            Cend   = Cstart + (where(strmatch(cfg[Cstart:*],'*}*'),Nmatch))[0]

            undefine, discard
            for iline = Cstart, Cend do begin
              cfgline = strtrim(cfg[iline], 2)
              if ~strmatch(cfgline,'*{') and ~strmatch(cfgline,'*}') then begin
                cfgsplit = strsplit(cfgline, '=', /extract)
                case cfgsplit[0] of
                  'ALIGN_CLIP' : red_make_prpara, prpara, cfgsplit[0] $
                                                  , red_expandrange(cfgsplit[1])
                  'FILENAME_TEMPLATE' : begin
                    filename_template = cfgsplit[1]
                    red_make_prpara, prpara, cfgsplit[0], cfgsplit[1]
                  end
                  'GAIN_FILE'         : begin
                    gain_file = cfgsplit[1]
                    red_make_prpara, prpara, cfgsplit[0], cfgsplit[1] 
                  end
                  'IMAGE_DATA_DIR'    : begin
                    image_data_dir = cfgsplit[1]
                    red_make_prpara, prpara, cfgsplit[0], cfgsplit[1]
                  end
                  'DISCARD'           : begin
                    discard = long(cfgsplit[1])
                    red_make_prpara, prpara, cfgsplit[0], cfgsplit[1]
                  end
                  else : begin
                    case n_elements(cfgsplit) of
                      1: red_make_prpara, prpara, cfg[iline]
                      2: red_make_prpara, prpara, cfgsplit[0], cfgsplit[1]
                      else: stop
                    endcase
                  end
                endcase
              endif
            endfor              ; iline

            if ichan eq 0 then begin
              ;; Make a list of files actually involved in the sum.
              undefine, fnames
              for ifile = 0, n_elements(file_nums)-1 do begin
                fname = image_data_dir + '/' $
                        + string(file_nums[ifile], format='(%"'+filename_template+'")')
                if file_test(fname) then red_append, fnames, fname
              endfor            ; ifile
              Nfiles = n_elements(fnames)
              if Nfiles eq 0 then stop
              
            endif

          endfor                ; ichan

          if n_elements(fnames) eq 0 then begin
            print, inam+' : WARNING, no files for '+filename_template
            continue
          endif

          ;; Make header corresponding to the sum.
          head = red_sumheaders(fnames, discard = discard)

          self -> headerinfo_addstep, head $
                                      , prstep = 'MOMFBD image restoration' $
                                      , prpara = prpara $
                                      , prmode = prmode $
                                      , addlib = 'momfbd/redux' 
          ;; At this point we don't know which program will be used,
          ;; momfbd or redux. Select one (if possible) and add version
          ;; number when reading the output.

          ;; Additional keywords that should be set after momfbd
          ;; processing.
          fxaddpar, head, 'FILLED', 1, 'Missing pixels have been filled.' ; Check gain file for missing?
          fxaddpar, head, 'FILENAME', file_basename(output_file), 'MOMFBD restored data'
          
          fxaddpar, head, before='DATE', 'SOLARNET', 0.5
          fxaddpar, head, before='DATE', 'BTYPE', 'Intensity'
          fxaddpar, head, before='DATE', 'BUNIT', 'DU' ; Digital unit?
          ;; DATE_OBS should be getting the value including decimals from
          ;; the raw data headers, not just integer seconds as here:
          fxaddpar, head, before='DATE', 'DATE_OBS', datestamp
          fxaddpar, head, before='DATE', 'STARTOBS', datestamp ; IS STARTOBS needed?

          ;; The CDELTn keywords should not change to HPLN-TAN/HPLT-TAN
          ;; until we know the position and orientation? /MGL
          fxaddpar, head, before='DATE', 'CTYPE1', 'x',                     '[arcsec]'
          fxaddpar, head, before='DATE', 'CTYPE2', 'y',                     '[arcsec]'
          fxaddpar, head, before='DATE', 'CDELT1', float(self.image_scale), '[arcsec] x-coordinate increment'
          fxaddpar, head, before='DATE', 'CDELT2', float(self.image_scale), '[arcsec] y-coordinate increment'
          fxaddpar, head, before='DATE', 'CUNIT1', 'arcsec', 'Unit along axix 1'
          fxaddpar, head, before='DATE', 'CUNIT2', 'arcsec', 'Unit along axix 2'

          ;; Write the header file
          fxwrite, header_file, head 

          
        endfor                  ; iobj
        
      endfor                    ; icfg

    endfor                      ; ipref

  endfor                        ; idir

end
