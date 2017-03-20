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
; 
;-
pro red::prepmomfbd_fitsheaders, dirs = dirs $
                                 , momfbddir = momfbddir $
                                 , pref = pref $
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
      momfbddir = 'momfbd' 
    endif else begin
      momfbddir = 'momfbd_pd' 
    endelse
  endif

  ;; Loop directories
  for idir=0L, Ndirs-1 do begin
    
    dir = dirs[idir]+'/'
    folder_tag = file_basename(dir)
    
    ;; base output location
    cfg_base_dir = self.out_dir + PATH_SEP() + momfbddir + PATH_SEP() + folder_tag + PATH_SEP()
    
    prefs = file_basename(file_search(cfg_base_dir+'*', count = Nprefs))
    
    for ipref = 0, Nprefs-1 do begin
      
      cfg_dir = cfg_base_dir + PATH_SEP() + prefs[ipref] + PATH_SEP() + 'cfg/'
      cfg_files = file_search(cfg_dir+'momfbd_reduc_'+prefs[ipref]+'_?????.cfg', count = Ncfg)
      
      progress_msg = 'Making fitsheaders for ' + strreplace(strreplace(cfg_dir,'//','/'),self.out_dir,'')

      ;; Parse all config files, make fitsheaders for all output files
      for icfg = 0, Ncfg-1 do begin
        
        red_progressbar, icfg, Ncfg, progress_msg, clock=clock, /predict

        spawn, 'cat '+cfg_files[icfg], cfg
        Nlines = n_elements(cfg)

        ;; Get prpara (processing parameters) from the global config keywords
        Gstart = (where(strmatch(cfg,'}*'),Nmatch))[Nmatch-1] + 1
        undefine, Gprpara       ; Start fresh
        for iline = Gstart, Nlines-1 do begin
          cfgsplit = strsplit(cfg[iline], '=', /extract)
          case cfgsplit[0] of
            ;; Make list of possible global keywords complete!
            'ARCSECPERPIX'    : red_append, Gprpara,  cfgsplit[0] + '=' + cfgsplit[1]
            'BASIS'           : red_append, Gprpara,  cfgsplit[0] + '=' + cfgsplit[1]
            'DATA_TYPE'       : red_append, Gprpara,  cfgsplit[0] + '=' + cfgsplit[1]
            'DATE_OBS'        : red_append, Gprpara,  cfgsplit[0] + '=' + cfgsplit[1]
            'FAST_QR'         : red_append, Gprpara,  cfgsplit[0]
            'FILE_TYPE'       : begin
              file_type = cfgsplit[1]
              red_append, Gprpara, cfgsplit[0] + '=' + file_type
            end
            'FIT_PLANE'       : red_append, Gprpara, cfgsplit[0]
            'FPMETHOD'        : red_append, Gprpara, cfgsplit[0] + '=' + cfgsplit[1]
            'GETSTEP'         : red_append, Gprpara, cfgsplit[0] + '=' + cfgsplit[1]
            'GET_PSF'         : red_append, Gprpara, cfgsplit[0]
            'GET_ALPHA'       : red_append, Gprpara, cfgsplit[0]
            'GET_PSF_AVG'     : red_append, Gprpara, cfgsplit[0]
            'GRADIENT'        : red_append, Gprpara, cfgsplit[0] + '=' + cfgsplit[1]
            'IMAGE_NUMS'      : red_append, Gprpara, cfgsplit[0] + '=' + '['+cfgsplit[1]+']'
            'MAX_LOCAL_SHIFT' : red_append, Gprpara, cfgsplit[0] + '=' + cfgsplit[1]
            'MODES'           : red_append, Gprpara, cfgsplit[0] + '=[' + cfgsplit[1]+']'
            'NEW_CONSTRAINTS' : red_append, Gprpara, cfgsplit[0]
            'NUM_POINTS'      : red_append, Gprpara, cfgsplit[0] + '=' + cfgsplit[1]
            'PIXELSIZE'       : red_append, Gprpara, cfgsplit[0] + '=' + cfgsplit[1]
            'SIM_X'           : red_append, Gprpara, cfgsplit[0] + '=[' + cfgsplit[1]+']'
            'SIM_Y'           : red_append, Gprpara, cfgsplit[0] + '=[' + cfgsplit[1]+']'
            'TELESCOPE_D'     : red_append, Gprpara, cfgsplit[0] + '=' + cfgsplit[1]
            ;; Ignored config lines:
            'PROG_DATA_DIR' : 
;            '' : 
;            '' : 
;            '' : 
            else : begin
              print, inam + ' : Unknown cfg parameter:'
              print, cfg[iline]
              stop
            end
          endcase
        endfor                  ; iline
        
        
        ;; Parse objects
        Ostarts = where(strmatch(cfg,'object{*'),Nobj)
        for iobj = 0, Nobj-1 do begin

          ;; Make a generic header. Data type and dimensions can be
          ;; added when we read the output. 
          red_mkhdr, head, 0        

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
                'WAVELENGTH'  : red_append, prpara, cfgsplit[0] + '=' + cfgsplit[1]
                'WEIGHT'      : red_append, prpara, cfgsplit[0] + '=' + cfgsplit[1]
                'OUTPUT_FILE' : begin
                  output_file = cfgsplit[1]
                  red_append, prpara, cfgsplit[0] + '=' + output_file
                end
                else : begin
                  print, inam + ' : Unknown cfg parameter:'
                  print, cfg[iline]
                  stop
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
          fxaddpar, head, 'FILENAME', file_basename(output_file), ' MOMFBD restored data'


          ;; Channels
          if Nchan gt 0 then prmode = 'Phase Diversity' else prmode = ''
          for ichan = 0, Nchan-1 do begin
            
            Cstart = Cstarts[ichan]
            Cend   = Cstart + (where(strmatch(cfg[Cstart:*],'*}*'),Nmatch))[0]
            
            for iline = Cstart, Cend do begin
              cfgline = strtrim(cfg[iline], 2)
              if ~strmatch(cfgline,'*{') and ~strmatch(cfgline,'*}') then begin
                cfgsplit = strsplit(cfgline, '=', /extract)
                case cfgsplit[0] of
                  'ALIGN_CLIP'        : red_append, prpara, cfgsplit[0] + '=[' + cfgsplit[1]+']'
                  'DARK_NUM'          : red_append, prpara, cfgsplit[0] + '=' + cfgsplit[1]
                  'DARK_TEMPLATE'     : red_append, prpara, cfgsplit[0] + '=' + cfgsplit[1]
                  'FILENAME_TEMPLATE' : red_append, prpara, cfgsplit[0] + '=' + cfgsplit[1]
                  'GAIN_FILE'         : red_append, prpara, cfgsplit[0] + '=' + cfgsplit[1]
                  'IMAGE_DATA_DIR'    : red_append, prpara, cfgsplit[0] + '=' + cfgsplit[1]
                  'INCOMPLETE'        : red_append, prpara, cfgsplit[0]
                  'XOFFSET'           : red_append, prpara, cfgsplit[0] + '=' + cfgsplit[1]
                  'YOFFSET'           : red_append, prpara, cfgsplit[0] + '=' + cfgsplit[1]
                  'DISCARD'           : red_append, prpara, cfgsplit[0] + '=' + cfgsplit[1]
                  'DIVERSITY'         : red_append, prpara, cfgsplit[0] + '=' + cfgsplit[1]
                  else : begin
                    print, inam + ' : Unknown cfg parameter:'
                    print, cfg[iline]
                    stop
                  end
                endcase
              endif
            endfor              ; iline

          endfor                ; ichan

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
          fxaddpar, head, 'FILLED', 1, 'Missing pixels have been filled.'

          

          ;; Write the header file
          fxwrite, header_file, head 

          
        endfor                  ; iobj
        
      endfor                    ; icfg

    endfor                      ; ipref

  endfor                        ; idir

end

