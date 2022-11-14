; docformat = 'rst'

;+
; Turn regular momfbd cfg files made with prepmomfbd into multiple
; files for automatic mosaic observations.
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
; 
; 
; 
; 
; 
; :Keywords:
; 
;   
;   
;   
; 
; 
; :History:
; 
;   2022-11-06 : MGL. First version.
; 
;-
pro red::prepmomfbd_mosaic, dirs=dirs $
                            , momfbddir=momfbddir $
                            , pref = pref $
                            , no_delete = no_delete  $
                            , no_pd = no_pd 
  
  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)  
  
  if n_elements(dirs) gt 0 then begin
    dirs = [dirs] 
  endif else begin
    print, inam+' : ERROR : undefined data_dir'
    return
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
    
    ;; Base output location
    cfg_base_dir = self.out_dir + PATH_SEP() + momfbddir + PATH_SEP() + folder_tag + PATH_SEP()


    if n_elements(pref) gt 0 then begin
      prefs = pref
      Nprefs = n_elements(prefs)
    endif else begin
      prefs = file_basename(file_search(cfg_base_dir+'*', count = Nprefs))
    endelse
    
    for ipref = 0, Nprefs-1 do begin
      
      cfg_dir = cfg_base_dir + PATH_SEP() + prefs[ipref] + PATH_SEP() + 'cfg/'
      cfg_files = file_search(cfg_dir + 'momfbd_reduc_' + prefs[ipref] $
                              + '_?????.cfg', count = Ncfg)
      for icfg = 0, Ncfg-1 do begin

        red_progressbar, icfg, Ncfg, cfg_files[icfg], /predict

        scanno = (strsplit(file_basename(cfg_files[icfg],'.cfg'),'_',/extract))[-1]
        
        cfg = redux_readcfg(cfg_files[icfg])

        ;; Find out how many mosaic positions there are
        mdir = redux_cfggetkeyword(cfg, 'OBJECT1.CHANNEL0.IMAGE_DATA_DIR')
        files = file_basename(file_search(mdir+'/*_'+scanno+'_*fits', count = Nfiles))
        pos = strpos(files[0], '_mos')
        mosnums = strmid(files, pos+4, 2)
        mosnums = mosnums[uniq(mosnums, sort(mosnums))]
        Nmos = n_elements(mosnums)

        Nobjects = n_elements(redux_cfggetkeyword(cfg, 'OBJECT*'))

        ;; Make new results directories
        file_mkdir, cfg_dir + 'results_mos'+mosnums
        
        for imos = 0, Nmos-1 do begin

          ;; Make a cfg file for mosnums[imos]

          cfg_file_mos = red_strreplace(cfg_files[icfg], 'momfbd_reduc', 'momfbd_reduc_mos'+mosnums[imos])
          cfg_mos = redux_readcfg(cfg_files[icfg]) ; Read again rather than copy (by reference)

          ;; Change IMAGE_NUMS
          files = file_search(mdir+'/*_'+scanno+'_*_mos'+mosnums[imos]+'_*fits', count = Nfiles)
          self -> extractstates, files, states
          undefine, image_nums
          for istate = 0, Nfiles-1 do $
             red_append, image_nums, states[istate].framenumber + indgen(states[istate].nframes)
          image_nums = rdx_ints2str(image_nums)
          redux_cfgaddkeyword, cfg_mos, 'IMAGE_NUMS', image_nums
          
          for iobject = 0, nobjects-1 do begin

            ;; Change OUTPUT_FILE 
            ofile = redux_cfggetkeyword(cfg_mos, 'OBJECT'+strtrim(iobject, 2)+'.OUTPUT_FILE')
            ofile = red_strreplace(ofile, 'results', 'results_mos'+mosnums[imos])
            redux_cfgaddkeyword, cfg_mos, 'OBJECT'+strtrim(iobject, 2)+'.OUTPUT_FILE', ofile

            ;; Change FILENAME_TEMPLATE (all channels!) (OBJECT0 not
            ;; needed, mosNN not included in [WD]_nostate file names.
            if iobject gt 0 && imos gt 0 then begin
              ftemp = redux_cfggetkeyword(cfg_mos, 'OBJECT'+strtrim(iobject, 2)+'.CHANNEL0.FILENAME_TEMPLATE')
              ftemp = red_strreplace(ftemp, 'mos00', 'mos'+mosnums[imos])
              redux_cfgaddkeyword, cfg_mos, 'OBJECT'+strtrim(iobject, 2)+'.CHANNEL0.FILENAME_TEMPLATE', ftemp
            endif

          endfor                ; iobject

          redux_writecfg, cfg_file_mos, cfg_mos
          
        endfor                  ; imos

        if ~keyword_set(no_delete) then file_delete, cfg_files[icfg]
        
        self -> prepmomfbd_fitsheaders, /mosaic $
                                        , dirs = dirs[idir] $
                                        , momfbddir = momfbddir $
                                        , pref = pref[ipref] $
                                        , no_pd = no_pd 

        ;; Possibly add some mosaic-specific info to the fitsheader
        ;; files in prepmomfbd_fitsheaders.
        
      endfor                    ; icfg ;

    endfor                      ; ipref
    
  endfor                        ; idir
  
end



;; TODO: Do something about the .fitsheader files. Either run
;; prepmomfbd_fitsheaders after (or from!) this method or modify
;; fitsheader files corresponding to original cfg files and write to
;; the results_mosNN directories. In either case, the fitsheader files
;; should have some added information about the mosNN number. And
;; ideally about the grid positions because we don't know if a
;; six-tile mosaic is 2x3 or 3x2.

