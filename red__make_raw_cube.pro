; docformat = 'rst'

;+
; Make a de-rotated and de-stretched time-series FITS data cube with
; raw data.
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
; 
; 
; 
; 
; 
; :Keywords:
; 
;
;    integer : in, optional, type=boolean
;
;       Store as integers instead of floats. Uses the BZERO and BSCALE
;       keywords to preserve the intensity scaling.
;
;    overwrite : in, optional, type=boolean
;
;       Don't care if cube is already on disk, overwrite it
;       with a new version.
;   
;   
;   
; 
; 
; :History:
; 
;   2019-07-04 : MGL. First version.
; 
;-
pro red::make_raw_cube, oldname $
                        , integer = integer $
                        , mtf_deconvolve = mtf_deconvolve $
                        , neuralnet = neuralnet $
                        , nexp = nexp $
                        , overwrite = overwrite

  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  if n_elements(tiles) eq 0 or n_elements(clips) eq 0 then begin
    tiles = [8, 16, 32, 64, 128]
    clips = [8, 4,  2,  1,  1  ]*2
  endif

  if ~file_test(oldname) then begin
    print, inam + ' : File does not exist: ',  oldname
  endif

  oldhdr = headfits(oldname)
  oldbase = file_basename(oldname)
  olddir  = (strsplit(file_dirname(oldname), '/', /extract))[-1]

  dimensions = long(fxpar(oldhdr, 'NAXIS*'))
  Nx      = dimensions[0]
  Ny      = dimensions[1]
  Ntuning = dimensions[2]
  Nstokes = dimensions[3]
  Nscans  = dimensions[4]
  
  
  ;; Output
  odir = self.out_dir + '/cubes_raw/'
  file_mkdir, odir
  oname = red_strreplace(oldbase, 'corrected', 'raw')

  filename = odir+oname
  
  ;; Already done?
  if file_test(filename) then begin
    if keyword_set(overwrite) then begin
      print, inam + ' : Overwriting existing data cube:'
      print, filename
    endif else begin
      print, inam + ' : This data cube exists already:'
      print, filename
      return
    endelse
  endif
  
  hdr = oldhdr
  fxaddpar, hdr, 'BITPIX', -32
  newcube = fltarr(Nx, Ny, Ntuning, Nstokes, Nscans)
  oldcube = fltarr(Nx, Ny, Ntuning, Nstokes, Nscans)
  

  
  case olddir of
    'cubes_wb' : begin

      ;; Modify header

      ;; Get rotated cube from quicklook. Trust that it gets
      ;; orientation and rotation right.
      cam = fxpar(oldhdr, 'CAMERA')
      detector = fxpar(oldhdr, 'DETECTOR')
      detectorinfo = red_camerainfo(detector)   
      dataset = (strsplit(fxpar(oldhdr, 'STARTOBS'),'T',/extract))[1]
      self -> extractstates, oldname, states
      state = states.prefilter+'_'+states.tuning
      sc = red_fitsgetkeyword(oldname, 'SCANNUM', variable_values = variable_values)
      scannos = reform(variable_values.values)

      self -> quicklook, cam = cam $
                         , cube = cube $
                         , datasets = dataset $
                         , /derotate $
                         , mtf_deconvolve = mtf_deconvolve $
                         , neuralnet = neuralnet $
                         , nexp = nexp $
                         , /no_plot_r0 $
                         , scannos = scannos $
                         , use_states = state $
                         , x_flip = x_flip $
                         , y_flip = y_flip

      cube = reform(cube, detectorinfo.xsize, detectorinfo.ysize, 1, 1, Nscans, /overwrite)
      for iscan = 0, Nscans-1 do begin
        red_progressbar, iscan, Nscans, 'Align and destretch with old cube images'
        red_fitscube_getframe, oldname, oldframe, iscan = iscan
        oldcube[0, 0, 0, 0, iscan] = oldframe
        newframe = red_centerpic(cube[*, *, 0, 0, iscan], xsize = Nx, ysize = Ny)
        dxdy = red_alignoffset(oldframe/median(oldframe), newframe/median(newframe))
        newframe = red_centerpic(cube[*, *, 0, 0, iscan], xsize = Nx, ysize = Ny $
                                 , xoffset = round(dxdy[0]), yoffset = round(dxdy[1]))
        
        grid1 = red_dsgridnest(oldframe/median(oldframe), newframe/median(newframe), tiles, clips)
        newcube[0, 0, 0, 0, iscan] = red_stretch(temporary(newframe), grid1)

      endfor                    ; iframe

      ;; WCS stuff
      red_fitscube_getwcs, oldname $
                           , coordinates = wcs $
                           , distortions = distortions
    
    end
    else : stop
  endcase

  self -> fitscube_initialize, filename, hdr, lun, fileassoc, dimensions

  ;; Write the data
  for ituning = 0, Ntuning-1 do begin
    for istokes = 0, Nstokes-1 do begin
      for iscan = 0, Nscans-1 do begin
        self -> fitscube_addframe, fileassoc, newcube[*, *, 0, 0, iscan] $
                                   , ituning = ituning $
                                   , istokes = istokes $
                                   , iscan = iscan
      endfor                    ; iscan
    endfor                      ; istokes
  endfor                        ; ituning
  free_lun, lun
  print, inam + ' : Wrote file '+filename

  ;; Add WCS
  red_fitscube_addwcs, filename, wcs, dimensions = dimensions

  ;; Copy variable-keywords from the oldcube file. Approximate
  ;; DATE-XXX but shouldn't matter for this?
  self -> fitscube_addvarkeyword, filename, 'SCANNUM',  old_filename = oldname $
                                  , anchor = anchor 
  self -> fitscube_addvarkeyword, filename, 'ATMOS_R0', old_filename = oldname $
                                  , anchor = anchor 
  self -> fitscube_addvarkeyword, filename, 'AO_LOCK', old_filename = oldname $
                                  , anchor = anchor 
  self -> fitscube_addvarkeyword, filename, 'ELEV_ANG', old_filename = oldname $
                                  , anchor = anchor 
  self -> fitscube_addvarkeyword, filename, 'DATE-BEG', old_filename = oldname $
                                  , anchor = anchor
  self -> fitscube_addvarkeyword, filename, 'DATE-END', old_filename = oldname $
                                  , anchor = anchor
  self -> fitscube_addvarkeyword, filename, 'DATE-AVG', old_filename = oldname $
                                  , anchor = anchor

  if keyword_set(integer) then begin
    self -> fitscube_integer, filename $
                              , /delete $
                              , outname = outname $
                              , overwrite = overwrite
    filename = outname
  endif


  ;; Add statistics - or remove existing entries from header? Get
  ;; keywords from quicklook?
  red_fitscube_statistics, filename, /write ;$
;                           , angles = ang $
;                           , full = wcFF $
;                           , grid = wcGRID $
;                           , origNx = Nxx $
;                           , origNy = Nyy $
;;                             , percentiles = percentiles $
;                           , shifts = wcSHIFT 

  print, inam + ' : Raw data cube stored in:'
  print, filename

  
end

a = chromisred(/dev)

oldname = 'cubes_wb/wb_3950_2016-09-19T09:28:36_scans=67-73_corrected_im.fits'

a -> make_raw_cube, oldname, /over, /mtf

end
