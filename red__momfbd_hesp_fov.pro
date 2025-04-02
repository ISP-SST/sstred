; docformat = 'rst'

;+
; Setup for MOMFBD-processing CRISP/CHROMIS in the HeSP FOV as a
; single subfield.
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
;    cfgfiles : in, type=string
; 
;      Existing MOMFBD config files to modify.
; 
; 
; :Keywords:
; 
;   
;    coord : in, optional, type="lonarr(2)"
;   
;      The center pixel coordinates of the HeSP FOV. Use one of coord
;      or cubefile.
; 
;    cubefile : in, optional, type=string
;
;      The path to a fitscube to be used to select the position of the
;      HeSP FOV. Use one of coord or cubefile.
; 
;    fov_size : in, optional, type=float
; 
;      The size of the momfbd subfield in arcsec. Default 9 arcsec
;      unless num_points is given.
; 
;    num_points : in, optional, type=integer
; 
;      The size of the momfbd subfield in pixels. 
; 
;    cfg_tag : in, optional, type=string, default="[xc,yc,size]"
; 
;      A short string to be added to the "cfg" part of the path. A
;      default tag is constructed from the center coordinates of the
;      FOV and the FOV size.
; 
; :History:
; 
;    2025-01-28 : MGL. First version.
; 
;-
pro red::momfbd_hesp_fov, cfgfiles $
                          , cfg_tag = cfg_tag $
                          , cubefile = cubefile $
                          , coord = coord $
                          , fov_size = fov_size $ 
                          , iframe = iframe $
                          , iscan = iscan $
                          , istokes = istokes $
                          , ituning = ituning $
                          , nmodes = nmodes $
                          , num_points = num_points
  
  ;; Name of this subprogram
  inam = red_subprogram(/low, calling = inam1)  

  cfg = redux_readcfg(cfgfiles[0])

  if n_elements(fov_size) ne 0 && n_elements(num_points) ne 0 then begin
    print, inam + ' : Please give only one of fov_size or num_points.'
    return
  endif
  
  if n_elements(num_points) eq 0 then begin
    ;; Set NUM_POINTS based on fov_size and ARCSECPERPIX from the
    ;; file.
    if n_elements(fov_size) eq 0 then fov_size = 9.0 ; Deafult approximate HeSP FOV size
    arcsecperpix = float(redux_cfggetkeyword(cfg, 'ARCSECPERPIX'))
    num_points = round(fov_size/arcsecperpix)

    ;; Need num_points to be a multiple of 4
    num_points += red_odd(num_points)
    num_points += red_odd(num_points/2) * 2
  endif

  case 1 of

    n_elements(cubefile) eq 0 && n_elements(coord) eq 0 : begin
      print, inam + ' : Please give one of cubefile or coord.'
      return
    end

    n_elements(cubefile) ne 0 && n_elements(coord) ne 0 : begin
      print, inam + ' : Please give only one of cubefile or coord.'
      return
    end

    else :
    
  endcase

  if n_elements(cubefile) ne 0 then begin

    ;; Display an image from the cube and let the user click to select
    ;; the center coordinates. Note: need to figure out how the cube
    ;; was cropped from the raw FOV.
    
    ;; Possibly define more keywords to specify what kind of image
    ;; from the cube to display. Assume HeSP FOV is fairly well
    ;; centered in the CRISP/CHROMIS FOVs.

    ;; If user provides also a suitable HeSP image, it might be
    ;; possible to automatically find the best fit FOV.

    red_fitscube_getframe, cubefile, frame $
                           , iframe = iframe $
                           , iscan = iscan $
                           , istokes = istokes $
                           , ituning = ituning

    tighttv, frame, /scroll
    print, 'Click the center of the wanted FOV'

    dims_orig = size(frame, /dim)

    ;; We can infer cropping info from SIM_X/SIM_Y/SIM_XY and
    ;; NUM_POINTS. Note we can check alignment between cube frame and
    ;; raw data.
    hdr = headfits(cubefile)

    num_points_orig = redux_cfggetkeyword(cfg, 'NUM_POINTS')
    sim_xy = redux_cfggetkeyword(cfg, 'SIM_XY', count = cnt)
    if cnt eq 0 then begin
      sim_x = rdx_str2ints(redux_cfggetkeyword(cfg, 'SIM_X'))
      sim_y = rdx_str2ints(redux_cfggetkeyword(cfg, 'SIM_Y'))
      llx = min(sim_x) - num_points_orig/2
      lly = min(sim_y) - num_points_orig/2
      urx = max(sim_x) + num_points_orig/2
      ury = max(sim_y) + num_points_orig/2
    endif else begin
      sim_xy = rdx_str2ints(sim_xy)
      Nsubfields = n_elements(sim_xy)/2
      sim_xy = reform(sim_xy, 2, Nsubfields)
      llx = min(sim_xy[0,*]) - num_points_orig/2
      lly = min(sim_xy[0,*]) - num_points_orig/2
      urx = max(sim_xy[1,*]) + num_points_orig/2
      ury = max(sim_xy[1,*]) + num_points_orig/2
    endelse
   
    stop
    
    coord = [0, 0]
    stop
    
  endif

  
  ;; Directory for the new cfg files
  splt = strsplit(file_dirname(cfgfiles[0]), '/', /extract)
  if n_elements(cfg_tag) eq 0 then cfg_tag = strjoin(strtrim([coord, num_points], 2), ',')
  splt[-1] = splt[-1] + '_' + cfg_tag
  outdir = strjoin(splt, '/') + '/'
  file_mkdir, outdir+'/results'

  ;; For calling prepmomfbd_fitsheaders below:
  momfbddir = splt[0]
  dir = splt[1]
  pref = splt[2]
  
  for icfg = 0, n_elements(cfgfiles)-1 do begin

    ;; cfgfiles[0] was already read above.
    if icfg gt 0 then cfg = redux_readcfg(cfgfiles[icfg])

    ;; Overwrite the old num_points
    redux_cfgaddkeyword, cfg, 'NUM_POINTS', num_points

    ;; Modes
    if n_elements(Nmodes) ne 0 then begin
      redux_cfgaddkeyword, cfg, 'MODES', '2-'+strtrim(Nmodes+1, 2)
      redux_cfgaddkeyword, cfg, 'SORT_MODES'
    endif
    
    ;; Remove old SIM_X, SIM_Y, SIM_XY and set SIM_XY to coord.
    redux_cfgdelkeyword, cfg, 'SIM_X'
    redux_cfgdelkeyword, cfg, 'SIM_Y'
    redux_cfgdelkeyword, cfg, 'SIM_XY'
    redux_cfgaddkeyword, cfg, 'SIM_XY', strjoin(strtrim(coord, 2), ',')
    
    ;; Write the new cfg file
    outfile = outdir + file_basename(cfgfiles[icfg])
    redux_writecfg, outfile, cfg

  endfor                        ; icfg
  
  ;; Make .fitsheaders files
  self -> prepmomfbd_fitsheaders, dirs = dir, momfbddir = momfbddir, cfg_tag = cfg_tag
  
end

cd, '/scratch/mats/2024-06-05/CRISP/'
a = crisp2red("config.txt", /dev, /no)

cubefile = 'cubes_scan/nb_6173_2024-06-05T10:46:34_scan=0_stokes_corrected.fits'

cfgfiles = file_search('momfbd_nopd/10:46:34/6173/cfg/*.cfg', count = Ncfg)
if Ncfg eq 0 then stop
;a -> momfbd_hesp_fov, cfgfiles, cubefile = cubefile
a -> momfbd_hesp_fov, cfgfiles, coord = [512, 512]

end
