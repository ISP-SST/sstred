; docformat = 'rst'

;+
; Make a fitscube that is corrected for wavelength distortions (cavity
; errors).
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
;     fname : in, type=string
; 
;       The name of the input file.
; 
; 
; :Keywords:
; 
;     flip : in, optional, type=boolean
;   
;       Produce a flipped version if this keywoed is set.
; 
;     interpolmethod  : in, optional, type=string, default='interpol'
;   
;       The interpolation method to use.
; 
;     outname : in, out, optional, type=string, default="Constructed from input file name"
;   
;       The name of the output file.
; 
;     overwrite : in, optional, type=boolean
;
;       Don't care if the output cube is already on disk, overwrite it
;       with a new version.
; 
;     _extra : in, optional 
;   
;      Any extra keywords are used when calling the chosen
;      interpolation method.
;
; 
; 
; :History:
; 
;   2018-11-02 : MGL. First version.
; 
;   2018-12-06 : MGL. Additional interpolation methods.
; 
;   2019-09-27 : MGL. Implement correction of individual cavity maps
;                for multiple prefilters.
; 
;-
pro red::fitscube_cmapcorr, fname $
                            , flip = flip $
                            , interpolmethod = interpolmethod $
                            , outname = outname $
                            , overwrite = overwrite $
                            , _extra = extra 
  
  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)                        

  dir_in = file_dirname(fname)

  if n_elements(interpolmethod) eq 0 then begin
    interpolmethod = 'interpol'
    if n_elements(extra) eq 0 then extra = {quadratic:1}
  endif
  
  ;; Make prpara
  red_make_prpara, prpara, interpolmethod
  if n_elements(extra) gt 0 then begin
    extra_tags = tag_names(extra)
    for itag = 0, n_elements(extra_tags)-1 do begin
      red_make_prpara, prpara, extra.(itag), paraname = extra_tags[itag]
    endfor                      ; itag
  endif
  
  if ~file_test(fname) then begin
    print, inam + ' : File does not exist: '
    print, fname
    return
  endif

  if n_elements(outname) eq 0 then begin
;    outname = dir_in + '/' + file_basename(fname, '_im.fits') + '_cmapcorr_im.fits'
    outname = dir_in + '/' + red_strreplace(file_basename(fname), '_corrected', '_corrected_cmapcorr')
  endif
  
  if file_test(outname) and ~keyword_set(overwrite) then begin
    print, inam + ' : The output file already exists: '
    print, outname
    return
  endif 

  dir_out = file_dirname(outname)
  file_mkdir, dir_out

  
  hdr = headfits(fname)
  outhdr = hdr

  dims = fxpar(hdr, 'NAXIS*', Naxis)
  Nx      = dims[0]
  Ny      = dims[1]
  Nwav    = dims[2]
  Nstokes = dims[3]
  Nscans  = dims[4]

  bitpix = fxpar(hdr, 'BITPIX')
  if bitpix gt 0 then stop      ; We have to change the header if we want to handle integer input.
  
  ;; Add info about this step
  self -> headerinfo_addstep, outhdr $
                              , prstep = 'SPECTRAL-DISTORTION-CORRECTION' $
                              , prpara = prpara $
                              , prproc = inam
  
  ;; Read WCS coordinates and distortions (the latter of which are the
  ;; cavity map wavelength distortions).
  red_fitscube_getwcs, fname $
                       , coordinates = coordinates $
                       , distortions = distortions $
                       , iscan = 0

  Ndist = n_elements(distortions)
  lambda = reform(coordinates[*,0].wave[0,0,*])

  ;; Are there any tunings for which we don't have distortions?
  for idist = 0, Ndist-1 do red_append, dist_indx, red_expandrange(distortions[idist].tun_index)
  nodist_index = cgSetDifference(indgen(Nwav),dist_indx, count = Nnodist)


  
  ;; Initialize the output file
  self -> fitscube_initialize, outname, outhdr, lun, fileassoc, dims 


  ;; Do one scan at a time
  cube_in  = fltarr(Nx, Ny, Nwav)
  cube_out = fltarr(Nx, Ny, Nwav)

  iprogress = 0
  Nprogress = Nscans*Nstokes
  for iscan = 0, Nscans-1 do begin

    ;; Read distortions for this scan
    red_fitscube_getwcs, fname $
                         , distortions = distortions $
                         , iscan = iscan

    for istokes = 0, Nstokes-1 do begin

      red_progressbar, iprogress, Nprogress, /predict  $
                       , 'Interpolating (iscan,istokes)=(' $
                       + strtrim(iscan, 2) + ',' + strtrim(istokes, 2) + ')'

      ;; Read an [Nx,Ny,Nwav] cube
      for iwav = 0, Nwav-1 do begin
        red_fitscube_getframe, fname, im $
                               , ituning = iwav $
                               , iscan = iscan $
                               , istokes = istokes
        cube_in[0, 0, iwav] = im
      endfor                    ; iwav

      for idist = 0, Ndist-1 do begin
        ;; Loop over multiple distortions, interpolate subcubes.

  
        dlambda = distortions[idist].wave
        tun_index = red_expandrange(distortions[idist].tun_index)
        
        for ix = 0, Nx-1 do for iy = 0, Ny-1 do begin
          ;; Interpolate
          case interpolmethod of
            'interpol' : begin
              cube_out[ix, iy, tun_index] = interpol(cube_in[ix, iy, tun_index] $
                                                     , lambda[tun_index] + dlambda[ix, iy] $
;;                                                     , lambda[tun_index] + dlambda[ix, iy, 0, 0, iscan] $
                                                     , lambda[tun_index] $
                                                     , _strict_extra = extra)
            end
            'spline' : begin
;;              cube_out[ix, iy, *] = spline(lambda + dlambda[ix,iy,iscan] $
              cube_out[ix, iy, *] = spline(lambda + dlambda[ix,iy] $
                                           , cube_in[ix, iy, *] $
                                           , lambda $
                                           , _strict_extra = extra)
            end
            'quadterp' : begin
;;              quadterp, lambda + dlambda[ix,iy,iscan] $
              quadterp, lambda + dlambda[ix,iy] $
                        , reform(cube_in[ix, iy, *]) $
                        , lambda $
                        , yint $
                        , _strict_extra = extra
              cube_out[ix, iy, *] = yint
            end
            else: begin
              print, inam + ' : Method not implemented: '+interpolmethod
              stop
            end
          endcase
        endfor                  ; ix,iy
        
      endfor                    ; idist

      ;; Just copy frames for tunings for which there are no
      ;; distortions (if any).
      for inodist = 0, Nnodist-1 do begin
        cube_out[0, 0, nodist_index[inodist]] = cube_in[*, *, nodist_index[inodist]] 
      endfor                    ; iwav
      
      for iwav = 0, Nwav-1 do begin
        ;; Write
        red_fitscube_addframe, fileassoc $
                               , cube_out[*, *, iwav] $
                               , ituning = iwav $
                               , iscan = iscan $
                               , istokes = istokes
      endfor                    ; iwav

      iprogress++
      
    endfor                      ; istokes
  endfor                        ; iscan

  self -> fitscube_finish, lun
  
  ;; Copy the variable keywords
  var_keys = red_fits_var_keys(hdr, count = Nkeys)

  for ikey = 0, Nkeys-1 do begin

    self -> fitscube_addvarkeyword, outname $
                                    , var_keys[ikey] $
                                    ,  old_filename = fname
  endfor                        ; ikey 

  ;; Copy WCS extension
  red_fits_copybinext, fname, outname, 'WCS-TAB'

  ;; Copy WB image extension (of a scan cube)
  if Nscans eq 1 then red_fitscube_copyextensions, fname, outname, ext_list = 'WBIMAGE'
  
  if keyword_set(flip) then begin
    ;; Make a flipped version
    red_fitscube_flip, outname $
                       , flipfile = flipfile $
                       , overwrite = overwrite
  endif

  print, inam + ' : Cavity-corrected cube stored in:'
  print, outname
  if keyword_set(flip) then print, flipfile
  
end

; Should make new statistics? In principle it's needed but then we
; have no way of avoiding the padding areas outside the FOV. Note that
; the output is not meant for archiving, so probably not needed.

a = chromisred(/dev)

fname = 'cubes_nb/nb_3950_2016-09-19T10:42:01_scans=0-4_corrected_im.fits'

a -> fitscube_cmapcorr, fname, /over
;, flip = flip $
;         , interpolmethod = interpolmethod $
;         , outname = outname $
;         , overwrite = overwrite $
;         , _extra = extra 

end
