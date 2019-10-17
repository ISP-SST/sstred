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
;-
pro red::fitscube_cmapcorr, fname $
                            , flip = flip $
                            , interpolmethod = interpolmethod $
                            , outname = outname $
                            , overwrite = overwrite $
                            , _extra = extra 
  
  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

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
    outname = dir_in + '/' + file_basename(fname, '_im.fits') + '_cmapcorr_im.fits'
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
  if bitpix gt 0 then stop    ; We have to change the header if we want to handle integer input.
  
  ;; Change what needs to be changed in the header
  ;; Add info about this step
  self -> headerinfo_addstep, outhdr $
                              , prstep = 'Interpolate' $
                              , prpara = prpara $
                              , prproc = inam


  
  wcs = mrdfits(fname,'WCS-TAB',whdr)
  coords = wcs.hpln_hplt_wave_time
;  hpln   = reform(coords[0,*,*,0,0,*])
;  hplt   = reform(coords[1,*,*,0,0,*])
  lambda = reform(coords[2,0,0,*,0,0])
;  pol    = reform(coords[3,0,0,0,*,0])
;  time   = reform(coords[4,0,0,*,0,*])

  
  dlambda = mrdfits(fname, 'WCSDVARR', chdr, status = status, /silent)
  
;  lambda_corr[ix,iy,itun,ipol,iscan] = lambda[itun] + dlambda[ix,iy,iscan]

  ;; Initialize the output file
  self -> fitscube_initialize, outname, outhdr, lun, fileassoc, dims 


  
  cube_in  = fltarr(Nx, Ny, Nwav)
  cube_out = fltarr(Nx, Ny, Nwav)

  iprogress = 0
  Nprogress = Nscans*Nstokes
  for iscan = 0, Nscans-1 do begin
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

      for ix = 0, Nx-1 do for iy = 0, Ny-1 do begin
        ;; Interpolate
        case interpolmethod of
          'interpol' : begin
            cube_out[ix, iy, *] = interpol(cube_in[ix, iy, *] $
                                           , lambda + dlambda[ix,iy,iscan] $
                                           , lambda $
                                           , _strict_extra = extra)
          end
          'spline' : begin
            cube_out[ix, iy, *] = spline(lambda + dlambda[ix,iy,iscan] $
                                         , cube_in[ix, iy, *] $
                                         , lambda $
                                         , _strict_extra = extra)
          end
          'quadterp' : begin
            quadterp, lambda + dlambda[ix,iy,iscan] $
                      , reform(cube_in[ix, iy, *]) $
                      , lambda $
                      , yint $
                      ,_strict_extra = extra
            cube_out[ix, iy, *] = yint
          end
          else: begin
            print, inam + ' : Method not implemented: '+interpolmethod
            stop
          end
        endcase
      endfor                    ; ix,iy
    
      for iwav = 0, Nwav-1 do begin
        ;; Write
        self -> fitscube_addframe, fileassoc $
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
