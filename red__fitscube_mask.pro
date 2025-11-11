; docformat = 'rst'

;+
; Define and apply a mask on the FOV of a fitscube and then crop to fit.
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
;    infile : in, type=string
; 
;      The path to the input fitscube.
; 
; 
; :Keywords:
;
;    ituning : in, optional, type=integer
;
;      The tuning index of the Stokes image displayed when defining
;      the mask. If not given, the core tuning will be used.
;
;    noflip : in, optional, type=boolean
;   
;      If set, don't make any spectral, "flipped" data cubes.
; 
;    overwrite : in, type=boolean
;   
;      Overwrite an existing file.
; 
;    tag : in, type=string, default='masked'
; 
;      If not overwriting the original file, add this tag to the
;      original filename.
; 
; 
; :History:
; 
;    2025-11-11 . MGL. First version.
; 
;-
pro red::fitscube_mask, infile $
                        , ituning = ituning $
                        , noflip = noflip $
                        , overwrite = overwrite $
                        , tag = tag

  if n_elements(tag) eq 0 then tag = 'masked'
  outfile = red_strreplace(infile, '_im.fits', '_' + tag + '_im.fits')

  if ~keyword_set(overwrite) && file_test(outfile) then begin
    red_message, 'Output file already exists. Use other tag or use /overwrite.'
    return
  endif
  
  indir = file_dirname(infile)
  tmpfile = indir + '/tmp'+cgTimeStamp(Random_Digits=5, /UTC)+'.fits'

  print, infile
  print, outfile
  print, tmpfile

  
  
  hdr = headfits(infile)
  naxis = fxpar(hdr, 'NAXIS*')

  Nx = naxis[0]
  Ny = naxis[1]
  Ntuning = naxis[2]
  Nstokes = naxis[3]
  Nscans = naxis[4]

  polarimetric = Nstokes gt 1

  red_fitscube_getwcs, infile, coordinates = coordinates
  lambda = coordinates[*, 0].wave[0,0] ;Wavelengths in nm

  if n_elements(ituning) eq 0 then begin
    red_message, 'Find core tuning point'
    red_fitscube_statistics, infile, frame_statistics, cube_statistics
    mn = min(frame_statistics[*,0,0].datamedn, ituning)
  endif
  
  ;; Base mask on total Stokes V for the selected tuning point 
  dispim = 0.0
  iStokes = 3
  for iscan = 0, Nscans-1 do begin
    red_fitscube_getframe, infile, thisframe $
                           , ituning = ituning $
                           , istokes = istokes $
                           , iscan = iscan
    dispim += thisframe
  endfor                        ; iscan
  origmask = finite(dispim)
  dispim = bytscl(red_histo_opt(dispim))
  
  red_message, 'Define the mask'
  XROI, dispim, Regions_Out=ROIs, /Block, title = 'Select area.'
  mask = 0*origmask
  for iroi = 0, n_elements(ROIs)-1 do $
     mask or= ROIs[iroi] -> ComputeMask(Dimensions=[Nx, Ny], Mask_Rule=2) 
  mask and= origmask

  bb = red_boundingbox(mask)
  Obj_Destroy, ROIs

  ;; Copy the original cube to a temporary file name and then
  ;; overwrite the frames with the masked frames. This way we maintain
  ;; all the metadata.
  file_copy, infile, tmpfile
  Nframes = long(Nscans)*long(NStokes)*long(Ntuning)
  for iframe = 0L, Nframes-1 do begin
    red_progressbar, iframe, Nframes, 'Copying frames'
    red_fitscube_getframe, infile, thisframe, iframe = iframe 
    red_fitscube_addframe, tmpfile, mask*thisframe, iframe = iframe 
  endfor                        ; iframe

  ;; Add info about this step
  hdr = headfits(tmpfile)
  inam = red_subprogram(/low, calling = inam1)
  red_make_prpara, prpara, ROIs[-1]
  self -> headerinfo_addstep, hdr $
     , prstep = 'MASKING' $
     , prpara = prpara $
     , prproc = inam
  red_fitscube_newheader, tmpfile, hdr

  if ~keyword_set(nomissing_nans) then begin
    ;; Set padding pixels to missing-data, i.e., NaN.
    self -> fitscube_missing, tmpfile $
       , /noflip $
       , missing_type = 'nan' 
  endif

  ;; Then run fitscube_crop on the temporary file with a sensible
  ;; output file name. Either with a standard tag or possibly
  ;; overwriting the original file. Use the bb as input to
  ;; fitscube_crop.
  self -> fitscube_crop, tmpfile $
     , nospectral = noflipping $
     , outfile = outfile $
     , /overwrite $
     , roi = bb[[0, 2, 1, 3]]

  ;; Do we need to generate also a matching WB file? Could use the
  ;; same bb but not care about the masking --> more context. Use same
  ;; tag in the new wb file name, or overwrite if that's what you want.

  file_delete, tmpfile
  
  print, outfile
  
end

dir = '/scratch_beegfs/mats/NEW/2025-08-18/CRISP2/'
cd, dir
a = crisp2red(/dev, /no)

infile = 'cubes_nb/nb_6173_2025-08-18T08:52:59_08:52:59=1-4_stokes_corrected_im.fits'

a -> fitscube_mask, infile, ituning = 6, /overwrite

end
