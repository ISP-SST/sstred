; docformat = 'rst'

;+
; Make an integer version of a floating point fitscube.
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
;     delete : in, optional, type=boolean
;
;       Delete the input file if this keyword is set.
; 
;     flip : in, optional, type=boolean
;   
;       Produce a flipped version if this keyword is set. 
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
; 
; 
; :History:
; 
;   2018-11-02 : MGL. First version.
; 
;-
pro red::fitscube_integer, fname $
                           , delete = delete $
                           , flip = flip $
                           , outname = outname $
                           , overwrite = overwrite

  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)
 
  dir_in = file_dirname(fname)
  
  if ~file_test(fname) then begin
    print, inam + ' : File does not exist: '
    print, fname
    return
  endif

  hdr = headfits(fname)
  outhdr = hdr
  
  if fxpar(hdr, 'BITPIX') gt 0 then begin
    print, inam + ' : The input file is already an integer cube: '
    print, fname
    return
  endif
  
  if n_elements(outname) eq 0 then begin
    outname = dir_in + '/' + file_basename(fname, '_im.fits') + '_int_im.fits'
  endif

  if file_test(outname) and ~keyword_set(overwrite) then begin
    print, inam + ' : The output file already exists: '
    print, outname
    return
  endif 

  dir_out = file_dirname(outname)
  file_mkdir, dir_out

  
  dims = fxpar(hdr, 'NAXIS*', Naxis)
  Nx      = dims[0]
  Ny      = dims[1]
  Nwav    = dims[2]
  Nstokes = dims[3]
  Nscans  = dims[4]
  
  ;; Calculate BZERO and BSCALE
  red_fitspar_getminmax, fname $
                         , datamax = datamax $
                         , datamin = datamin
  
  arraymin = -32000.
  arraymax =  32000.

  bscale = (datamax-datamin)/(arraymax-arraymin)
  bzero = datamin - arraymin*bscale
  
  ;; Set keywords for rescaled integers
  red_fitsaddkeyword, outhdr, 'BITPIX', 16
  red_fitsaddkeyword, outhdr, 'BZERO',  bzero
  red_fitsaddkeyword, outhdr, 'BSCALE', bscale

  ;; Added later if needed
  red_fitsdelkeyword, outhdr, 'VAR_KEYS'

  ;; Initialize the output file
  self -> fitscube_initialize, outname, outhdr, lun, fileassoc, dims 

  ;; Copy the data
  Nframes = long(Nscans) * long(Nstokes) * long(Nwav)
  for iframe = 0, Nframes - 1 do begin

    red_progressbar, iframe, Nframes $
                     , /predict $
                     , 'Copy the data'

    ;; Read
    red_fitscube_getframe, fname, im, iframe = iframe

    ;; Write
    self -> fitscube_addframe, fileassoc $
                               , fix(round((temporary(im)-BZERO)/BSCALE)) $
                               , iframe = iframe
    
  endfor                        ; iframe

  ;; Close
  self -> fitscube_finish, lun

  ;; Copy extensions
  red_fitscube_copyextensions, fname, outname

  ;; Add info about this step
  ;; Put BZERO, BSCALE, and BITPIX in prpara or prref?
  self -> headerinfo_addstep, hdr $
                              , prstep = 'BZERO-BSCALE-TRUNCATION' $
;                              , prpara = prpara $
                              , prproc = inam


  if keyword_set(flip) then begin
    ;; Make a flipped version
    red_fitscube_flip, outname $
                       , flipfile = flipfile $
                       , overwrite = overwrite
  endif

  if keyword_set(delete) then begin
    print, inam + ' : Deleting float cube:'
    print, fname
    file_delete, fname
    ;; Possibly delete a flipped version of fname?
  endif

  print, inam + ' : Integer cube stored in:'
  print, outname
  if keyword_set(flip) then print, flipfile

end

;; WHere do the errors happen?

;; [====================] 100.0% in 08m33s (00s remaining): Write the flipped cube, istokes,iscan,ix=0,3,973                                                                        
;; % MRD_HREAD: Warning-Invalid characters in header
;; % Type conversion error: Unable to convert given STRING to Integer.
;; % Detected at: RED_FITSGETKEYWORD  130 /home/mats/idl/bin/crispred/red_fitsgetkeyword.pro
;; % Type conversion error: Unable to convert given STRING to Integer.
;; % Detected at: RED_FITSGETKEYWORD  130 /home/mats/idl/bin/crispred/red_fitsgetkeyword.pro
;; /scratch/mats/2016.09.19/CRISP-aftersummer//cubes_nb/nb_6302_2016-09-19T09:30:20_scans=25-28_stokes_corrected_int_sp.fits
;; 
