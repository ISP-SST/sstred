; docformat = 'rst'

;+
; Wrapper method for make_nb_cube with or without demodulation.
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
;    nopolarimetry : in, optional, type=boolean
;
;       For a polarimetric dataset, don't make a Stokes cube.
;       Instead combine all LC states for both cameras into a single
;       NB image per tuning, producing a cube similar to that for a
;       data set without polarimetry. (For a nonpolarlimetric dataset,
;       no effect.)
;
  
;   
;   
; 
; 
; :History:
; 
; 
; 
;-
pro red::make_nb_cube, wcfile $
                       , clips = clips  $
                       , cmap_fwhm = cmap_fwhm $
                       , fitpref_time = fitpref_time $
                       , integer = integer $
                       , intensitycorrmethod = intensitycorrmethod $ 
                       , nearest = nearest $
                       , noaligncont = noaligncont $  
                       , nocavitymap = nocavitymap $
                       , nocrosstalk = nocrosstalk $
                       , noflipping = noflipping $
                       , nomissing_nans = nomissing_nans $
                       , nopolarimetry = nopolarimetry $
                       , noremove_periodic = noremove_periodic $
                       , nostretch = nostretch $
                       , notimecor = notimecor $
                       , nthreads = nthreads $
                       , odir = odir $
                       , overwrite = overwrite $
                       , redemodulate = redemodulate $
                       , tiles = tiles $
                       , wbsave = wbsave
  
  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)

  ;; Deprecated keyword:
  if n_elements(notimecor) gt 0 then begin
    print, inam + ' : Keyword notimecor is deprecated. Use intensitycorrmethod="none" instead.'
    return
  endif

  ;; Default keywords
  if n_elements(cmap_fwhm) eq 0 then fwhm = 7.0 else fwhm = cmap_fwhm
  if n_elements(tiles) eq 0 or n_elements(clips) eq 0 then begin
    tiles = [8, 16, 32, 64, 84]
    clips = [8, 12,  4,  2,  1]
  endif
  if n_elements(nthreads) eq 0 then nthreads = 1 ; Default single thread

  ;; Make prpara
  red_make_prpara, prpara, clips         
  red_make_prpara, prpara, integer
  red_make_prpara, prpara, intensitycorrmethod
  red_make_prpara, prpara, cmap_fwhm
  red_make_prpara, prpara, nearest
  red_make_prpara, prpara, noaligncont 
  red_make_prpara, prpara, nocavitymap 
  red_make_prpara, prpara, nomissing_nans
  red_make_prpara, prpara, nostretch 
  red_make_prpara, prpara, np           
  red_make_prpara, prpara, overwrite
  red_make_prpara, prpara, tiles        
  red_make_prpara, prpara, wcfile

  
  ;; Camera/detector identification
  self -> getdetectors
  wbindx      = where(strmatch(*self.cameras,'*-W'))
  wbcamera    = (*self.cameras)[wbindx[0]]
  wbdetector  = (*self.detectors)[wbindx[0]]
  nbindx      = where(strmatch(*self.cameras,'*-[NTR]')) 
  nbcameras   = (*self.cameras)[nbindx]
  nbdetectors = (*self.detectors)[nbindx]
  Nnbcams     = n_elements(nbcameras)

;  nbtindx     = where(strmatch(*self.cameras,'*-T')) 
;  nbtcamera   = (*self.cameras)[nbtindx[0]]
;  nbtdetector = (*self.detectors)[nbtindx[0]]
;  nbrindx     = where(strmatch(*self.cameras,'*-R')) 
;  nbrcamera   = (*self.cameras)[nbrindx[0]]
;  nbrdetector = (*self.detectors)[nbrindx[0]]

  instrument = (strsplit(wbcamera, '-', /extract))[0]

  polarimetric_data = self -> polarimetric_data()

  ;; How to handle small scale variations in cavity maps.
  case instrument of
    'Crisp' : begin
      ;; We do correct for the small scale cavity map in CRISP data.
      ;; (We should get this from earlier meta data?)
      remove_smallscale = 1
    end
    'Chromis' : begin
      remove_smallscale = 0
    end
  endcase

    ;; Read the header from the corrected WB cube. Variables begin with
  ;; WC for Wideband Cube. 
  if ~file_test(wcfile) then begin
    print, 'WB cube missing, please run make_wb_cube.'
    print, wcfile
    retall
  endif
  wchead = red_readhead(wcfile)
  ;; Read parameters from the WB cube
  fxbopen, bunit, wcfile, 'MWCINFO', bbhdr
  fxbreadm, bunit, row = 1 $
            , ['ANG', 'CROP', 'FF', 'GRID', 'ND', 'SHIFT', 'TMEAN', 'X01Y01', 'DIRECTION'] $
            ,   ANG, wcCROP, wcFF, wcGRID, wcND, wcSHIFT, wcTMEAN, wcX01Y01,   direction
  ;; Note that the strarr wfiles cannot be read by fxbreadm! Put it in
  ;; wbgfiles (WideBand Global).
  fxbread, bunit, wbgfiles, 'WFILES', 1
  fxbclose, bunit
  if self.filetype eq 'MIXED' then wbgfiles = strtrim(wbgfiles, 2)

    ;; Default for wb cubes without direction parameter
  if n_elements(direction) eq 0 then direction = 0

  ;; CRISP or CHROMIS with W(D)TR cameras
  else : self -> make_nb_cube_stokes, wcfile, _extra = _extra

  endcase
  
end
