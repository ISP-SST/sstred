; docformat = 'rst'

;+
; Get info from flats/spectral_flats/.
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
; :Keywords:
; 
;    cmap : out, optional, type=array
;
;      The cavity error map, units [nm] with the sign set so
;      lambda_correct = lambda + cmap.
;
;    detector : in, optional, type=string
;
;      Limit results to this detector.
;
;    fit_results_struct : out, optional, type=struct
;
;      The output from fitgains as a struct with members PARS (the fit
;      parameters as a 3D cube), XL, YL (), and ONAME (the output
;      cavity free flats file names). Only returned if the length of
;      fit_results_files is 1.
;
;    fit_results_files : out, optional, type=strarr
;
;      The path to the fitgains results save file(s)
;
;    flats_cube : out, optional, type=fltarr
;
;      A 4D data cube [Ntuning,Nstokes,Nx,Ny].
;
;    flats_files : out, optional, type=array
;
;      The path to the prepflatscube results save file(s)
;
;    flats_names : out, optional, type=array
;
;      The paths to the summed flats files that go into the flats
;      cube.
;
;    flats_wav : out, optional, type=array
;
;      The tunings from the nominal line center.
;
;    prefilter : in, optional, type=string
;
;      Limit results to this prefilter.
;
;   
;   
; 
; 
; :History:
; 
;   2025-12-01 . MGL. First version.
; 
;-
pro red::get_spectral_flats_info, cmap = cmap $
                                  , detector = detector $
                                  , fit_results_struct = fit_results_struct $
                                  , fit_results_files = fit_results_files $
                                  , flats_cube = flats_cube $
                                  , flats_files = flats_files $
                                  , flats_names = flats_names $
                                  , flats_wav = flats_wav $
                                  , prefilter = prefilter

  if n_elements(detector) eq 0 then detector = '*'
  if n_elements(prefilter) eq 0 then prefilter = '*'

  search_prefix = self.out_dir + '/flats/spectral_flats/' $
                  + detector $
                  + ['_', '_*'] $ ; Could be exposure time and detector gain or just an underscore
                  + prefilter 
  

  if arg_present(flats_cube) || arg_present(flats_file) then begin

    ;; The flats.sav file contains the variables CUB, NAMELIST, and
    ;; WAV. These variables are also available in individual fits/text
    ;; files. CUB is much larger than the other two so unless we
    ;; *only* want the other two, we might as well read (or return the
    ;; name of) the sav file.
    
    flats_files = red_uniquify(file_search(search_prefix + '_flats.sav') $
                               , count = Nfiles)
    
    if Nfiles eq 0 then begin
      red_message, 'No flats cube files found. Run prepflatcubes first!'
      ;; Might as well give up here because without the flats cube
      ;; there are no fits.
      retall                    
    endif

    if arg_present(flats_cube) || arg_present(flats_names) ||  arg_present(flats_wav) then begin

      if Nfiles gt 1 then begin
        red_message, 'More than one flats.sav file found.'
      endif else begin
      
        if arg_present(flats_cube) then begin
          restore, flats_files[0]
          flats_cube = cub
          flats_wav = wav
;        flats_names = namelist ; Is empty?
        endif else begin
          ;; Assume now that if the sav file exists, then the wav file
          ;; and the filenames files also exist. If this is uncertain,
          ;; we need a test and fall back to the sav file if that is all
          ;; that is there.
          flats_wav = readfits(red_strreplace(flats_files[0], 'flats.sav', 'flats_wav.fits'))
        endelse
        
        ;; The namelist variable in the save file was empty in one test
        ;; so read the flats_names from the filenames.txt file instead.
        spawn, 'cat ' + red_strreplace(flats_files[0], 'flats.sav', 'filenames.txt'), flats_names
      endelse
      
    endif

  endif

  if arg_present(cmap) || arg_present(fit_results_struct) || arg_present(fit_results_files) then begin

    fit_results_files = red_uniquify(file_search(search_prefix + '_fit_results.sav') $
                                     , count = Nfiles)

    if Nfiles eq 0 then begin
      red_message, 'No fit results files found. Run fitgains first!'
      retall
    endif

    if arg_present(cmap) || arg_present(fit_results_struct) then begin
      if Nfiles gt 1 then begin
        red_message, 'More than one fit results file found.'
        print, fit_results_files
        red_message, 'No fit results returned.'
        stop
      endif else begin
        restore, fit_results_files[0]
        if arg_present(fit_results_struct) then fit_results_struct = fit
        ;; The cavity map in nm and sign set so lambda_correct =
        ;; lambda + cmap:
        if arg_present(cmap) then cmap = -reform(fit.pars[1,*,*]) / 10.
      endelse
    endif
    
  endif

end

a = crisp2red(/dev, /no)

pref = '6173'
detector = 'camXXXIII'

a -> get_spectral_flats_info, cmap = cmap $
                                     , detector = detector $
                                     , fit_results_struct = fit_results_struct $
                                     , fit_results_files = fit_results_files $
                                     , flats_cube = flats_cube $
                                     , flats_files = flats_files $
                                     , flats_names = flats_names $
                                     , flats_wav = flats_wav $
                                     , prefilter = pref

help, cmap, detector $
      , fit_results_struct  $
      , fit_results_files  $
      , flats_cube  $
      , flats_files  $
      , flats_names  $
      , flats_wav  $
      , pref


end

