; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CHROMIS pipeline
; 
; 
; :Author:
; 
; 
; 
; 
; :Returns:
; 
; 
; :Params:
; 
; 
; :Keywords:
; 
;    flatdir  : 
;   
;   
;   
;    verbose  : 
;   
;   
;    pref : in, optional, type=string
; 
;       Select a narrowband prefilter.
; 
; :History:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2013-08-27 : MGL. Added support for logging. Let the subprogram
;                find out its own name.
; 
;   2016-02-15 : MGL. Use loadbackscatter. Remove keyword descatter,
;                new keyword no_descatter.
;
;   2016-03-22 : JLF. Fixed a bug in which red::prepflatcubes_lc4 would
;		 overwrite the results of red::prepflatcubes if there are
;		 both lc4 and lc0-3 datasets at the same prefilter. If
;		 lc0-3 datasets are detected it will output files with
;		 .lc4 added to the state name.
; 
;   2016-05-10 : THI. New keyword pref.
;
;   2016-08-23 : THI. Rename camtag to detector and channel to camera,
;                so the names match those of the corresponding SolarNet
;                keywords.
;
;   2016-10-04 : MGL. Adapted red::preflatcubes_lc4 for CHROMIS data.
;
;   2017-04-13 : MGL. Do not read cavityfree flats! Make FITS header
;                for the cube.
;
;-
pro chromis::prepflatcubes, flatdir = flatdir $
                            , pref = pref $
                            , verbose = verbose


  ;; Prepare for logging (after setting of defaults).
  ;; Set up a dictionary with all parameters that are in use
  red_make_prpara, prpara, flatdir
  red_make_prpara, prpara, pref
  
  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

;  ;; Logging
;  help, /obj, self, output = selfinfo 
;  red_writelog, selfinfo = selfinfo

  ;; Check keywords
  if(~keyword_set(flatdir)) then flatdir = self.out_dir+'/flats/'

  outdir = self.out_dir + '/flats/spectral_flats/'
  file_mkdir, outdir

  ;; Find the narrowband camera.
  cams = *self.cameras
  isnb = strmatch(*self.cameras,'*-N')
  if total(isnb) eq 0 then return  
  nbcam = (*self.cameras)[where(isnb)]

  self -> getdetectors
  detector = (*self.detectors)[where(isnb)]

  ;; Find the files (make sure not to get the cavity free flats!)
  files = file_search(flatdir+'/'+detector+'*[0-9].flat.fits', count = Nfiles)


  ;; Check files
  if Nfiles eq 0 then begin
    print, inam+' : No NB files found, returning'
    return
  endif else begin
    print, inam + ' : '+red_stri(Nfiles)+' NB flats found'
  endelse

  self -> extractstates, files, states

  ;; Prefilters
  uprefs = states[uniq(states.prefilter, sort(states.prefilter))].prefilter
  Nprefs = n_elements(uprefs)

  ;; Loop nb prefilters 
  for ipref = 0L, Nprefs - 1 do begin

    if keyword_set(pref)  then begin
      if uprefs[ipref] ne pref then begin
        print, inam + ' : skipping prefilter -> '+uprefs[ipref]
        continue
      endif
    endif 

    upref = uprefs[ipref]

    print, inam + ' : working on scans with wb prefilter -> '+upref

    ;;    pos = where((states.prefilter eq upref), Nstates)

    ;; Flats come with different camera settings
    usettings = states[uniq(states.cam_settings, sort(states.cam_settings))].cam_settings
    Nsettings = n_elements(usettings)
    print, inam + ' : '+red_stri(Nsettings)+' different camera settings found'
    
    for isetting = 0L, Nsettings-1 do begin

      ;; Select states with the current settings
      usetting = usettings[isetting]
      sindx = where(states.prefilter eq upref and states.cam_settings eq usetting, Nstates)
;        Nstates = n_elements(sindx)
      if Nstates gt 0 then begin
        sstates = states[sindx]

        for istate = 0L, Nstates - 1 do begin
          
          if(istate eq 0) then begin
            head = red_readhead(sstates[istate].filename)
            dim = [fxpar(head, 'NAXIS1'), fxpar(head, 'NAXIS2')]
            cub = fltarr([Nstates, dim])
            wav = dblarr(Nstates)
            tomask = bytarr(dim)
          endif

          ;; Load flats 
          tmp = red_readdata(sstates[istate].filename)
          idx = where(~finite(tmp) OR (tmp LT -0.001), nnan)

          if nnan gt 0 then begin
            tmp[idx] = 0.0
            tomask[idx] = 1B
          endif

          ;; Result
          cub[istate, *, *] = tmp
;           wav[istate] = stat.wav[pos[istate]] - double(stat.pref[pos[istate]])
;           wav[istate] = sstates[istate].tun_wavelength
          wav[istate] = double((strsplit(sstates[istate].tuning,'_',/extract))[1])*1d-13

        endfor                  ; istate

        for jj = 0L, dim[1]-1 do for ii = 0L, dim[0]-1 do begin
          if(tomask[ii,jj] eq 1B) then cub[*,ii,jj] = 0.0
        endfor

        ;; Sort states (so far they are sorted as strings -> incorrect order)
        print, inam + ' : sorting wavelengths ... ', FORMAT = '(A,$)'
        ord = sort(wav)
        cub = (temporary(cub))[ord, *,*]
        wav = wav[ord]
        sstates = sstates[ord]
        print, 'done'

        
        ;; Save results (separate fits files for old fortran fitgains_ng)
        doutname = detector + '_' + usetting + '_' + upref + '_flats_data.fits'
        woutname = detector + '_' + usetting + '_' + upref + '_flats_wav.fits'
        noutname = detector + '_' + usetting + '_' + upref + '_filenames.txt'
        soutname = detector + '_' + usetting + '_' + upref + '_flats.sav'
        
        ;; Make FITS header
        check_fits, cub, head, /UPDATE, /SILENT  
        fxaddpar, head, 'DATE', red_timestamp(/iso), 'UTC creation date of FITS header'
        fxaddpar, head, 'FILENAME', file_basename(doutname[0]), after = 'DATE'
        self -> headerinfo_addstep, head, prstep = 'Flat cubes' $
                                    , prproc = inam, prpara = prpara

        ;; More things need to be added...
        ;; * What's on the axes?
        

        writefits, outdir + doutname, cub, head

        ;; Should make header also for the wavelength file. (Or add
        ;; the wavelengths to the cube file, but reader program has to
        ;; deal with that!)
        writefits, outdir + woutname, wav
        namelist = strarr(Nstates)

        ;; Print file names in order into a text file. Note that this
        ;; is no longer a list of file names for the cavity error
        ;; corrected flats (produced by the fitgains method), it's the
        ;; filenames of the input (average) flats.
        openw, lun, outdir + noutname, /get_lun
        printf, lun, sstates.filename, format = '(a0)'
        free_lun, lun

        ;; Save as structure for new routines (red::fitgain_ng)
        print, inam + ' : saving -> ' + soutname
        
        save, file = outdir+soutname, cub, wav, namelist
        
      endif                     ; Nstates?      
    endfor                      ; isetting
  endfor                        ; ipref

end
