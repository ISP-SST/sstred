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
pro crisp::prepflatcubes, flatdir = flatdir $
                            , pref = pref $
                            , verbose = verbose


  ;; Prepare for logging (after setting of defaults).
  ;; Set up a dictionary with all parameters that are in use
  red_make_prpara, prpara, flatdir
  red_make_prpara, prpara, pref
  
  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  ;; Check keywords
  if(~keyword_set(flatdir)) then flatdir = self.out_dir+'/flats/'

  outdir = self.out_dir + '/flats/spectral_flats/'
  file_mkdir, outdir

  ;; Find the narrowband cameras.
  self -> getdetectors
  cams = *self.cameras
  detectors = *self.detectors
  Ncams = n_elements(cams)

  for icam = 0, Ncams-1 do begin

    if strmatch(cams[icam],'*-[DW]') then continue ; Don't do this for WB cameras
    
    ;; Find the files (make sure not to get the cavity free flats!)
    files = file_search(flatdir+'/'+detectors[icam]+'*[0-9].flat.fits', count = Nfiles)

    ;; Check files
    if Nfiles eq 0 then begin
      print, inam+' : No NB files found, returning'
      return
    endif else begin
      print, inam + ' : '+red_stri(Nfiles)+' NB flats found for '+cams[icam]
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

      
      print, inam + ' : working on scans with WB prefilter -> '+upref
      
      sindx = where(states.prefilter eq upref, Nstates)

      if Nstates gt 0 then begin

        sstates = states[sindx]
        
        ;; LC
        ulc = sstates[uniq(sstates.lc, sort(sstates.lc))].lc
        Nlc = n_elements(ulc)
        polarized = Nlc gt 1
        
        if polarized then begin

          pindx = where(sstates.lc eq 0, Nstates)
          sstates = sstates[pindx]
          
          ;; Get polcal data
          pname = self->filenames('polc',sstates)
          pname = pname[uniq(pname, sort(pname))]
          if n_elements(pname) gt 1 then stop
;          pname = file_search(self.out_dir+'/polcal/'+detectors[icam]+'_'+upref+'_polcal.fits', count = Npolcal)
          
          if ~file_test(pname) then begin
;            print, inam + ' : ERROR, '+detectors[icam]+'.'+upref+' -> files not found in '+self.out_dir+'/polcal/'
            print, inam + ' : ERROR, file not found:'
            print, pname
            continue
          endif
          print, inam +' : loading polcal data -> '+pname
          immt = red_readdata(pname)
          immt = red_invert_mmatrix(temporary(immt))
        endif

        ;; Load backscatter data?
        if ~keyword_set(no_descatter) AND (upref eq '8542' or upref eq '7772') then begin
          self -> loadbackscatter, detectors[icam], upref, bg, psf
        endif
        
        for istate = 0L, Nstates - 1 do begin
          
          if istate eq 0  then begin
            head = red_readhead(sstates[istate].filename)
            dim = fxpar(head, 'NAXIS*')
            cub = fltarr([Nstates, dim])
            wav = dblarr(Nstates)
            tomask = bytarr(dim)
          endif

          if polarized then begin
            
            ;; Load flats and demodulate

            fname0 = self->filenames('flat',sstates[istate])           
            ;;flatdir+'/'+detectors[icam]+'_'+sstates[istate].fullstate+'.flat.fits'
                
            lc0 = file_search(fname0, count = nlc0)
            lc1 = file_search(red_strreplace(fname0, 'lc0', 'lc1'), count = nlc1)
            lc2 = file_search(red_strreplace(fname0, 'lc0', 'lc2'), count = nlc2)
            lc3 = file_search(red_strreplace(fname0, 'lc0', 'lc3'), count = nlc3)
            
            if (nlc0 ne 1) or (nlc1 ne 1) or (nlc2 ne 1) or (nlc3 NE 1) then begin
              print, inam+' : ERROR, there is not exactly 1 flat per state!'
              print, inam+' : lc0 -> '+strtrim(nlc0,2)
              print, inam+' : lc1 -> '+strtrim(nlc1,2)
              print, inam+' : lc2 -> '+strtrim(nlc2,2)
              print, inam+' : lc3 -> '+strtrim(nlc3,2)
              stop
            endif
            
            ;; Print info
            print, inam+' : demodulating images: '
            print, '   -> '+lc0
            print, '   -> '+lc1
            print, '   -> '+lc2
            print, '   -> '+lc3
            
            ;; Load data
            lc0 = red_readdata(lc0, /silent)
            lc1 = red_readdata(lc1, /silent)
            lc2 = red_readdata(lc2, /silent)
            lc3 = red_readdata(lc3, /silent)
            
            ;; Descatter ?
            if ~keyword_set(no_descatter) AND (upref EQ '8542' or upref eq '7772') then begin
              lc0 = rdx_descatter(temporary(lc0), bg, psf, /verbose, nthreads = nthreads)
              lc1 = rdx_descatter(temporary(lc1), bg, psf, /verbose, nthreads = nthreads)
              lc2 = rdx_descatter(temporary(lc2), bg, psf, /verbose, nthreads = nthreads)
              lc3 = rdx_descatter(temporary(lc3), bg, psf, /verbose, nthreads = nthreads)
            endif
            
            ;; Demodulate flats
            tmp = reform((red_demodulate_simple(immt, lc0, lc1, lc2, lc3))[*,*,0]) 
            
          endif else begin
                          
            ;; Load non-polarized flats 
            tmp = red_readdata(sstates[istate].filename, /silent)
            if ~keyword_set(no_descatter) AND (upref EQ '8542' or upref eq '7772') then begin
              tmp = rdx_descatter(temporary(tmp), bg, psf, /verbose, nthreads = nthreads)
            endif
            
          endelse
          
          idx = where(~finite(tmp) OR (tmp LT -0.001), nnan)
          if nnan gt 0 then begin
            tmp[idx] = 0.0
            tomask[idx] = 1B
          endif
          
          ;; Result
          cub[istate, *, *] = tmp
          ;;wav[istate] = sstates[istate].tun_wavelength
          wav[istate] = double((strsplit(sstates[istate].tuning,'_',/extract))[1])*1d-13

        endfor                  ; istate

          
        for jj = 0L, dim[1]-1 do for ii = 0L, dim[0]-1 do begin
          if(tomask[ii,jj] eq 1B) then cub[*,ii,jj] = 0.0
        endfor                  ; jj, ii

        ;; Sort states (so far they are sorted as strings -> incorrect order)
        print, inam + ' : sorting wavelengths ... ', FORMAT = '(A,$)'
        ord = sort(wav)
        cub = (temporary(cub))[ord, *,*]
        wav = wav[ord]
        sstates = sstates[ord]
        print, 'done'
        
        ;; Save results (separate fits files for old fortran fitgains_ng)
        firstpart = detectors[icam] + '_' + upref + '_'
        doutname = firstpart + 'flats_data.fits'
        woutname = firstpart + 'flats_wav.fits'
        noutname = firstpart + 'filenames.txt'
        soutname = firstpart + 'flats.sav'
        
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
    endfor                      ; ipref
  endfor                        ; icam

end
