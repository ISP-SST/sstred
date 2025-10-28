; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    SST pipeline
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
;   2022-08-01 : MGL. New version for CRISP2 (and old CRISP with new
;                cameras) based on the old CRISP version.
;
;   2024-11-02 : JdlCR. Modifications for new
;                demodulation/flat-fielding scheme
; 
;   2025-10-27 : MGL. Exclude polcal flats.
;
;-
pro red::prepflatcubes, flatdir = flatdir $
                        , pref = pref $
                        , verbose = verbose $
                        , nthreads = nthreads $
                        , no_descatter = no_descatter


  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)  

  ;; Check keywords
  if(~keyword_set(flatdir)) then flatdir = self.out_dir+'/flats/'

  outdir = self.out_dir + '/flats/spectral_flats/'
  file_mkdir, outdir

  ;; Find the narrowband cameras.
  self -> getdetectors
  cams = *self.cameras
  detectors = *self.detectors
  Ncams = n_elements(cams)

  ;; Descatter only needed for CRISP with old Sarnoff cameras,
  ;; processed with the CRISP class.
  if ((typename(self)).tolower()) ne 'crisp' then no_descatter = 1 
  
  ;; Prepare for logging (after setting of defaults).
  ;; Set up a dictionary with all parameters that are in use
  red_make_prpara, prpara, flatdir
  red_make_prpara, prpara, pref
  red_make_prpara, prpara, no_descatter
  
  for icam = 0, Ncams-1 do begin

    if strmatch(cams[icam],'*-[DW]') then continue ; Don't do this for WB cameras
    
    ;; Find the files (make sure not to get the cavity free flats!)
    files = file_search(flatdir+'/'+detectors[icam]+'_*[0-9].flat.fits', count = Nfiles)

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
          
          if ~file_test(pname) then begin
            print, inam + ' : ERROR, file not found:'
            print, pname
            continue
          endif

          polcal_flatfielding = mrdfits(pname,'PF')
          
          if(~polcal_flatfielding) then begin
            print, inam +' : loading polcal data -> '+pname
                      
            immt = red_readdata(pname)
            immt = red_invert_mmatrix(temporary(immt))
          endif


          if polcal_flatfielding then begin
            ;; We probably want to exclude polcal flats from the
            ;; flatcube, as it is usually far outside the line. This
            ;; requires checking the polcal tuning against the tunings
            ;; of the flats, and - if there are matches - against all
            ;; science data.

            ;; The polcal cube (polc) has no tuning info so check the
            ;; polcal sums.
            pfiles = file_search('polcal_sums/*-R/cam*_'+upref+'_*pols.fits', count = Npolcal)
            if Npolcal gt 0 then begin
              phdr = headfits(pfiles[0])
              pstate = fxpar(phdr, 'STATE')
              red_extractstates, pstate, wav = ptuning 
              tmp = where(sstates.tuning eq ptuning[0], Nmatch, complement = kindx)
              if Nmatch gt 0 then begin
                ;; There are flats for the polcal tuning. We want to
                ;; exclude them from the flatcube if there are no
                ;; science data with that tuning.
                e_such_data = 0 ; Assume there are no such science data until we find them
                red_message, 'Checking for science data with polcal tuning ('+ptuning+').'
                for idir = 0L, n_elements(*self.data_dirs)-1 do begin
                  ;; The raw_search method does not seem to support
                  ;; filtering on tuning at the moment so just search
                  ;; for scan 0 and check the tunings.
                  dfiles = self -> raw_search((*self.data_dirs)[idir] + '/' $
                                              + ((typename(self)).tolower()).capwords() + '-R/' $
                                              , scannos = 0)
                  self -> extractstates, dfiles, dstates
                  Nmatch = round(total(dstates.tuning eq ptuning[0] $
                                       and self -> match_prefilters(upref, dstates.prefilter)))
                  ;; If Nmatch is non-zero we do have science data
                  ;; with the polcal tuning. So we do not want to
                  ;; exclude the tuning from the flatcube.
                  e_such_data OR= (Nmatch gt 0)
                  if e_such_data then break ; No need to check more data
                endfor          ; idir
                if ~e_such_data then begin
                  red_message, 'Found no science data with polcal tuning, removing polcal flats from the cubes.'
                  sstates = sstates[kindx]
                  Nstates = n_elements(kindx)
                endif
              endif             ; Nmatch
            endif               ; Npolcal
          endif                 ; polcal_flatfielding
          
        endif                   ; polarized

        ;; Load backscatter data?
        if ~keyword_set(no_descatter) AND self.dodescatter AND (upref eq '8542' or upref eq '7772') then begin
          self -> loadbackscatter, detectors[icam], upref, bg, psf
        endif
        
        for istate = 0L, Nstates - 1 do begin
          
          if istate eq 0  then begin
            head = red_readhead(sstates[istate].filename)
            dim = fxpar(head, 'NAXIS*')
            cub = fltarr([Nstates, min([Nlc, 4]), dim])
            wav = dblarr(Nstates)
            tomask = bytarr(dim)
          endif

          if polarized then begin
            
            ;; Load flats and demodulate

            fname0 = self->filenames('flat',sstates[istate])           
            
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
            print, inam+' : processing images: '
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
            if ~keyword_set(no_descatter) AND self.dodescatter AND (upref EQ '8542' or upref eq '7772') then begin
              lc0 = rdx_descatter(temporary(lc0), bg, psf, /verbose, nthreads = nthreads)
              lc1 = rdx_descatter(temporary(lc1), bg, psf, /verbose, nthreads = nthreads)
              lc2 = rdx_descatter(temporary(lc2), bg, psf, /verbose, nthreads = nthreads)
              lc3 = rdx_descatter(temporary(lc3), bg, psf, /verbose, nthreads = nthreads)
            endif


            if(~polcal_flatfielding) then begin
              ;; Demodulate flats
              tmp = reform((red_demodulate_simple(immt, lc0, lc1, lc2, lc3))[*,*,0]) 
            endif else begin

              dim = size(lc0, /dim)
              tmp = fltarr(4,dim[0], dim[1])
              
              tmp[0,*,*] = lc0
              tmp[1,*,*] = lc1
              tmp[2,*,*] = lc2
              tmp[3,*,*] = lc3
            endelse
            
          endif else begin
                          
            ;; Load non-polarized flats 
            tmp = red_readdata(sstates[istate].filename, /silent)
            if ~keyword_set(no_descatter) AND self.dodescatter AND (upref EQ '8542' or upref eq '7772') then begin
              tmp = rdx_descatter(temporary(tmp), bg, psf, /verbose, nthreads = nthreads)
            endif
            
          endelse

          tmp1 = total(tmp,1)*0.25
          idx = where(~finite(tmp1) OR (tmp1 LT -0.001), nnan)
          if nnan gt 0 then begin
            tomask[idx] = 1B
          endif
          
          ;; Result
          cub[istate, *, *, *] = tmp
          ;;wav[istate] = sstates[istate].tun_wavelength
          wav[istate] = double((strsplit(sstates[istate].tuning,'_',/extract))[1])*1d-13

        endfor                  ; istate


       
        for jj = 0L, dim[1]-1 do for ii = 0L, dim[0]-1 do begin
          if(tomask[ii,jj] eq 1B) then cub[*,*,ii,jj] = 0.0
        endfor                  ; jj, ii
        
        ;; Sort states (so far they are sorted as strings -> incorrect order)
        print, inam + ' : sorting wavelengths ... ', FORMAT = '(A,$)'
        ord = sort(wav)
        cub = (temporary(cub))[ord,*,*,*]
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
        self -> headerinfo_addstep, head, prstep = 'CONCATENATE' $
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
