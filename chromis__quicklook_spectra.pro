; docformat = 'rst'

;+
; Plot average spectra for CHROMIS.
;
; This method will produce one plot per (WB) prefilter of the median
; intensity as a function of wavelength. 
;
; To use this for science data, it should probably be rewritten so it
; can do something intelligent with multiple scans.
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
; :Keywords:
; 
;   dir : in, type=string, default='flat'
; 
;     The input directory in which to find the data.
; 
;   outdir : in, optional, type=string, default="Constructed from the input directory"
;
;     The directory in which to write the plots.
;
;   regex : in, optional, type=string, default='cam*[0-9].flat.fits'
;
;
; :History:
; 
;   2018-11-19 : MGL. First version.
; 
;-
pro chromis::quicklook_spectra, dir = dir $
                                , outdir = outdir $
                                , regex = regex 

  if n_elements(regex) eq 0 then regex = 'cam*[0-9].flat.fits'
  if n_elements(dir) eq 0 then dir = 'flats'
  indir = red_strreplace(dir, self.out_dir, '')
  indir = red_strreplace(dir, self.out_dir+'/', '')
  if n_elements(outdir) eq 0 then outdir = 'quicklook_spectra/'+indir+'/'
  

  file_mkdir, outdir
  
  files = file_search(indir+'/'+regex, count = Nfiles)

  self -> extractstates, files, states

  ;; For CHROMIS we need to select first based on WB prefilter, then
  ;; get the narrowband files with the same FPI_STATE. Alternatively,
  ;; somewhere the pipeline knows what WB prefilter that any NB
  ;; prefilter is combined with. So we can select for that.
  indx = where(strmatch(states.camera,'*-N'))
  states = states[indx]
  files = files[indx]
  
  upref = states[uniq(states.prefilter, sort(states.prefilter))].prefilter
  Nprefs = n_elements(upref)

  ucam = states(uniq(states.camera, sort(states.camera))).camera
  Ncams = n_elements(ucam)
  if Ncams gt 1 then stop
  symbols = [1]

  
  for ipref = 0, Nprefs-1 do begin
 
    pfile = outdir+'/spectrum_' + upref[ipref] + '.pdf'

    self -> selectfiles, files = files, states = states $
                         , pref = upref[ipref] $
                         , selected = sel, count = Nsel

    if Nsel gt 0 then begin
 
      selstates = states[sel]
      selfiles = files[sel]

      intensity = fltarr(Nsel)

      for isel = 0, Nsel-1 do begin
        frame = red_readdata(selfiles[isel])
        intensity[isel] = median(frame)
      endfor

      cgwindow
      xrange = [min(selstates.tun_wavelength*1e9) $
                , max(selstates.tun_wavelength*1e9)]
      margin = (xrange[1]-xrange[0]) * 0.05 >0.05

      cgplot, /add, [0], [0], /nodata $
              , xrange = xrange + margin*[-1, 1] $
              , yrange = [0, max(intensity)]*1.05 $
              , xtitle = '$\lambda$/1 nm', ytitle = 'intensity / 1 D.N.' $
              , title = self.isodate + ' ' + indir + ' ' + upref[ipref]

      cgplot, /add, /over $
              , [selstates.tun_wavelength*1e9], [intensity] $
              , psym = symbols[0]
      
      cglegend, /add, titles = ucam $
                , color = 'black', psym = symbols $
                , length = 0, align = 2, location = [.33, .1] $
                , vspace = 2.

      cgcontrol, out = pfile
      

    endif
  
  endfor                        ; ipref

end

a = chromisred(/dev)
a -> quicklook_spectra



end

