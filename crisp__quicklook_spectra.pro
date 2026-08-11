; docformat = 'rst'

;+
; Plot average spectra for CRISP.
;
; This method will produce one plot per (WB) prefilter of the median
; intensity as a function of wavelength. Colors and symbols are used
; to plot different cameras and LC states in the same diagram.
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
;   regex : in, optional, type=string, default='cam*_lc?.flat.fits'
;
;
; :History:
; 
;   2018-11-19 : MGL. First version.
; 
;-
pro crisp::quicklook_spectra, dir = dir $
                              , outdir = outdir $
                              , regex = regex 

  if n_elements(regex) eq 0 then regex = 'cam*_lc?.flat.fits'
  if n_elements(dir) eq 0 then dir = 'flats'
  indir = red_strreplace(dir, self.out_dir, '')
  indir = red_strreplace(dir, self.out_dir+'/', '')
  if n_elements(outdir) eq 0 then outdir = 'quicklook_spectra/'+indir+'/'
  

  file_mkdir, outdir
  
  files = file_search(indir+'/'+regex, count = Nfiles)

  self -> extractstates, files, states

  upref = states[uniq(states.prefilter, sort(states.prefilter))].prefilter

  for ipref = 0, n_elements(upref)-1 do begin

    pfile = outdir+'/spectrum_' + upref[ipref] + '.pdf'


    self -> selectfiles, files = files, states = states $
                         , pref = upref[ipref] $
                         , selected = sel, count = Nsel
 
    selstates = states[sel]
    selfiles = files[sel]

    intensity = fltarr(Nsel)

    for isel = 0, Nsel-1 do begin
      frame = red_readdata(selfiles[isel])
      intensity[isel] = median(frame)
    endfor

    ulc = selstates[uniq(selstates.lc, sort(selstates.lc))].lc
    Nlc = n_elements(ulc) <4
    colors = ['forest green', 'red', 'blue', 'dark goldenrod']
    ulc = ulc[0:Nlc-1]
    colors = colors[0:Nlc-1]
    
    ucam = selstates(uniq(selstates.camera, sort(selstates.camera))).camera
    Ncams = n_elements(ucam)
    symbols = [1, 7]
    
    cgwindow
    xrange = [min(selstates.tun_wavelength*1e9) $
              , max(selstates.tun_wavelength*1e9)]
    margin = (xrange[1]-xrange[0]) * 0.05

    cgplot, /add, [0], [0], /nodata $
            , xrange = xrange + margin*[-1, 1] $
            , yrange = [0, max(intensity)]*1.05 $
            , xtitle = '$\lambda$/1 nm', ytitle = 'intensity / 1 D.N.' $
            , title = self.isodate + ' ' + indir + ' ' + upref[ipref]

    for icam = 0, Ncams-1 do begin
      for ilc = 0, Nlc-1 do begin
        indx = where(selstates.lc eq ulc[ilc] and selstates.camera eq ucam[icam])
        cgplot, /add, /over $
                , selstates[indx].tun_wavelength*1e9, intensity[indx] $
                , psym = symbols[icam] $
                , color = colors[ilc]
      endfor                    ; ilc
    endfor                      ; icam
    
    
    cglegend, /add, titles = ucam $
              , color = 'black', psym = symbols $
              , length = 0, align = 2, location = [.3, .1] $
              , vspace = 2.

    if Nlc gt 1 then cglegend, /add, titles = 'LC'+strtrim(long(ulc),2) $
                               , psym = 16, color = colors $
                               , length = 0, align = 3, location = [.8, .1] $
                               , vspace = 2.

    cgcontrol, out = pfile
    
  endfor                        ; ipref

  
end


a = crispred(/dev)
a -> quicklook_spectra


end

