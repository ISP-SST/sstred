; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CRISP pipeline
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
;    fixcav  : 
;   
;   
;   
;    w0  : 
;   
;   
;   
;    w1  : 
;   
;   
;   
;    pref  : 
;   
;   
;   
;    noasy  : 
;   
;   
;   
;    shift  : 
;   
;   
;   
; 
; 
; :History:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2013-07-11 : MGL. Use red_satlas rather than satlas. Use
;                red_intepf rather than intepf.
; 
;   2013-08-27 : MGL. Added support for logging. Let the subprogram
;                find out its own name.
;
;   2015-05-05 : MGL. If keyword pref is given, use it when looking
;                for files.
;
;   2016-03-22 : JLF. Added support for lc4 flat cubes (see
;                modifications to red::prepflatcubes_lc4). Changed
;                'pref' keyword to 'state' and made its behavior
;                similar to red::fitgains. Added some error checking.
; 
;   2016-04-01 : THI. Changed 'state' keyword back to 'pref' for now,
;                so that old scripts do not need to be modified.
;
;   2016-10-13 : MGL. Adapted to CHROMIS and new pipeline mechanisms. 
;
;   2016-10-27 : MGL. Special case for prefilters with only a single
;                wavelength point.
; 
; 
;-
pro red::fitprefilter, fixcav = fixcav $
                       , w0 = w0, w1 = w1 $
                       , pref = pref $
                       , noasy = noasy $
                       , shift = shift $
                       , init=init $
                       , stretch=stretch $
                       , weight = weight

  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  ;; Logging
  help, /obj, self, output = selfinfo 
  red_writelog, selfinfo = selfinfo

  ;; Load idlsave file with the results from fitgains

  ;; file =  self.out_dir + '/flats/spectral_flats/' + cam + '.fit_results.sav'
  ;; file1 = self.out_dir + '/flats/spectral_flats/' + cam + '_flats.sav'
  ;; if(~file_test(file)) then begin
  ;;    print, inam + ' : ERROR, file not found -> ' + file
  ;;    return
  ;; endif
  ;; if(~file_test(file1)) then begin
  ;;    print, inam + ' : ERROR, file not found -> ' + file1
  ;;    return
  ;; endif

  if n_elements(pref) eq 0 then begin
    sfiles = file_search(self.out_dir + '/flats/spectral_flats/*_fit_results.sav' $
                         , count=Nsfiles)
    xfiles = file_search(self.out_dir + '/flats/spectral_flats/*_flats.sav' $
                         , count=Nxfiles)
  endif else begin
    sfiles = file_search(self.out_dir +'/flats/spectral_flats/*' + pref $
                         + '*_fit_results.sav', count=Nsfiles)
    xfiles = file_search(self.out_dir +'/flats/spectral_flats/*' + pref $
                         + '*_flats.sav', count=Nxfiles)
  endelse

  print, sfiles
  print
  print, xfiles
  
  if Nsfiles eq 0 then begin
     message, 'No fit results files found. Run the fitgains method first!'
     return
  endif

  ;; Select data
  selectionlist = file_basename(xfiles,'_flats.sav') 
  if keyword_set(all) or Nsfiles eq 1 then begin
    Nselect = Nsfiles
    sindx = indgen(Nselect)
  endif else begin
    tmp = red_select_subset(selectionlist $
                            , qstring = 'Select data to be processed' $
                            , count = Nselect, indx = sindx)
  endelse
  print, inam + ' : Will process the following data:'
  print, selectionlist[sindx], format = '(a0)'

  for iselect = 0L, Nselect-1 do begin

    isfile = sindx[iselect]
    
    xfile = xfiles[isfile]                                      ; Input flats results 
    sfile = strreplace(xfile, '_flats.sav', '_fit_results.sav') ; Input fit results 
    nfile = strreplace(xfile, '_flats.sav', '_filenames.txt')   ; Input flat filenames
;    ffile = strreplace(xfile, '_flats.sav', '_flats_data.fits') ; Input flat intensities
;    wfile = strreplace(xfile, '_flats.sav', '_flats_wav.fits')  ; Input flat wavelengths
    
    if ~file_test(nfile) then begin
      print, inam+' : Cannot find file '+nfile
      continue
    endif
;    if ~file_test(ffile) then begin
;      print, inam+' : Cannot find file '+ffile
;      continue
;    endif
;    if ~file_test(wfile) then begin
;      print, inam+' : Cannot find file '+wfile
;      continue
;    endif


    ;; Load data
    print, file_basename(nfile, '_filenames.txt')
    spawn, 'cat ' + nfile, namelist
    Nwav = n_elements(namelist)

    self -> extractstates, namelist, states

    ;;camera = states[0].camera
    detector = states[0].detector

    ;; Get prefilter
    prefs = states[0].prefilter
;    prefs = (strsplit(file_basename(fit.oname[0]), '.',/extract))[1]
    print, inam + ' : Processing prefilter ' + prefs
;    dpr = double(prefs)
    ;; This is the reference point of the fine tuning:
    dpr = double((strsplit(states[0].tuning,'_',/extract))[0])

    if ~file_test(sfile) then begin
       print, inam+' : Cannot find file '+sfile
       continue
    endif

    restore, sfile
    restore, xfile
 
    fac = median(fit.pars[0,*,*])
    fit.yl *= fac
    if(~keyword_set(w0)) then ww0 = 0
    if(~keyword_set(w1)) then ww1 = n_elements(wav) - 1

    ;; Load satlas
    red_satlas, min(fit.xl) + dpr - 1.0, max(fit.xl) + dpr + 1.0, xs, ys
    xs -= dpr
    
    ;; Get FPI transmission profile
    fpi = red_get_fpi_par(line = prefs)
    dw = xs[1] - xs[0]
    np = long((max(xs) - min(xs)) / dw) - 2
    if(np/2*2 eq np) then np -= 1L
    tw = (dindgen(np) - np/2) * dw
    tr = red_get_fpi_trans(fpi, tw + fpi.w0, ecl = 0.0, ech = 0.0, erh = -0.01)

    dum = max(tr, p)
    cc = poly_fit(tw[p-1:p+1] * 100.d0, tr[p-1: p+1], 2)
    off = -0.005d0 * cc[1] / cc[2]
    tr = red_get_fpi_trans(fpi, tw + fpi.w0 + off, ecl = 0.0, ech = 0.0, erh = -0.01)

    ;; Convolve Solar Atlas with the FPI transmission profile
    ys = red_convl(ys, tr, /usefft)

    if Nwav gt 1 then begin

      ;; Usual case with several wavelength points.
       
      ;; Pack variables for mpfit                                
      yl1 = red_intepf(fit.xl, fit.yl, wav[ww0:ww1])

      if(n_elements(weight) eq 0) then weight=dblarr(n_elements(fit.yl))+1.0d0
      
      mm = {xl:fit.xl, yl:fit.yl, wav:wav[ww0:ww1], yl1:yl1}
      functargs = {xs:xs, ys:ys, dpr:dpr, mm:mm, w:weight}
      
      ;; Init guess model
      pp = dblarr(8)
      if(n_elements(init) gt 0) then pp[0]=init else begin
        pp[0] = fac             ; Scale factor
        pp[1] = -0.5d0          ; Pref. shift
        pp[2] = 6.0d0           ; Pref. FWHM
        pp[3] = 2.d0            ; Pref. ncav
        pp[4] = -0.001d0        ; Line shift (satlas-obs)
        pp[5] = 0.0001d0
        pp[6] = 0.0001d0
        pp[7] = 1.0d0
      endelse

      fitpars = replicate({mpside:2, limited:[0,0], limits:[0.0d, 0.0d], fixed:0}, 8)
          
      ;; Limits and settings for the various model parameters

      ;; Prefilter FWHM
      fitpars[2].limited = [1,1]
      fitpars[2].limits = [1.d0, 10.2d0]

      ;; Filter cavities
      if keyword_set(fixcav) then begin
        pp[3] = fixcav
        fitpars[3].fixed = 1
      endif else begin
        pp[3] = 2.0d0           ; Default 2 cavities
        fitpars[3].limits = [1.4d0, 2.4d0]
        fitpars[3].limited = [1,1]
      endelse

      ;; Solar atlas wavelength shift
      if(keyword_set(shift)) then begin
        pp[4] = shift
      endif

      ;; Asymmetry
      if(keyword_set(noasy)) then begin
        fitpars[5].fixed = 1B
        fitpars[6].fixed = 1B
        pp[5] = 0.
        pp[6] = 0.
      endif
      fitpars[6].fixed = 1B

      ;; Wavelength stretch
      if(~keyword_set(stretch)) then fitpars[7].FIXED = 1B


      print, inam + ' : Initial model parameters:'
      print, inam + ' : p[0] -> ', pp[0], ' (scale factor)'
      print, inam + ' : p[1] -> ', pp[1], ' (prefilter shift)'
      print, inam + ' : p[2] -> ', pp[2], ' (prefilter FWHM)'
      print, inam + ' : p[3] -> ', pp[3], ' (prefilter number of cavities)'
      print, inam + ' : p[4] -> ', pp[4], ' (solar atlas shift)'
      print, inam + ' : p[5] -> ', pp[5], ' (asymmetry term 1)'
      print, inam + ' : p[6] -> ', pp[6], ' (asymmetry term 2)'
      print, inam + ' : p[7] -> ', pp[7], ' (Wavelength stretch)'

      ;; Call mpfit
      pp = mpfit('red_fit_prefilter', pp, functargs = functargs, parinfo = fitpars, /quiet)
      dum = red_fit_prefilter(pp, xs = xs, ys = ys, dpr = dpr, mm = mm, pref = pref, w=weight)

      print, inam + ' : Fitted model parameters:'
      print, inam + ' : p[0] -> ', pp[0], ' (scale factor)'
      print, inam + ' : p[1] -> ', pp[1], ' (prefilter shift)'
      print, inam + ' : p[2] -> ', pp[2], ' (prefilter FWHM)'
      print, inam + ' : p[3] -> ', pp[3], ' (prefilter number of cavities)'
      print, inam + ' : p[4] -> ', pp[4], ' (solar atlas shift)'
      print, inam + ' : p[5] -> ', pp[5], ' (asymmetry term 1)'
      print, inam + ' : p[6] -> ', pp[6], ' (asymmetry term 2)'
      print, inam + ' : p[7] -> ', pp[7], ' (Wavelength stretch)'

    endif else begin

      ;; Special treatment of CaH-cont with only a single wavelength
      ;; point.

      ;; Pack variables for mpfit                                
      if fit.xl[0] ne wav[0] then stop
      yl1 = fit.yl
      
      if n_elements(weight) eq 0 then weight = [1.0d0]
      
      mm = {xl:fit.xl, yl:fit.yl, wav:wav[ww0:ww1], yl1:yl1}
      functargs = {xs:xs, ys:ys, dpr:dpr, mm:mm, w:weight}
      
      ;; Init guess model
      Nparam = 8
      pp = dblarr(Nparam)
      if(n_elements(init) gt 0) then pp[0]=init else begin
        pp[0] = fac             ; Scale factor
        pp[1] = 0d;-0.5d0          ; Pref. shift
        pp[2] = 6.0d0           ; Pref. FWHM
        pp[3] = 2.d0            ; Pref. ncav
;        pp[4] = -0.001d0        ; Line shift (satlas-obs)
;        pp[5] = 0.0001d0
;        pp[6] = 0.0001d0
        pp[7] = 1.0d0
      endelse

      fitpars = replicate({mpside:2 $
                           , limited:[0,0] $
                           , limits:[0.0d, 0.0d] $
                           , fixed:0} $
                          , Nparam)
      
      ;; Limits and settings for the various model parameters

      fitpars[1:7].fixed = 1

;      ;; Prefilter FWHM
;      fitpars[2].limited = [1,1]
;      fitpars[2].limits = [1.d0, 10.2d0]
;
;      ;; Filter cavities
;;      If(keyword_set(fixcav)) then begin
;;        pp[3] = fixcav
;        fitpars[3].fixed = 1
;;      endif else begin
;;        pp[3] = 2.0d0           ; Default 2 cavities
;;        fitpars[3].limits = [1.4d0, 2.4d0]
;;        fitpars[3].limited = [1,1]
;;      endelse
;
;      ;; Asymmetry 2
;      fitpars[6].fixed = 1
;
;      ;; Wavelength stretch
;      if(~keyword_set(stretch)) then fitpars[7].FIXED = 1B
;
;;      if(keyword_set(noasy)) then begin
;        fitpars[5].fixed = 1
;        fitpars[6].fixed = 1
;        pp[5] = 0.
;        pp[6] = 0.
;;      endif
;      if(keyword_set(shift)) then begin
;        pp[4] = shift
;      endif

      print, inam + ' : Initial model parameters:'
      print, inam + ' : p[0] -> ', pp[0], ' (scale factor)'
      print, inam + ' : p[1] -> ', pp[1], ' (prefilter shift)'
      print, inam + ' : p[2] -> ', pp[2], ' (prefilter FWHM)'
      print, inam + ' : p[3] -> ', pp[3], ' (prefilter number of cavities)'
      print, inam + ' : p[4] -> ', pp[4], ' (solar atlas shift)'
      print, inam + ' : p[5] -> ', pp[5], ' (asymmetry term 1)'
      print, inam + ' : p[6] -> ', pp[6], ' (asymmetry term 2)'
      print, inam + ' : p[7] -> ', pp[7], ' (Wavelength stretch)'

      ;; Call mpfit
      pp = mpfit('red_fit_prefilter', pp, functargs = functargs, parinfo = fitpars, /quiet)
      dum = red_fit_prefilter(pp, xs = xs, ys = ys, dpr = dpr, mm = mm, pref = pref, w=weight)

      print, inam + ' : Fitted model parameters:'
      print, inam + ' : p[0] -> ', pp[0], ' (scale factor)'
      print, inam + ' : p[1] -> ', pp[1], ' (prefilter shift)'
      print, inam + ' : p[2] -> ', pp[2], ' (prefilter FWHM)'
      print, inam + ' : p[3] -> ', pp[3], ' (prefilter number of cavities)'
      print, inam + ' : p[4] -> ', pp[4], ' (solar atlas shift)'
      print, inam + ' : p[5] -> ', pp[5], ' (asymmetry term 1)'
      print, inam + ' : p[6] -> ', pp[6], ' (asymmetry term 2)'
      print, inam + ' : p[7] -> ', pp[7], ' (Wavelength stretch)'


      ;; Fit the single point to the convolved spectrum.
      
      ;; Must make the data corresponding to pref, mm.wav, and pp, to
      ;; be written to file.

    endelse
    
    odir = self.out_dir + '/prefilter_fits/'
    file_mkdir, odir
    
    ofile = detector + '.'+prefs+'.prefilter.f0'
    print, inam + ' : saving prefilter to file -> ' + odir + ofile
    fzwrite, float(pref), odir + ofile, ' '
    ofile = detector + '.'+prefs+'.prefilter_wav.f0'
    fzwrite, float(mm.wav), odir + ofile, ' '
    
    ofile = detector + '.'+prefs+'.prefilter_pars.f0'
    print, inam + ' : saving fit-results to file -> ' + odir + ofile
    fzwrite, pp, odir + ofile, ' '
    
  endfor                         ; iselect
  
end
