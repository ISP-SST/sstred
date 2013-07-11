; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :author:
; 
; 
; 
; 
; :returns:
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
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2013-07-11 : MGL. Use red_satlas rather than satlas. Use
;                red_intepf rather than intepf.
; 
; 
; 
;-
pro red::fitprefilter,  fixcav = fixcav, w0 = w0, w1 = w1, pref = pref, noasy = noasy, shift = shift
  inam = 'red::fitprefilter : '
  
                                ;
  ;; Load idlsave file with the results from fitgains_ng
                                ;
  ;; file =  self.out_dir + '/flats/spectral_flats/' + cam + '.fit_results.sav'
  ;; file1 = self.out_dir + '/flats/spectral_flats/' + cam + '_flats.sav'
  ;; if(~file_test(file)) then begin
  ;;    print, inam + 'ERROR, file not found -> ' + file
  ;;    return
  ;; endif
  ;; if(~file_test(file1)) then begin
  ;;    print, inam + 'ERROR, file not found -> ' + file1
  ;;    return
  ;; endif
  
  files = file_search(self.out_dir + '/flats/spectral_flats/*.fit_results.sav', count=count)
  files1 = file_search(self.out_dir + '/flats/spectral_flats/*.flats.sav', count=count1)
  
  ;; Select file
   
  stat = strarr(count)
  stat1 = strarr(count1)
  idd = intarr(count)
   
  for ii = 0, count-1 do stat[ii] = strjoin((strsplit(file_basename(files[ii]),'.',/extract))[0:1],'.')
  for ii = 0, count1-1 do stat1[ii] = strjoin((strsplit(file_basename(files1[ii]),'.',/extract))[0:1],'.')
   
  print, inam + 'found valid states:'
  k =-1
  for ii = 0, count -1 do begin
     dum = where(stat1 eq stat[ii], cc)
     if(cc eq 0) then begin
        continue
     endif
     idd[ii] = dum
     k += 1
     print, ii,' -> ' + stat[ii]
  endfor
  idx = 0
  if(k gt 0) then read, idx, prompt = inam + 'select state: '
  file = files[idx]
  file1 = files1[idd[idx]]
  cam = (strsplit(stat[idx],'.',/extract))[0]

  restore, file
  restore, file1
  fac = median(fit.pars[0,*,*])
  fit.yl *= fac
  if(~keyword_set(w0)) then w0 = 0
  if(~keyword_set(w1)) then w1 = n_elements(wav) - 1
  
  ;; Get prefilter
   
  prefs = (strsplit(file_basename(fit.oname[0]), '.',/extract))[1]
  print, inam + 'Processing prefilter at ' + prefs
  dpr = double(prefs)

  ;; load satlas
   
  red_satlas, min(fit.xl) + dpr - 1.0, max(fit.xl) + dpr + 1.0, xs, ys
  xs -= dpr
  
  ;; Get CRISP transmission profile
   
  fpi = red_get_fpi_par(line = prefs)
  dw = xs[1] - xs[0]
  np = long((max(xs) - min(xs)) / dw) - 2
  if(np/2*2 eq np) then np -= 1L
  tw = (dindgen(np) - np/2) * dw
  tr = red_get_fpi_trans(fpi, tw + fpi.w0, ecl = 0.0, ech = 0.0, erh = -0.01)
                                ; 
  dum = max(tr, p)
  cc = poly_fit(tw[p-1:p+1] * 100.d0, tr[p-1: p+1], 2, /double)
  off = -0.005d0 * cc[1] / cc[2]
  tr = red_get_fpi_trans(fpi, tw + fpi.w0 + off, ecl = 0.0, ech = 0.0, erh = -0.01)
  
  ;; Convolve Solar Atlas with the FPI transmission profile
   
  ys = red_convl(ys, tr, /usefft)
  
  ;; Pack variables for mpfit
                                
  yl1 = red_intepf(fit.xl, fit.yl, wav[w0:w1])
  mm = {xl:fit.xl, yl:fit.yl, wav:wav[w0:w1], yl1:yl1}
  functargs = {xs:xs, ys:ys, dpr:dpr, mm:mm}
  
  ;; Init guess model
   
  pp = dblarr(8)
  pp[0] = fac                   ; Scale factor
  pp[1] = -0.5d0                ; Pref. shift
  pp[2] = 6.0d0                 ; Pref. FWHM
  pp[3] = 2.d0                  ; Pref. ncav
  pp[4] = -0.001d0              ; Line shift (satlas-obs)
  pp[5] = 0.001d0
  pp[6] = 0.001d0
  pp[7] = 1.0d0

  fitpars = replicate({mpside:2, limited:[0,0], limits:[0.0d, 0.0d], fixed:0}, 8)
   
  fitpars[2].LIMITED = [1,1]
  fitpars[2].LIMITS = [1.d0, 11.d0]
  pp[3] = 2.0d0
  fitpars[3].LIMITS = [1.4d0, 2.4d0]
  fitpars[3].LIMITED = [1,1]
  fitpars[7].FIXED =1B
   
  If(keyword_set(fixcav)) then begin
     pp[3] = fixcav
     fitpars[3].fixed = 1
  endif
  If(keyword_set(noasy)) then begin
     fitpars[5].fixed = 1
     fitpars[5].fixed = 1
     pp[5] = 0.
     pp[6] = 0.
  endif
  if(keyword_set(shift)) then begin
     pp[4] = shift
  endif

  ;; call mpfit
   
  pp = mpfit('red_fit_prefilter', pp, functargs = functargs, parinfo = fitpars, /quiet)
  dum = red_fit_prefilter(pp, xs = xs, ys = ys, dpr = dpr, mm = mm, pref = pref)
   
  print, inam + 'p[0] -> ', pp[0], ' (scale factor)'
  print, inam + 'p[1] -> ', pp[1], ' (prefilter shift)'
  print, inam + 'p[2] -> ', pp[2], ' (prefilter FWHM)'
  print, inam + 'p[3] -> ', pp[3], ' (prefilter number of cavities)'
  print, inam + 'p[4] -> ', pp[4], ' (solar atlas shift)'
  print, inam + 'p[5] -> ', pp[5], ' (asymmetry term 1)'
  print, inam + 'p[6] -> ', pp[6], ' (asymmetry term 2)'
  print, inam + 'p[7] -> ', pp[7], ' (Wavelength stretch)'

   
  odir = self.out_dir + '/prefilter_fits/'
  file_mkdir, odir
   
  ofile = cam + '.'+prefs+'.prefilter.f0'
  print, inam + 'saving prefilter to file -> ' + odir + ofile
  fzwrite, float(pref), odir + ofile, ' '
  ofile = cam + '.'+prefs+'.prefilter_wav.f0'
  fzwrite, float(mm.wav), odir + ofile, ' '
  
  ofile = cam + '.'+prefs+'.prefilter_pars.f0'
  print, inam + 'saving fit-results to file -> ' + odir + ofile
  fzwrite, pp, odir + ofile, ' '
  
  return
end
