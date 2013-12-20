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
;    npar  : 
;   
;   
;   
;    niter  : 
;   
;   
;   
;    rebin  : 
;   
;   
;   
;    xl  : 
;   
;   
;   
;    yl  : 
;   
;   
;   
;    densegrid  : 
;   
;   
;   
;    res  : 
;   
;   
;   
;    thres  : 
;   
;   
;   
;    initcmap  : 
;   
;   
;   
;    x0  : 
;   
;   
;   
;    x1  : 
;   
;   
;   
;    state  : 
;   
;   
;   
;    nosave  : 
;   
;   
;   
;    myg  : 
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
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2013-08-27 : MGL. Added support for logging. Let the subprogram
;                find out its own name.
; 
;   2013-12-20 : PS join with _ng version:  reflectivity as keyword,
;                   npar is number of *additional* terms (default: 1 = linear)
; 
;-
pro red::fitgains, npar = npar, niter = niter, rebin = rebin, xl = xl, yl = yl, $
                   densegrid = densegrid, res = res, thres = thres, initcmap = initcmap, $
                   fit_reflectivity = fit_reflectivity, $
                   x0 = x0, x1 = x1, state = state, nosave = nosave, myg = myg, $
                   w0 = w0, w1 = w1, nthreads = nthreads

  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  ;; Logging
  help, /obj, self, output = selfinfo 
  red_writelog, selfinfo = selfinfo

    ;;; get_imean uses colortables.  Save old graphics setting
  device, get_decompose = odc
  device, decompose = 0

  if(~keyword_set(npar)) then npar = 1L
  npar_t = npar + (keyword_set(fit_reflectivity) ? 3 : 2)

  if(~keyword_set(niter)) then niter = 3L
  if(~keyword_set(rebin)) then rebin = 100L
  
  ;; Restore data
  files = self.out_dir + 'flats/spectral_flats/cam*.flats.sav'

  files = file_search(files, count = ct)
  prefs = strarr(ct)
  cams = strarr(ct)
  states = strarr(ct)

  print, inam + ' : found states:'
  for ii = 0, ct - 1 do begin
     tmp = strsplit(file_basename(files[ii]), '.',/extract)
     cams[ii] = tmp[0]
     prefs[ii] = tmp[1]
     states[ii] = strjoin(tmp[0:1],'.')
     print, ii, ' -> ',states[ii], FORMAT='(I3,A,A)'
  endfor

  if(~keyword_set(state)) then begin
     idx = 0
     if(ct gt 1) then read, idx, prompt =  inam + ' : select state ID to be processed: '
  endif else begin
     idx = where(states eq state, count)
     if count eq 0 then begin
        print, inam + ' : ERROR, external state -> '+state+' not found'
        return
     endif
  endelse
   
  ;; Load data
  cam = cams[idx]
  pref = prefs[idx]
  print, inam + ' : selected -> '+states[idx]
  restore, files[idx]
  if(keyword_set(w0)) then begin
     cub = (temporary(cub))[w0:*,*,*]
     wav = wav[w0:*]
     namelist = namelist[w0:*]
  endif else w0=0
  if(keyword_set(w1)) then begin
     cub = (temporary(cub))[0:w1-w0,*,*]
     wav = wav[0:w1-w0]
     namelist = namelist[0:w1-w0]
  endif

  wav = float(wav)
  if(n_elements(dat) eq 0) then begin 
     dat = temporary(cub)
  endif

  ;; Init output vars
  dim = size(dat, /dim)
  res = dblarr(npar_t, dim[1], dim[2])
  ratio = fltarr(dim[0], dim[1], dim[2])
  nwav = dim[0]
   
  res[0,*,*] = total(dat,1) / nwav
   ;;; is this still needed?
  IF keyword_set(fit_reflectivity) THEN IF npar_t GT 3 THEN res[3:*,*,*] = 1.e-3
  
  ;; Init cavity map?   
  if(keyword_set(initcmap)) then begin
     print, inam + ' : Initializing cavity-errors with parabola-fit'
     res[1,*,*] = red_initcmap(wav, dat, x0 = x0, x1 = x1)
  endif
  
  ;; Loop niter
  for it = 0L, niter - 1 do begin
   
     ;; Get mean spectrum using Hermitian Spline
     yl = red_get_imean(wav, dat, res, npar_t, it $
                        , xl = xl, rebin = rebin, densegrid = densegrid, thres = thres, $
                        myg = myg, reflec = fit_reflectivity)
      
     ;; Pixel-to-pixel fits using a C++ routine to speed-up things
     if(it eq 0) then begin
        res1 = res[0:1,*,*]
        red_cfitgain, res1, wav, dat, xl, yl, ratio, nthreads=nthreads
        res[0:1,*,*] = temporary(res1)
     ENDIF ELSE BEGIN
         IF keyword_set(fit_reflectivity) THEN $
           red_cfitgain2, res, wav, dat, xl, yl, ratio, pref, nthreads = nthreads $
         ELSE $
           red_cfitgain, res, wav, dat, xl, yl, ratio, nthreads = nthreads
     ENDELSE
  endfor
  yl = red_get_imean(wav, dat, res, npar_t, it, xl = xl, rebin = rebin, densegrid = densegrid, $
                       thres = thres, myg = myg, reflec = fit_reflectivity)

   ;;; reset the color setting
  device, decompose = odc
  
  ;; Create cavity-error-free flat (save in "ratio" variable)
  print, inam + ' : Recreating cavity-error-free flats ... ', FORMAT='(A,$)'
   
  for ii = 0L, nwav - 1 do ratio[ii,*,*] *= reform(res[0,*,*]) * $
     reform(red_get_linearcomp(wav[ii], res, npar_t, reflec = fit_reflectivity))
   
  print, 'done'
   
  ;; Save results
  if(~keyword_set(nosave)) then begin

     outdir = self.out_dir + '/flats/'
     for ii = 0L, nwav - 1 do begin
        namout =  outdir + namelist[ii]
        fzwrite, reform(ratio[ii,*,*]), namout, 'npar='+red_stri(npar_t)
        print, inam + ' : saving file -> '+namout
     endfor
   
     fit = {pars:res, yl:yl, xl:xl, oname:namelist}
     save, file = outdir+'spectral_flats/'+cam+'.'+pref+'.fit_results.sav', fit
     fit = 0B
  endif

  return
end
