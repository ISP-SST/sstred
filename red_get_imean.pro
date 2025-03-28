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
;    wav : 
;   
;   
;   
;    dat : 
;   
;   
;   
;    pp : 
;   
;   
;   
;    npar : 
;   
;   
;   
;    iter : 
;   
;   
;   
; 
; :Keywords:
; 
;    xl  : 
;   
;   
;   
;    rebin  : 
;   
;   
;   
;    densegrid  : 
;   
;   
;   
;    thres  : 
;   
;   
;   
;    myg  : 
;   
;   
;   
; 
; 
; :History:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2013-07-11 : MGL. Use red_intepf, not intepf.
;
;   2016-10-04 : MGL. Use cgplot. Let the subprogram find out its own
;                name. 
;
;   2016-10-26 : MGL. Add progress bar, don't redraw plots so often.  
;
;   2017-04-20 : MGL. Plot spline nodes on the model curve.
;
;   2023-11-02 : JdlCR. Modifications for multiple LC fitting with
;                a single cavity map value for the new demodulation
;                and flat-fielding scheme.
; 
;-
function red_get_imean, wav, dat, pp, npar, iter $
                        , xl = xl $
                        , rebin = rebin $
                        , densegrid = densegrid $
                        , thres = thres $
                        , myg = myg $
                        , reflect = reflect $
                        , bezier=bezier $
                        , title = title
  

  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)            

  if n_elements(title) eq 0 then title = ''

  ddims = size(dat, /dim)
  Nlc = ddims[1]
  Nwav = ddims[0]

  npar_one = (npar-1)/Nlc
  npol = npar_one - 1
  
  
  ;; Check for NaNs

  dum = (size(pp, /dim))[0]
  for ii=0,dum-1 do begin
    tmp = reform(pp[ii,*,*])
    pos = where(~finite(tmp), nancount, complement=pos1)
    if(nancount gt 0) then begin
      tmp[pos] = median(tmp[pos1])
      pp[ii,*,*] = temporary(tmp)
    endif
  endfor

  If(n_elements(rebin) eq 0) then rebin = 10L
  print, inam + ' : Rebin factor -> ' + red_stri(rebin)
  
  dim = size(dat, /dimension)
  rdim = long(dim / rebin)
  
  ;; Build huge wavelength array (rebin by a factor)

  imean = reform(fltarr(dim[0],dim[1]), [dim[0], dim[1]])
  ;; iwav  = wav - median(pp[1,*,*])
  ;; if(keyword_set(myg) AND iter GT 0) then  iwav  = myg - median(pp[1,*,*]) else iwav  = wav - median(pp[1,*,*])
  if keyword_set(myg) and iter ge 0 then iwav = myg - median(pp[1,*,*]) else iwav = wav - median(pp[1,*,*])

  wl = fltarr(dim[0], dim[1], dim[2], dim[3])

  print, inam + ' : Ordering and correcting data to fit mean spectrum ... ', FORMAT = '(A,$)'
  print

  for ii = 0L, dim[0] - 1 do begin
    red_progressbar, ii, dim[0], 'Ordering and correcting wavelength array to fit mean spectrum.'
    for ss =0, dim[1]-1 do begin
      wl[ii,ss,*,*] = wav[ii] - pp[1,*,*]
    endfor
  endfor

  wl = reform(temporary(wl), dim[0]*dim[1]*dim[2]*dim[3])

  ;; Same with the intensity

  fl = dat                      ;fltarr(dim[0], dim[1], dim[2])


  tmp_par = fltarr(npar_one+1, dim[2], dim[3])
  tmp_par[1,*,*] = pp[1,*,*]

  ntot = dim[1]*dim[0]
  for ss = 0L, dim[1]-1 do begin
    
    if(ss ne 0) then begin
      off = 1 + ss*npar_one ;;npar_one + 1 + (ss-1)*npar_one

      tmp_par = fltarr(npar_one+1, dim[2], dim[3])
      tmp_par[0,*,*] = pp[off,*,*]
      if(npol ge 1) then tmp_par[2:*,*,*] = pp[off+1:off+npol,*,*]
    endif else begin
      tmp_par = pp[0:npar_one,*,*]
    endelse
    
    for ii = 0L, dim[0] - 1 do begin
      red_progressbar, ii+ss*dim[0], ntot, 'Correcting the intensity array'
      ;;mask = reform(dat[ii,ss,*,*]) ge 1.e-3
      lcom = red_get_linearcomp(wav[ii], tmp_par, npar_one+1,reflect=reflect)
      ;; Remove mean shape from the linear component so we don't
      ;; introduce distortion in the average profile
      lcom /= median(lcom) 
      fl[ii,ss,*,*] /= reform(tmp_par[0,*,*]) * lcom
      imean[ii,ss] = median(reform(fl[ii,ss,*,*])) 
    endfor
  endfor
  
  imean = total(imean, 2)/Nlc
  
  if(keyword_set(myg) AND iter GE 0) then  begin
    if(~keyword_set(bezier)) then begin
      imean = red_intepf(wav-median(pp[1,*,*]),imean, iwav)
    endif else begin
      imean = red_bezier3(wav-median(pp[1,*,*]),imean, iwav)
    endelse
  endif

  fl = reform(temporary(fl), dim[0]*dim[1]*dim[2]*dim[3])
  count = dim[0]*dim[1]*dim[2]*dim[3]

  ;; Remove pixels that are set to zero

  idx = where(dat gt 0.01, count)
  if(count ge 1L) then begin
    wl = (temporary(wl))[idx]
    fl = (temporary(fl))[temporary(idx)]
  endif else begin
    print, inam+' : Error -> all pixels are zero!'
    stop
  endelse

  ;; Sort elements by increasing wav

  idx = sort(wl)
  wl = (temporary(wl))[idx]
  fl = (temporary(fl))[temporary(idx)]
;  print, 'done'

  ;; Rebin
  
  nrebin = long64(count / rebin)
  wl = rebin((temporary(wl))[0:nrebin*rebin-1], nrebin)
  fl = rebin((temporary(fl))[0:nrebin*rebin-1], nrebin)

  ;; Use denser grid to fit spline?

  bla = red_minmax_filter(wl)
  iwav[0] = min(wl[bla])
  iwav[n_elements(iwav)-1] = max(wl[bla])

  if keyword_set(densegrid) AND (iter ge 1) then begin
    iwav1 = red_densegrid(iwav, thres = thres)
    imean = red_intepf(iwav, imean , iwav1)
    iwav = iwav1
    print, inam + ' : Using denser grid of points to fit mean spectrum: ' $
           + red_stri(dim[0]) + ' -> '+ red_stri(n_elements(iwav1))
  endif

  ;; Fit to hermitian spline

  mmmi = min(fl) > 0
  mmma = max(fl) < 3

  cgcontrol, /delete, /all
  cgwindow, 'cgplot', /load, wl, fl, psym = 3, title = title $
            , xtitle = 'Wavelength / 1 $\Angstrom$', xstyle = 3 $
            , ytitle = 'Normalized intensity', ystyle = 1 $
            , yrange = [mmmi-abs(mmma-mmmi)*0.1,mmma+abs(mmma-mmmi)*0.1]
  
;  cgwindow, 'cgplot', /load, /over, iwav, imean, psym = 9, color = 'yellow'
;  cgwindow, 'cgplot', /load,  /over, iwav, imean, psym = 1, color = 'red'
  
  print, inam + ' : Fitting data to Hermitian spline with '+red_stri(n_elements(imean)) $
         + ' node points (this might take a while) ... ', FORMAT = '(A,$)'

  if(keyword_set(bezier)) then dobezier=1 else dobezier = 0

  functargs = {wl:temporary(wl), fl:temporary(fl), iwav:iwav, bezier:dobezier}
  if(iter gt 0) then yl = mpfit('red_fit_hspline', imean, functargs = functargs, /quiet) else yl = imean
  xl = iwav
  print, 'done'
  
  ;; Plot result (use coarser grid, just for displaying)
  
  np = 201L
  pr = (max(iwav) - min(iwav))
  pwl = findgen(np) / (np - 1.0) * pr + min(iwav)
  functargs = 0B 
  
                                ;
  if(keyword_set(bezier)) then begin
    cgwindow, 'cgplot', /load, /over, pwl, red_bezier3(xl, yl, pwl, /linear), color = 'blue'
    yyy = red_bezier3(xl, yl, iwav, /linear)
  endif else  begin
    cgwindow, 'cgplot', /load, /over, pwl, red_intepf(xl, yl, pwl, /linear), color = 'blue'
    yyy = red_intepf(xl, yl, iwav, /linear)
  endelse
  
  cgwindow, 'cgplot', /load, /over, iwav, yyy, psym = 9, color = 'yellow'
  cgwindow, 'cgplot', /add,  /over, iwav, yyy, psym = 1, color = 'red'
  
;  wait, 0.2                     ; otherwise IDL does not update the plot (?)

  return, yl

end
