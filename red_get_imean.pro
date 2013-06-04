function red_get_imean, wav, dat, pp, npar, iter, xl = xl, rebin = rebin, densegrid = densegrid, thres = thres, myg = myg
  inam = 'red_get_imean : '
                                ;
  If(n_elements(rebin) eq 0) then rebin = 10L
  print, inam + 'Rebin factor -> ' + red_stri(rebin)
                                ;
  dim = size(dat, /dimension)
  rdim = long(dim / rebin)
                                ;
                                ; Build huge wavelength array (rebin by a factor)
                                ;
  imean = fltarr(dim[0])
                                ; iwav  = wav - median(pp[1,*,*])
                                ; if(keyword_set(myg) AND iter GT 0) then  iwav  = myg - median(pp[1,*,*]) else iwav  = wav - median(pp[1,*,*])
  if(keyword_set(myg) AND iter GE 0) then  iwav  = myg - median(pp[1,*,*]) else iwav  = wav - median(pp[1,*,*])

                                ;
  wl = fltarr(dim[0], dim[1], dim[2])
                                ;
  print, inam + 'ordering and correcting data to fit mean spectrum ... ', FORMAT = '(A,$)'
  for ii = 0L, dim[0] - 1 do begin
     wl[ii,*,*] = wav[ii] - pp[1,*,*]
  endfor

                                ;
  wl = reform(temporary(wl), dim[0]*dim[1]*dim[2])
                                ;
                                ; Same with the intensity
                                ;
  fl = dat                      ;fltarr(dim[0], dim[1], dim[2])
  for ii = 0L, dim[0] - 1 do begin
     mask = reform(dat[ii,*,*]) ge 1.e-3
     fl[ii,*,*] /= reform(pp[0,*,*]) * red_get_linearcomp(wav[ii], pp, npar)
     imean[ii] = median(reform(fl[ii,*,*])) 
  endfor
  
  if(keyword_set(myg) AND iter GE 0) then  imean = intepf(wav-median(pp[1,*,*]),imean, iwav)


                                ;
  fl = reform(temporary(fl), dim[0]*dim[1]*dim[2])
  count = dim[0]*dim[1]*dim[2]
                                ;
                                ; Remove pixels that are set to zero
                                ;
  idx = where(dat gt 0.01, count)
  if(count ge 1L) then begin
     wl = (temporary(wl))[idx]
     fl = (temporary(fl))[temporary(idx)]
  endif else begin
     print, inam+'Error -> all pixels are zero!'
     stop
  endelse
                                ;
                                ; sort elements by increasing wav
                                ;
  idx = sort(wl)
  wl = (temporary(wl))[idx]
  fl = (temporary(fl))[temporary(idx)]
  print, 'done'
                                ;
                                ; Rebin
                                ;
  nrebin = long(count / rebin)
  wl = rebin((temporary(wl))[0:nrebin*rebin-1], nrebin)
  fl = rebin((temporary(fl))[0:nrebin*rebin-1], nrebin)
                                ;
                                ; Use denser grid to fit spline?
                                ;
  bla = red_minmax_filter(wl)
  iwav[0] = min(wl[bla])
  iwav[n_elements(iwav)-1] = max(wl[bla])
                                ;
  if(keyword_set(densegrid) AND (iter ge 1)) then begin
     iwav1 = red_densegrid(iwav, thres = thres)
     imean = intepf(iwav, imean , iwav1)
     iwav = iwav1
     print, inam + 'Using denser grid of points to fit mean spectrum: '+red_stri(dim[0]) +$
            ' -> '+ red_stri(n_elements(iwav1))
  endif
                                ;

                                ;
                                ; fit to hermitian spline
                                ;
  loadct, 3, /silent
  mmmi = min(fl) > 0
  mmma = max(fl) < 3
  plot, wl, fl, psym = 3, xtitle = 'Wavelength', ytitle = 'Normalized intensity', $
        ystyle = 1, xstyle = 3,yrange=[mmmi,mmma]
                                ;
  oplot, iwav, imean, psym = 1 , color = 175
                                ;
  print, inam + 'fitting data to Hermitian spline with '+red_stri(n_elements(imean))+$
         ' node(s) points (this might take a while) ... ', FORMAT = '(A,$)'
  functargs = {wl:temporary(wl), fl:temporary(fl), iwav:iwav}
  if(iter gt 0) then  yl = mpfit('red_fit_hspline', imean, functargs = functargs, /quiet) else yl = imean
  xl = iwav
  print, 'done'
                                ;
                                ; plot result (use coarser grid, just for displaying)
                                ;
  loadct, 1, /silent
  np = 201L
  pr = (max(iwav) - min(iwav))
  pwl = findgen(np) / (np - 1.0) * pr + min(iwav)
  functargs = 0B 
                                ;
  oplot, pwl, intepf(xl, yl, pwl, /linear), color = 180
  loadct, 0, /silent
  wait, 0.2                     ; otherwise IDL does not update the plot (?)
                                ;
  return, yl
end
