; docformat = 'rst'

;+
; Plot convergence of pinhole calibrations.
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, Institute for Solar Physics, 2012
; 
; 
; :Keywords:
; 
;    it : in, optional
; 
;       # of iterations so far
; 
;    mtol : in, optional
; 
;       Metric tolerance.
; 
;    dtol : in, optional
; 
;       Diversity tolerance.
; 
;    ttol : in, optional
; 
;       Tilt tolerance. 
; 
;    xconvergence : out, optional
; 
;       X tilt convergence
;
;    yconvergence : out, optional
; 
;       Y tilt convergence
;
; 
;    dconvergence : out, optional
; 
;       Diversity convergence
;
; 
;    mconvergence : out, optional
; 
;       Metric convergence
;
;    dvalues : out, optional
; 
;        Diversity convergence 
;    
; 
; 
; 
; :History:
; 
; 
; 
; 
;-
pro red_pinh_plot_convergence, IT = it $ ; # of iterations so far
                               , MTOL = mtol $
                               , DTOL = dtol $
                               , TTOL = ttol $
                               , XCONVERGENCE = xconvergence $ ; X tilt convergence (out)
                               , YCONVERGENCE = yconvergence $ ; Y tilt convergence (out)
                               , DCONVERGENCE = dconvergence $ ; Diversity convergence (out)
                               , DVALUES = dvalues $           ; Diversity convergence (out)
                               , MCONVERGENCE = mconvergence   ; Metric convergence (out)

  if n_elements(it) eq 0 then plotall = 1 else plotall = 0

  if n_elements(xconvergence) ne 0 then begin
     window, 10, xs = 400, ys = 300
     Nch =  (size(xconvergence, /dim))[0]
     if plotall then it = (size(xconvergence, /dim))[1]-1
     plot, xconvergence[0, 0:it], title = 'xconvergence' $
           , yrange = [min(xconvergence) <0, max(xconvergence) >0], /nodata
     for ich = 1, Nch-1 do oplot, xconvergence[ich, 0:it]
     if n_elements(ttol) ne 0 then begin
        oplot, [0, 1]*it*2,  ttol*[1, 1], linestyle = 1
        oplot, [0, 1]*it*2, -ttol*[1, 1], linestyle = 1
     endif
     oplot, [0, 1]*it*2, [1, 1]*0, linestyle = 2
  endif

  if n_elements(yconvergence) ne 0 then begin
     window, 11, xs = 400, ys = 300
     Nch =  (size(yconvergence, /dim))[0]
     if plotall then it = (size(yconvergence, /dim))[1]-1
     plot, yconvergence[0, 0:it], title = 'yconvergence' $
           , yrange = [min(yconvergence) <0, max(yconvergence) >0], /nodata
     for ich = 1, Nch-1 do oplot, yconvergence[ich, 0:it]
     if n_elements(ttol) ne 0 then begin
        oplot, [0, 1]*it*2,  ttol*[1, 1], linestyle = 1
        oplot, [0, 1]*it*2, -ttol*[1, 1], linestyle = 1
     endif
     oplot, [0, 1]*it*2, [1, 1]*0, linestyle = 2
  endif


  if n_elements(dconvergence) ne 0 then begin
     window, 12, xs = 400, ys = 300
     Nch =  (size(dconvergence, /dim))[0]
     if plotall then it = (size(dconvergence, /dim))[1]-1
     plot, dconvergence[0, 0:it], title = 'diversity updates' $
           , yrange = [min(dconvergence) <0, max(dconvergence) >0], /nodata
     for ich = 1, Nch-1 do oplot, dconvergence[ich, 0:it]
     if n_elements(dtol) ne 0 then begin
        oplot, [0, 1]*it*2,  dtol*[1, 1], linestyle = 1
        oplot, [0, 1]*it*2, -dtol*[1, 1], linestyle = 1
     endif
     oplot, [0, 1]*it*2, [1, 1]*0, linestyle = 2
  endif
     
  if n_elements(dvalues) ne 0 then begin
     window, 13, xs = 400, ys = 300
     Nch =  (size(dvalues, /dim))[0]
     if plotall then it = (size(dvalues, /dim))[1]-1
     plot, dvalues[0, 0:it], title = 'diversity values' $
           , yrange = [min(dvalues) <0, max(dvalues) >0], /nodata
     for ich = 1, Nch-1 do oplot, dvalues[ich, 0:it]
  endif

  if n_elements(mconvergence) ne 0 then begin
     window, 14, xs = 400, ys = 300
     if plotall then it = (size(mconvergence, /dim))[2]-1
     maxmconvergence = max(max(mconvergence, dim = 1), dim = 1)
     minmconvergence = min(min(mconvergence, dim = 1), dim = 1)
     dims = size(mconvergence, /dim)
     medmconvergence = median(reform(mconvergence, dims[0]*dims[1], dims[2]), dim = 1)
     mn = min(minmconvergence[where(minmconvergence ne 0)])
     mx = max(maxmconvergence)
     plot, maxmconvergence[0:it], title = 'MOMFBD metric', /ylog, yrange = [mn, mx]
     oplot, medmconvergence[0:it]
     oplot, minmconvergence[0:it]
  endif

end
