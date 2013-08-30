; docformat = 'rst'

;+
; Fit surfaces to xtilts and ytilts. 
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
; :Params:
; 
;    simx
; 
; 
; 
;    simy
; 
; 
; 
;    sx
; 
; 
; 
;    sy
; 
; 
; 
;
;
; :Keywords:
; 
;    xtilts : in, optional
; 
; 
; 
;    ytilts : in, optional
; 
; 
; 
;    foc : in, optional
; 
; 
; 
;    dxoffs : out, optional
; 
; 
; 
;    dyoffs : out, optional
; 
; 
; 
;    dfoc : out, optional
; 
; 
; 
;    ddiv : out, optional
; 
; 
; 
; :History:
; 
; 
; 
; 
;-
pro red_pinh_make_fits, simx, simy, sx, sy $
                        , XTILTS = xtilts $ 
                        , YTILTS = ytilts $ 
                        , FOC = foc $       
                        , DXOFFS = dxoffs $ 
                        , DYOFFS = dyoffs $ 
                        , DFOC = dfoc $     
                        , DDIV = ddiv       

  ;; Use outlier resistant quartic fits? Or see how Pit did it in the
  ;; python code? Then add those surfaces to xoffs and yoffs.

  use_robust = 1                ; Use robust means and plane fitting

  Nch = (size(xtilts, /dim))[0] 
  Npinhx = (size(simx, /dim))[0]
  Npinhy = (size(simy, /dim))[0]

  if n_elements(diversity) eq 0 then diversity = replicate(0.0, Nch)

  dxoffs = fltarr(sx, sy, Nch)
  dyoffs = fltarr(sx, sy, Nch)

  ;; Coordinates for the fitting
  simxx=transpose(reform(rebin(simx,Npinhy*Npinhx,/samp),Npinhy,Npinhx))
  simyy=reform(rebin(simy,Npinhy*Npinhx,/samp),Npinhx,Npinhy)

  ;; Coordinates for the offset files
  sxx=transpose(reform(rebin(findgen(sx),sx*sy,/samp),sy, sx))
  syy=reform(rebin(findgen(sy),sx*sy,/samp),sx, sy)

  dxtilts = fltarr(Nch, Npinhy, Npinhx)
  dytilts = fltarr(Nch, Npinhy, Npinhx)
  dfoc = fltarr(Nch, Npinhy, Npinhx)
  ddiv = fltarr(Nch)

  for ich = 1, Nch-1 do begin

     ;; X tilts
     dxtilts[ich, *, *] = xtilts[ich, *, *]-xtilts[0, *, *]
     if use_robust then begin
        dxtiltsc = robust_planefit(double(simxx), double(simyy), dxtilts[ich, *, *], dxtiltsfit)
     endif else begin
        dxtiltsc = planefit(simxx, simyy, dxtilts[ich, *, *], 0, dxtiltsfit)
     endelse
     dxoffs[*, *, ich] = dxtiltsc[0]+dxtiltsc[1]*sxx+dxtiltsc[2]*syy
     
     ;; Y tilts
     dytilts[ich, *, *] = ytilts[ich, *, *]-ytilts[0, *, *]
     if use_robust then begin
        dytiltsc = robust_planefit(double(simxx), double(simyy), dytilts[ich, *, *], dytiltsfit)
     endif else begin
        dytiltsc = planefit(simxx, simyy, dytilts[ich, *, *], 0, dytiltsfit)
     endelse 
     dyoffs[*, *, ich] = dytiltsc[0]+dytiltsc[1]*sxx+dytiltsc[2]*syy

;     mnn = min([dxtiltsfit, dytiltsfit])*3
;     mxx = max([dxtiltsfit, dytiltsfit])*3
;
;     tvscl, bytscl([rebin(reform(dxtilts[ich, *, *]), Npinhx*10, Npinhy*10, /samp) $
;                    , rebin(reform(dxtiltsfit, Npinhx, Npinhy), Npinhx*10, Npinhy*10, /samp) $
;                   ], mnn, mxx), 0
;
;     tvscl, bytscl([rebin(reform(dytilts[ich, *, *]), Npinhx*10, Npinhy*10, /samp) $
;                    , rebin(reform(dytiltsfit, Npinhx, Npinhy), Npinhx*10, Npinhy*10, /samp) $
;                   ], mnn, mxx), 1


     if n_elements(foc) ne 0 then begin
        
        ;; Focus diversity
        ;; The plane fit is just for display, the update is just
        ;; the mean value.
        dfoc[ich, *, *] = foc[ich, *, *]-foc[0, *, *]
        if use_robust then begin
           dfocc = robust_planefit(simxx, simyy, dfoc, dfocfit)
           ddiv[ich] = biweight_mean(dfoc[ich, *, *]) 
        endif else begin
           dfocc = planefit(simxx, simyy, dfoc, 0, dfocfit)
           ddiv[ich] = mean(dfoc[ich, *, *]) 
        endelse

;        mnn = min(dfocfit)*3
;        mxx = max(dfocfit)*3
;
;        tvscl, bytscl([rebin(reform(dfoc), Npinhx*10, Npinhy*10, /samp) $
;                       , rebin(reform(dfocfit, Npinhx, Npinhy), Npinhx*10, Npinhy*10, /samp) $
;                      ], mnn, mxx), 2

     endif                      ; if foc     
  endfor                        ; for ich
end
