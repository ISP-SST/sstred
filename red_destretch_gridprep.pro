; docformat = 'rst'

;+
; Detrends and unsharp masks a given destretch grid.
;
; :author:
;
;    Tom Berger, LMSAL
;
; :returns:
;
;    Delta prep'd for application to image time series
;
;
; :params:
;
;     delta ; in, type="fltarr(nt,2,ngx,ngy)" 
;
;       Destretch grid, e.g. as returned by dsgridnest.pro
;
;     tstep : in
;
;       Smooth factor applied in the unsharp mask step
; 
;
; :history:
;
;   2013-07-24 : Renamed for inclusion in the crispred pipeline.
;
;-
FUNCTION red_destretch_gridprep, delta, tstep

  sz = SIZE(delta)
  IF sz[2] NE 2 THEN BEGIN
     MESSAGE, 'Destretch grid has funny dimensions - returning.'
     RETURN, -1
  END
  IF N_ELEMENTS(tstep) EQ 0 THEN BEGIN
     PRINT, 'tstep not specified - assuming default of 15.'
     tstep =  15
  END

  ;;Zero out any IEEE NaNs coming from GRIDMATCH:
  badi =  WHERE(FINITE(delta) eq 0, num)
  IF num GT 0 THEN delta[badi] = 0.0

  delx = REFORM(delta[*,0,*,*])
  dely = REFORM(delta[*,1,*,*])
  dsz = SIZE(delx)
  nt = dsz[1]
  xvec = INDGEN(nt)
  ngx = dsz[2]
  ngy = dsz[3]

  PRINT,'Detrending and unsharping the displacements...'
  t0 = SYSTIME(1)
  FOR j=0,ngy-1 DO for i=0,ngx-1 DO BEGIN 

     xq = TOTAL(delx[*,i,j],/cumu)  

     ;;Note that we apply a cumulative displacement grid to the time
     ;;series.
     
     xq1 = xq[1]                ;anchor point for first stretched image
     cf = POLY_FIT(xvec,xq,1,yfit)
     xq = xq - yfit
     xq = xq - SMOOTH(xq,tstep, /EDGE_TRUNCATE)
     xq = xq - xq[1] + xq1
     delx[0,i,j] = xq

     yq = TOTAL(dely[*,i,j],/cumu)
     yq1 = yq[1]
     cf = POLY_FIT(xvec,yq,1,yfit)
     yq = yq - yfit
     yq = yq - SMOOTH(yq,tstep, /EDGE_TRUNCATE)
     yq = yq - yq[1] + yq1
     dely[0,i,j] = yq

  END
  PRINT,'Detrend and unsharp time = ',SYSTIME(1)-t0,' seconds' 

  ;; Return variable:

  deltap = delta
  deltap[*,0,*,*] = delx
  deltap[*,1,*,*] = dely
  deltap[0, *, *, *] = 0.0      ;return reference image displacements to 0.

  RETURN, deltap
END
