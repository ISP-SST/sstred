; docformat = 'rst'

;+
; Convolve a function with a kernel (1D), optionally using FFT.
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :author:
; 
;    Jaime de la Cruz Rodriguez (ISP-KVA 2010)
; 
; 
; :returns:
; 
; 
; :Params:
; 
;   A : in, 
;   
;      The "function" A. A is padded to n_elements(B)/2+2 points.
;   
;   B : 
;   
;      The "kernel". Must be on the same grid as A.
;   
; 
; :Keywords:
; 
;   plot : in, boolean
;   
;     Set this to make a plot.
;   
;   usefft : in, boolean
;   
;     Set this to use FFT, otherwise use convol() from IDL.
;   
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
; 
;-
function red_convl, A, B, plot=plot, usefft=usefft

  n0=size(a,/dim)
  n1=size(b,/dim)

  ;; Pads A with half of the size of B.
  ;; First/last value are repeated.

  if not keyword_set(usefft) then begin
     res=convol(a,b/total(b),/edge_truncate)
  endif else begin
     npad0 = n1/2+2         
     npad1 = 2L*npad0+n0-n1 
     aa = [replicate(a[0],npad0), a, replicate(a[n0-1],npad0)]
     bb = shift([b,replicate(0.d0,npad1)], -n1/2)
     bb/=total(bb)
     res=double((fft(fft(aa,1)*fft(bb,1),-1))[npad0:npad0+n0-1])
  endelse

  ;; Optionally plot original and result

  if keyword_set(plot) then begin
     plot,a,line=1,xtitle='X [grid point]',ytitle='Y [values]'
     oplot,res
     legend,['Convolved','Original'],line=[0,1],/right,/bottom
  endif

  return, res

end
