function red_convl,a,b,plot=plot,usefft=usefft
;
; FFT based convolver
; A and B must be in the same grid
; A is padded to n_elements(b)/2+2 points
;
; Jaime de la Cruz Rodriguez (ISP-KVA 2010)
;
; a: function
; b: kernel of the convolution
;
  n0=size(a,/dim)
  n1=size(b,/dim)

;
; Pads a with half of the size of b.
; First/last value are repeated.
; Can use fft based convolution or just convol() from idl
;
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
;
; Plot original and result
;
  if keyword_set(plot) then begin
     plot,a,line=1,xtitle='X [grid point]',ytitle='Y [values]'
     oplot,res
     legend,['Convolved','Original'],line=[0,1],/right,/bottom
  endif
;
  return,res
end
