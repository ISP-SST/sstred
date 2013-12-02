function fftconvol, yl, tr
  n = n_elements(yl)
  n1 = n_elements(tr)
  nn1 = n1
  npad = n + n1 - 1 
;  if(n1/2*2 ne n1) then npad -= 1

  ylp = dblarr(npad)
  ylp[0:n-1] = yl
  ylp[n:n+n1/2-1]=yl[n-1]
  ylp[n+n1/2:*]=yl[0]

  trr = shift([tr,replicate(0, npad-n1)], -n1/2)
  
  return,double(((fft(fft(ylp, 1) * fft(trr,1), -1)))[0:n-1])
end
