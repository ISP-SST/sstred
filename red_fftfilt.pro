function red_fftfilt, img, filter
  inam = 'red_fftfilt : '
                                ;
                                ; Square image!
                                ;
  dim = size(img,/dim)
                                ;if(dim[0] ne dim[1]) then begin
                                ;   print, inam + 'Error, square image needed -> returning'
                                ;   return, 0
                                ;endif
                                ;
  me = median(img)
                                ; win = red_taper(dim[0]*[1.,0.,1./24.])
  win = red_taper2(dim[0], dim[1],1./24.)
                                ;
;  show, (img-me)*win
  return, float(fft(fft((img - me) * win, 1) * filter, -1)) + me
end
