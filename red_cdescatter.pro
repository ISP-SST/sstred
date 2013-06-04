function red_cdescatter, img, fgain, fpsf, nthreads = nthreads, verbose = verbose
  dim = size(img, /dimension)
  nx = dim[0]
  ny = dim[1]
                                ; Ensure input is a float array
  dimg = fltarr(nx, ny)
                                ;
  if(n_elements(verbose) eq 0) then verbose = 0L
  if(n_elements(nthreads) eq 0) then nthreads = 4L
                                ;
  dir=getenv('CREDUC')
                                ;
  b = call_external(dir+'/creduc.so', 'cdescatter', long(nx), $
                    long(ny), float(img), float(fgain), float(fpsf), $
                    dimg, long(nthreads), long(verbose))
                                ;
  return, dimg
end
