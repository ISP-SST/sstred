; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :author:
; 
; 
; 
; 
; :returns:
; 
; 
; :Params:
; 
;    par : 
;   
;   
;   
;    wav : 
;   
;   
;   
;    dat : 
;   
;   
;   
;    xl : 
;   
;   
;   
;    yl : 
;   
;   
;   
;    ratio : 
;   
;   
;   
; 
; :Keywords:
; 
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
; 
;-
pro red_cfitgain2, par, wav, dat, xl, yl, ratio, line,nthreads = nthreads
  if(n_elements(nthreads) eq 0) then nthreads = 4L 
  npar = (size(par, /dim))[0]
  dim = size(dat, /dim)
  npix = dim[1] * dim[2]
  nwav = dim[0]
  nl = n_elements(xl)
                                ;
  libfile = red_libfile('creduc.so')
                                ;

  b = call_external(libfile, 'cfitgain2', long(nwav), long(nl), long(npar), $
                    long(npix), float(xl), float(yl), float(wav), dat, par, $
                    ratio, long(nthreads), long(line))
                                ;
  return
end
