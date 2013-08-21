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
pro red_cfitgain, par, wav, dat, xl, yl, ratio, nthreads = nthreads
  if(n_elements(nthreads) eq 0) then nthreads = 4L 
  npar = (size(par, /dim))[0]
  dim = size(dat, /dim)
  npix = dim[1] * dim[2]
  nwav = dim[0]
  nl = n_elements(xl)
                                ;
  dir=getenv('CREDUC')
                                ;
  b = call_external(dir+'/creduc.so', 'cfitgain', long(nwav), long(nl), long(npar), long(npix), float(xl), float(yl), float(wav), dat, par, ratio, long(nthreads))
                                ;
  return
end
