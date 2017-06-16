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
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
; 
;-
pro pol::loadimages 
   
  nt = n_elements(self.tfiles)
   
  ;; Clean-up the heap variable (by default a null pointer)
                                ;
  ptr_free, self.timg
  ptr_free, self.rimg
   
  ;; Make sure that the heap variable is a valid pointer 
   
  self.timg = ptrarr(nt, /allocate_heap)
  self.rimg = ptrarr(nt, /allocate_heap)
   
  ;; Load momfbd files
  if(self.ftype  eq 'momfbd') then begin
     for ii = 0L, nt - 1 do *self.timg[ii] = momfbd_read(self.tfiles[ii])
     for ii = 0L, nt - 1 do *self.rimg[ii] = momfbd_read(self.rfiles[ii])
  endif else begin
     for ii = 0L, nt - 1 do *self.timg[ii] = f0(self.tfiles[ii])
     for ii = 0L, nt - 1 do *self.rimg[ii] = f0(self.rfiles[ii])
  endelse
     
   
  return
end
