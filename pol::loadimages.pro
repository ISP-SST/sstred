pro pol::loadimages 
                                ;
  nt = n_elements(self.tfiles)
                                ;
                                ; Clean-up the heap variable (by default a null pointer)
                                ;
  ptr_free, self.timg
  ptr_free, self.rimg
                                ;
                                ; Make sure that the heap variable is a valid pointer 
                                ;
  self.timg = ptrarr(nt, /allocate_heap)
  self.rimg = ptrarr(nt, /allocate_heap)
                                ;
                                ; Load momfbd files
                                ;
  for ii = 0L, nt - 1 do *self.timg[ii] = momfbd_read(self.tfiles[ii])
  for ii = 0L, nt - 1 do *self.rimg[ii] = momfbd_read(self.rfiles[ii])
                                ;
  return
end
