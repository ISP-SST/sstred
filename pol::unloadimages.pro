pro pol::unloadimages
                                ;
  nt = n_elements(self.tfiles)
                                ;
                                ; Free the memory that is pointed at!
                                ;
  ptr_free, self.timg
  ptr_free, self.rimg
                                ;
  return
end
