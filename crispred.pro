function crispred, filename
                                ;
                                ; device, decompose=0
  done = {done, sumdark:0B, sumflat:0B, cleandata:0B, sumpolcal:0B, polcal:0B, makegains:0B}
                                ;
  struct = {red, $
            dark_dir:'',$
            flat_dir:'',$
            data_dir:'',$
            data_list:strarr(100),$
            ndir:0B, $
            pinh_dir:'',$
            prefilter_dir:'',$
            polcal_dir:'',$
            camt:'', $
            camr:'', $
            camwb:'', $
            out_dir:'',$
            filename:'',$
            dopolcal:0B, $
            dodata:0B, $
            doflat:0B, $
            dodark:0B, $
            dopinh:0B, $
            docamt:0B, $
            docamr:0B, $
            docamwb:0B,$
            camttag:'',$
            camrtag:'',$
            camwbtag:'',$
            descatter_dir:'',$
            root_dir:'',$
            dodescatter:0B,$
            telescope_d:'',$
            image_scale:'',$
            pixel_size:'',$
            camsz:'',$
            done:done}
                                ;
  tmp = obj_new('red')
                                ;
                                ; check input file
                                ;
  if((n_params() eq 0) AND file_test('config.txt')) then filename = 'config.txt'
  if(~file_test(filename)) then begin
     print, 'reduc : ERROR, cannot find '+filename
     return, 0
  endif
                                ;
  tmp -> initialize, filename
                                ;
  return, tmp
end
