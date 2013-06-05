;; Load class definitions
@pol::assign_states.pro
@pol::demodulate2.pro
@pol::demodulate.pro
@pol::fillclip.pro
@pol::getvar.pro
@pol::loadimages.pro
@pol::setvar.pro
@pol::state.pro
@pol::unloadimages.pro
@red::count2diskcenter.pro
@red::fitgains_ng.pro
@red::fitprefilter.pro
@red::getcamtags.pro
@red::getoffsets.pro
@red::get_scansquality.pro
@red::getstats.pro
@red::initialize.pro
@red::link_data.pro
@red::make_cmap_intdif2.pro
@red::make_cmap_intdif.pro
@red::make_cmaps.pro
@red::makegains.pro
@red::make_intdif_gains2.pro
@red::make_intdif_gains.pro
@red::make_pol_crispex.pro
@red::make_stokes_crispex2.pro
@red::make_stokes_crispex.pro
@red::make_time_cmap.pro
@red::make_time_series.pro
@red::make_unpol_crispex.pro
@red::polarim.pro
@red::polcalcube.pro
@red::polcal.pro
@red::polish_Tseries.pro
@red::prefilter_data.pro
@red::prepflatcubes_lc4.pro
@red::prepflatcubes.pro
@red::prepmomfbd2.pro
@red::prepmomfbd.pro
@red::quicklook_movie.pro
@red::setflatdir.pro
@red::sumdark.pro
@red::sumflat.pro
@red::sumpinh.pro
@red::sumpolcal.pro
@red::sumprefilter.pro
@red::whichoffset.pro

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
