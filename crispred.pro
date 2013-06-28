; docformat = 'rst'

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
@red::getalignclips
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

;+
; Class crispred and subroutines, class polarim and subroutines.
; 
; Reduction pipeline for the SST. Steps prior to momfbd. The "pol"
; class takes care of the demodulation of the momfbd-ed data.
;
; :Categories:
;
;    CRISP pipeline
; 
; 
; :author:
; 
;   Jaime de la Cruz Rodriguez (IFA-UU 2011),
;         Department of Physics and Astronomy, Uppsala University,
;         jaime.cruz@physics.uu.se, jaime@astro.uio.no
;
;   Michiel van Noort (Max Plank - Lindau):
;         red_matrix2momfbd (adapted),
;         C++ polcal curves (adapted),
;         MOMFBD image format DLMm
;         fillpix adapted from Mats' (IDL) and Michiel's (c++)
;
;   Pit Suetterlin (ISP-KVA):
;         red::polcal (interface for the C++ code),
;         red::getalignclips,
;         red::getoffsets,
;         and dependencies
;
;   Tomas Berger (LMSAL):
;         destretch IDL module and routines.
;
;   Mats Löfdahl (ISP-KVA):
;         red::taper,
;         red::offsets,
;         red_findpinholegrid,
;         fillpix adapted from Mats' (IDL) and Michiel's (c++),
;         shift-and-sum in red_sumfiles for summing pinholes
;
; 
; :returns:
; 
; 
; :Params:
; 
;   filename : in, optional, type=string, default="config.txt"
;   
;     The name of the configuration file, that describes the location
;     of the science and calibration data (and a few other things).
;   
; 
; :Keywords:
; 
; 
; :dependencies:
;
;    The external C++ module must be compiled and the system variable
;    CREDUC point to creduc.so.
;
; 
; :history:
; 
;    2012-03-08 : JdlCR, The c++ routine red_cconvolve seemed to shift
;                 the convolved data by dx,dy > 0.5 < 2 pixels.
;                 Re-using the old IDL convolution routine (much
;                 slower, yet it works better). Must check the c++
;                 version.
;
;   2012-03-21 : JdlCR, red_matrix2momfbd -> Corrected clips,
;                pol::demodulate2 -> some extra features like
;                (/noclip)
;  
;   2012-05-26 : JdlCR, added new IDL polcal routines. Using my C++
;                with Pit's modifications and his interface.
;
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
; 
;-
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
