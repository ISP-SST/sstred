; docformat = 'rst'
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
; :Author:
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
; :Returns:
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
; 
; :Keywords:
; 
;    develop : in, optional, type=boolean
; 
;       Run in developer mode.
;
;    no_db : in, optional, type=boolean
;
;       Do not use metadata database.
; 
; :History:
; 
;   2012-03-08 : JdlCR, The c++ routine red_cconvolve seemed to shift
;                the convolved data by dx,dy > 0.5 < 2 pixels.
;                Re-using the old IDL convolution routine (much
;                slower, yet it works better). Must check the c++
;                version.
;
;   2012-03-21 : JdlCR, red_matrix2momfbd -> Corrected clips,
;                pol::demodulate -> some extra features like
;                (/noclip)
;  
;   2012-05-26 : JdlCR, added new IDL polcal routines. Using my C++
;                with Pit's modifications and his interface.
;
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2013-06-09 : MGL. Added red::setpinhdir.pro to list of included
;                files.
; 
;   2013-09-04 : MGL. Added red::pinholecalib.pro to list of included
;                files.
; 
;   2013-09-06 : MGL. Added red::prepmfbd.pro to list of included
;                files.
; 
;   2013-12-10 : PS  adapt for multiple flat_dir
;
;   2016-04-29 : THI. Split class RED into a base-class (instrument
;                independent parts) and derived classes
;                (CRISP/CHROMIS).
;
;   2017-03-09 : MGL. New keyword "develop".
;
;   2021-03-03 : MGL. New keyword no_db.
; 
;-
function crispred, filename, develop = develop, no_db = no_db

  return, obj_new('crisp', filename, develop = develop, no_db = no_db)
  
end
