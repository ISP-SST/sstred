; docformat = 'rst'

;+
; Warn if any calibration data is missing or incomplete.
; 
; :Categories:
;
;    SST observations
; 
; 
; :Author:
; 
;    Mats Löfdahl, 2015-09-01
; 
; 
; :Params:
; 
; :Keywords:
; 
;
; :History:
;
;    2015-09-01 : MGL. Looking for different kinds of data ready.
;                 Checking science data for corresponding dark frame
;                 data ready.
;
;    2015-09-03 : MGL. Checking science data for corresponding flat
;                 field data ready.
;
;    2015-09-04 : MGL. Checking science data for corresponding polcal
;                 data ready. Checking science data for corresponding
;                 pinhole data ready.  
;
;    2015-09-07 : MGL. Moved most of the code from this method to a
;                 stand-alone procedure, red_check_calibrations. This
;                 method is now only responsible for getting info from
;                 the config file and then calling the procedure.
;
; 
;-
pro red::check_calibrations, all = all $
                             , darks = darks $
                             , flats = flats $
                             , polcal = polcal $
                             , pinholes = pinholes $
                             , logfile = logfile

  root_dir = self.root_dir

  red_check_calibrations, root_dir $
                          , all = all $
                          , darks = darks $
                          , flats = flats $
                          , polcal = polcal $
                          , pinholes = pinholes $
                          , logfile = logfile


end
