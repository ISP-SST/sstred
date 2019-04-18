; docformat = 'rst'

;+
; Find the zero point of CHROMIS scans with wavelength given in HRE
; digital units.
;
; The zero points are stored in files to be used by the extractstates
; method. 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl
; 
; 
; 
; :History:
;
;    2016-09-18 : MGL. First version.
;
;    2016-09-21 : MGL. Work in meters, not Å or mÅ. Bugfix in
;                 extracting the hrz tunings. Typos hzr --> hrz.
;                 Remove previously generated hrz_zeropoint files.
;                 Change filter tags to four-digit tags representing
;                 approximate filter wavelength.
;    
;    2017-08-18 : THI. Use the FPI calibration values to get hrz instead
;                 of parsing files.
; 
;    2019-04-18 : THI. Make the class/method version call the procedure version
;                 instead of keeping duplicate code in 2 files.
;
;-

pro chromis::hrz_zeropoint

  chromis_hrz_zeropoint, self.out_dir
