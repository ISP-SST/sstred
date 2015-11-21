; docformat = 'rst'

;+
; Compute the SST telescope Mueller matrix.
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :author:
; 
;    2012-04-01  Peter Sütterlin, ISP
; 
; 
; :returns:
; 
;    M_TEL:  (float) a 4x4 array with the Mueller Matrix of the telescope 
;
; :Params:
; 
;    lam : in, type=string
;   
;      Wavelength range identifier. Currently allowed values are: 5250
;      5380 6173 6302 7772 8542
;   
;    telpos : in, type="Structure{TIME, AZ, ELEV, TILT}"
;
;      Either an array (as e.g. returned by routine READ_AZEL), in
;      this case TIME also has to be specified and will be used to
;      interpolate the correct values, or an one-element strucure with
;      already interpolated data. In this case, TELPOS.TIME is not
;      used and may be empty
;   
;   
;   
;    time : in, optional, type="string or float"
;   
;      Either (string) in format 'HH:MM:SS.sss' or (float) in seconds
;      since midnight.
;   
; 
; :Keywords:
; 
;    no_zero_offset : in, optional, type=boolean
;   
;      If set, do not set the angle offset between polarimeter
;      coordinates and telescope coordinates to zero. See RESTRICTIONS
;      below.
;      
;    old_polcal :  in, optional, type=boolean
;
;      If the focal plane Stokes vector had been computed using the
;      old polcal whis has opposite handedness, set this keyword to
;      convert to the proper orientation.
;   
;    no_solar_north : in, optional, type=boolean
;   
;      If set, do not correct for the rotation angle between last
;      mirror and solar north direction.
;   
; :restrictions:
;
;    The mount of the linear polarizer (LP) of the calibration unit
;    had been lose, definitely in 2011, probably also in 2010.
;    Therefore all model calculations have been made using the quarter
;    wave plate as a reference, which results in a 10.7 degree offset
;    between polarimeter and telescope coordinates. For normal
;    operation with a well-aligned LP this angle is zero. But for 2011
;    (and maybe 2010) data one should use the QWP as reference and
;    then call red_telmat with NO_ZERO_OFFSET=1
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2013-07-12 : MGL. Use red_m_*.pro rather than m_*.pro.
;                These are Pit's Mueller matrix routines. Also
;                use red_time_conv rather than time_conv.
;
;   2015-11-19 : JdlCR added new telescope parameters from
;                Pit's telcal fits. These are preliminary values!
; 
;-
function red_telmat, lam, telpos, time, NO_ZERO_OFFSET=nzo, OLD_POLCAL=oldpc, $
                  NO_SOLAR_NORTH=nsn, YEAR=year
  
  on_error, 2
  If(n_elements(year) eq 0) then year=2015
  IF n_params() LT 2 THEN $
     message, 'Usage is M_TEL=RED_TELMAT(LAMBDA, TELPOS [, TIME]'
  
  IF size(telpos, /type) NE 8 THEN $
     message, 'TELPOS must be a structure {TIME, AZ, ELEV, TILT}'
  
  nt = n_elements(telpos)
  IF nt EQ 1 THEN $
     tp = telpos $
  ELSE BEGIN
     IF N_PARAMS() NE 3 THEN $
        message, 'If TELPOS is an array you also must specify the time'
     tp = {TIME: 0., AZ: 0., ELEV:0., TILT:0.}
     IF size(time, /type) EQ 7 THEN tim = red_time_conv(time) ELSE tim = time
     tp.time = tim
     tp.az = interpol(telpos.az, telpos.time, tim)
     tp.elev = interpol(telpos.elev, telpos.time, tim)
     tp.tilt = interpol(telpos.tilt, telpos.time, tim)
  ENDELSE
  
  lam_s = strtrim(lam, 2)
  if(year lt 2015) then begin
     CASE lam_s OF
        '5250': BEGIN
           ;; combined telcal of 2011.06.21 and 2011.10.13
           ;; 2011.05.16 data not used (bad calibration polarizer)
           par = [0.96833,  0.00393, -0.06554,  0.98679,  0.00903,  $
                  0.93720, 19.95967,  0.99967, -1.61346,  0.99923,  1.36401, $
                  45.00000,  0.00000, -37.10000, 10.70000]
           
        END
        '5380': BEGIN
           ;; 2011.05.16 data not used (bad calibration polarizer)
           ;; data of 2011.06.21 is not fittable!? use old model....
           par = [0.985871, 0.00304087, -0.00675178, 0.982882, -0.0149910, $
                  0.928876, 21.4228, 1.00000, -0.377807, 1.00000, 1.77330, $
                  45.00000,  0.00000, -37.10000, 10.70000]
        END
        '6173': BEGIN
           ;; data of 2011.10.13;  Used full day = bad fit, but only morning
           ;; is equally bad...
           par = [0.98267,  0.00512, -0.07116,  0.99895, -0.02139,  $
                  0.93017, 17.30311,  0.99657, -0.90633,  0.99882, -0.00390, $
                  45.00000,  0.00000, -37.10000, 10.70000]
           
        END
        '6302': BEGIN
           ;; combined telcal of 2011.05.16, 2011.06.21 and 2011.10.13
           par = [0.97584,  0.00279, -0.03748,  0.99020, -0.00982,  $
                  0.93092, 17.60193,  1.00000,  0.11770,  1.00000,  0.12265, $
                  45.00000,  0.00000, -37.10000, 10.70000]
        END
        '7772': BEGIN
           par = [0.97522,  0.00462, -0.01387,  0.98501, -0.06892,  $
                  0.93033, 20.24933,  0.99802, -0.25072,  0.99938, -0.29458, $
                  45.00000,  0.00000, -37.10000, 10.70000]
        END
        '8542': BEGIN
           ;; combined telcal of 2011.05.16, 2011.06.21 and 2011.10.13
           par = [0.80000,  0.00199, -0.01863,  0.80652, -0.01026,  $
                  0.89742, 14.94485,  1.00000,  1.32428,  1.00000, -1.35152, $
                  45.00000,  0.66346, -37.1, 10.7]
           
        END
        ELSE: BEGIN
           message, 'No information for wavelength '+lam, /continue
           return, diag_matrix([1.,1.,1.,1.])
        END
     ENDCASE
  endif else begin
     print, "red_telmat: Using new telescope parameters, year -> ", year
     CASE lam_s OF
        '5250': BEGIN
           par = [0.989604, 0.00323896, -0.0467535, 0.999643, $
                  0.00775884, 0.938265, 21.8072, 0.993518, -4.44205, $
                  0.999937, 4.37446, -3.40605, 0.0143140, -36.2724, 50.7530]
           
        END
        '5876': BEGIN
           par = [0.991903, 0.00292893, -0.0392169, 0.999995, $
                  0.000716989, 0.934077, 19.3976, 0.993662, -4.29518, $
                  0.997367, 3.73832, -39.9867, 0.00730105, -36.2691, 50.3306]
        END
        '5896': BEGIN
           par = [0.994311, 0.00177954, -0.0363199, 1.00000, $
                  0.00791549, 0.932689, 19.5200, 1.00000, -4.84623, $
                  0.993829, 4.30628, 33.6838, 0.0113876, -35.9496, 51.1442]
        END
        '6173': BEGIN
           par = [0.984056, 0.00251833, -0.0378189, 0.990263, $
                  0.00548041, 0.932719, 18.4015, 0.996951, -3.89294, $
                  0.998076, 3.60348, -24.4415, 0.00405037, -36.0837, 50.7645]
           
        END
        '6302': BEGIN
           par = [0.991819,  0.000789474, -0.0384591,  0.998238, $
                  0.00443999,  0.930898,  17.9770,  0.999954, -2.83140, $
                  0.996109,  2.49997,  -33.9575,  0.00674379, -36.1642,  50.3226]
        END
        '7772': BEGIN
           par = [0.840312, 0.000837646, -0.0223084, 0.842277, $
                  0.00190706, 0.903898, 14.9221, 1.00000, -9.74332, $
                  0.993433, 9.81547, -35.4966, 0.0710514, -35.9800, 50.3051]
        END
        '8542': BEGIN
           par = [ 0.800000, 0.00146548, -0.0182155,  0.802096, $
                   -0.00672012, 0.897333, 15.2085, 0.992795, -2.45377,  $
                   1.00000, 2.46583, -7.62938, 0.687400, -36.7221, 49.5955]
           
        END
        ELSE: BEGIN
           message, 'No information for wavelength '+lam, /continue
           return, diag_matrix([1.,1.,1.,1.])
        END
     ENDCASE



  endelse
  ;; comment: Parameter 14 is the orientation missmatch between the
  ;; polcal-corrected Q direction and the north direction. The way I
  ;; do the fit (use QWP as reference) gives ~10.7 degrees, but for
  ;; data done with a normal polcal one should set it to zero...
  
  IF NOT keyword_set(nzo) THEN par(14) = 0.
  
  m_tel = red_m_sst2(tp.az, tp.elev, par, /deg)
  
  IF keyword_set(oldpc) THEN $
     m_tel = m_tel#red_m_free_mirror(1., 0.)
  
  IF NOT keyword_set(nsn) THEN $
     m_tel = red_m_rot(tp.tilt, /deg)#m_tel
  
  return, m_tel
  
end
