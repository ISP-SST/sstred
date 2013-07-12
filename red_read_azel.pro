FUNCTION red_read_azel, file, date, ACTUAL_POSITION=act
;+
; NAME:
;       red_READ_AZEL
; PURPOSE:
;       Read turret logfile of SST, put values into structure
; CALLING SEQUENCE:
;       Res = red_READ_AZEL(Logfile, Date, Time [, /ACTUAL_POSITION]
; INPUTS:
;       Logfile:  (string) name of the turret logfile
;       Date:     (string) date of the day to extract.  Should be YYYY/MM/DD
; KEYWORDS:
;       ACTUAL_POSITION:  (bool) normally the routine reads the nominal
;                         position of the turret (i.e., for disk center).
;                         Deviation from actual is less than hundreths of a
;                         degree.  This data does not have weird jumps (when
;                         stopping the telescope, pointing away from sun etc.)
;                         and thus is safer for interpolation.
;                         Set this keyword if you want to read the actual
;                         turret pointing.
; OUTPUTS:
;       Res:        (struct) {Time, Az, Elev, Tilt}
;       Time:       (float) time axis in seconds since midnight
;       Az:         (float) turret azimuth in degrees
;       Elev:       (float) turret elevation in degrees
;       Tilt:       (float) image rotation angle in degrees
;
; MODIFICATION HISTORY:
;       11-Feb-2011  P.Suetterlin, ISP
;       15-Dec-2011  Change to function delivering a structure
;        9-Feb-2012  Better(?): array of structures
;
;       2013-07-12 : MGL. Renamed to red_read_azel for inclusion in
;                    crispred pipeline.
;
;-

IF keyword_set(act) THEN pos = 2 ELSE pos = 10
IF n_params() LT 2 THEN BEGIN
    date = file_basename(file)
    date = strmid(date, strpos(date, '-/_')+1)
    date = strjoin(strsplit(date, '.', /extr), '/')
ENDIF

openr, unit, file, /get
line = ''
tt = 0.
az = 0.
el = 0.
tilt = 0.
WHILE NOT eof(unit) DO BEGIN
    readf, unit, line
    IF strpos(line, date) NE 0 THEN CONTINUE
      ;;; filter out park position/night date
 ;   IF strpos(line, '05h') GT 0 THEN continue
    tags = strsplit(strcompress(line), ' ', /extract, count=nt)
    IF nt NE 12 THEN CONTINUE   ;;; something's fishy
    tt = [temporary(tt), float(red_time_conv(tags(1)))]
    az = [temporary(az), float(tags(pos))]
    el = [temporary(el), float(tags(pos+1))]
    tilt = [temporary(tilt), float(tags(9))]
ENDWHILE

tt = tt(1:*)
az = az(1:*)
el = el(1:*)
tilt = tilt(1:*)

free_lun, unit
n = n_elements(tt)
res = replicate({Time: 0., Az: 0., Elev: 0., Tilt: 0.}, n)
res.time = tt
res.az = az
res.elev = el
res.tilt = tilt

return, res

END
