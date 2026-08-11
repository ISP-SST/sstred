FUNCTION red_Time_conv, intime, f2s=f2s, frac=frac
;+
; NAME:
;	red_Time_conv
; PURPOSE:
;	Onvert a Time string HH:MM:SS to seconds or vice versa
; CALLING SEQUENCE:
;	Res = red_Time_conv( Time )
; INPUTS:
;	Time : Either integer value with number of seconds or string
;              with time of the day in format HH:MM:SS
; KEYWORDS:
;       F2S  : Default is convert string to integer, F2S reverts.
;       FRAC : (integer) if specified and non-zero, write fractions of 
;              seconds with FRAC digits
; OUTPUTS:
; 	Res : integer with number of seconds or string in format
;             HH:MM:SS, depending on Keyword F2S
; MODIFICATION HISTORY:
;	01-Apr-1991  P.Suetterlin, KIS
;       17-Dec-2010  PS output fractional seconds
;       25-Dec-2010  PS fix zero padding in fraction
;       2013-07-12 : MGL. Renamed to red_time_conv for inclusion in
;                    crispred pipeline.
;
;-

on_error, 2

IF keyword_set(f2s) THEN BEGIN
    tim1 = '00:00:00'
    z = strtrim(fix(intime/3600), 2)
    IF strlen(z) EQ 1 THEN z = '0'+z
    strput, tim1, z
    z1 = strtrim(fix((intime-z*3600l)/60), 2)
    IF strlen(z1) EQ 1 THEN z1 = '0'+z1
    strput, tim1, z1, 3
    z2 = strtrim(fix(intime-z*3600l-z1*60), 2)
    IF strlen(z2) EQ 1 THEN z2 = '0'+z2
    strput, tim1, z2, 6
    IF keyword_set(frac) THEN BEGIN
        tf = intime-long(intime)
        tf = string(tf, form=string(frac+2, frac, form='("(f",I0,".",I0,")")'))
        tim1 += strmid(tf, 1)
    ENDIF
ENDIF ELSE BEGIN
    tim1 = 3600d*double(strmid(intime, 0, $
                               2))+60d*double(strmid(intime, 3, 2))
    tim1 = tim1+double(strmid(intime, 6))
ENDELSE

return, tim1

END
