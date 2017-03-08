; docformat = 'rst'

;+
; Make a time stamp for the current time.
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, Institute for Solar Physics
; 
;    Parts stolen from The Coyote Graphics timestamp command, see
;    http://www.idlcoyote.com/.
;
; :Returns:
; 
;    A timestamp string formatted like: YYYYMMDD@hh:mm:ss. 
; 
; :Params:
; 
;   
; 
; :Keywords:
; 
;    utc : in, optional, type=boolean
;
;         If this keyword is set, the UTC time will be used for
;         the time stamp instead of the local time.
; 
;    iso : in, optional, type=boolean
; 
;         If this keyword is set, the timestamp will be formatted as
;         YYYY-MM-DDThh:mm:ss. 
; 
; :History:
; 
;   2013-09-02 : MGL. First version.
;
;   2016-05-30 : MGL. New keyword "iso".
; 
; 
;-
function red_timestamp, utc = utc, iso = iso

  time = Systime(UTC=Keyword_Set(utc))
  day = Strmid(time, 0, 3)
  date = String(StrMid(time, 8, 2), Format='(I2.2)') ; Required because UNIX and Windows differ in time format.
  month = Strmid(time, 4, 3)
  year = Strmid(time, 20, 4)
  stamp = Strmid(time, 11, 8)
  months = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']
  m = (Where(months EQ StrUpCase(month))) + 1
  
  if keyword_set(iso) then begin
    timestamp = year + '-' + string(m, FORMAT='(I2.2)') + '-' + date + 'T' + stamp
  endif else begin
    timestamp = year + String(m, FORMAT='(I2.2)') + date + '@' + stamp
  endelse

  return, timestamp

end
