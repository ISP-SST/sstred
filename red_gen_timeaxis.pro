; docformat = 'rst'

;+
; Generate nice hour:minute marks from seconds-since-midnight array.
;
; :Author:
; 
;    23-Jan-2006  P.Suetterlin, SIU
;
; :Returns:
;
;      Structure with needed values to generate plot axis: ticks,
;      minor, tickv, name 
;
; :Params:
;
;       X : in, type=array
;
;           Vector (x-axis) with time as seconds since 00:00 UT
;
;
; :History:
; 
;       2006-01-23 : P.Suetterlin, SIU
;
;       2007-06-15 : Rewrite as case-clause. Take care of series that
;                    are shorter than 2 binwidth
;
;       2009-05-04 : Take care of 24h wrap
;
;       2015-08-12 : MGL. Import into crispred namespace. Quick fix of
;                    bug making x1 < x0 in some cases. 
;
;       2015-08-17 : MGL. Removed the fixing of the number of minor
;                    tickmark intervals to five. The number given in
;                    the case statement is now in force.
;
;-
FUNCTION red_gen_timeaxis, x

;       Based on range set spacing between major tickmarks and number
;       of minor tickmarks.  Pack into a structure.  
;       Example of use:
;          Given a variable Y measured at times X
;          L = RED_GEN_TIMEAXIS(X)
;          PLOT, X, Y, XTICKV=l.tickv, XTICKS=l.ticks, XMIN=l.minor, XTICKNAM=l.name

  mi = min(x, max=ma)
  range = ma-mi
  print, 'range:', range
  CASE 1 OF
    range LE 2:  BEGIN
      step = 1l
      minor = 10 
    END
    range LE 5:  BEGIN
      step = 1l
      minor = 2 
    END
    range LE 15:  BEGIN
      step = 5l
      minor = 5 
    END
    range LE 30:  BEGIN
      step = 10l
      minor = 10 
    END
    range LE 120:  BEGIN
      step = 30l
      minor = 3 
    END
    range LE 300:  BEGIN
      step = 60l
      minor = 6 
    END
    range LE 600:  BEGIN
      step = 60l
      minor = 6 
    END
    range LE 1500:  BEGIN
      step = 300l
      minor = 5
    END
    range LE 2700:  BEGIN
      step = 600l
      minor = 5
    END
    range LE 5400:  BEGIN
      step = 900l
      minor = 3
    END
    range LE 10800: BEGIN
      step = 1200l
      minor = 4
    END
    range LE 21600: BEGIN
      step = 1800l
      minor = 3
    END
    Else:           BEGIN
      step = 3600l
      minor = 4
    END
  ENDCASE

  x0 = fix(mi/step)
  x0 += (x0 EQ mi/float(step) ? 0:1)
  x1 = fix(ma/step)
  if x1 lt x0 then begin
    temp = x0
    x0 = x1
    x1 = temp
  end
  ticks = x1-x0
  IF ticks EQ 0 THEN BEGIN
          ;;; IDL won't do it correctly with only one tickmark...
    IF abs(mi-step*(x0-1)) LT abs(step*(x1+1)-ma) THEN BEGIN
      x0 -= 1
      ticks = 1
    ENDIF ELSE BEGIN
      x1 += 1
      ticks = 1
    ENDELSE
  ENDIF
  tickv = (lindgen(ticks+1)+x0)*step
;minor = 5
;  IF step EQ 3600 THEN $
;     name = strtrim((tickv/3600) MOD 24, 2) $
;  ELSE $
;     name = strtrim((tickv/3600) MOD 24, 2)+':'+red_nnumber((tickv-tickv/3600*3600)/60, 2)

  h = tickv/3600
  m = (tickv-h*3600)/60
  s = tickv-h*3600-m*60

  h_tags = strtrim(h MOD 24, 2)
  m_tags = string(m, format = '(i02)')
  s_tags = string(s, format = '(i02)')
  
  case 1 of
    step eq 3600 : name = h_tags
    step ge 60   : name = h_tags + ':' + m_tags
    else         : name = h_tags + ':' + m_tags + ':' + s_tags
  endcase 

;  case 1 of
;    step eq 3600 : name = strtrim((tickv/3600) MOD 24, 2)
;    step ge 60 : name = strtrim((tickv/3600) MOD 24, 2)+':'+red_nnumber((tickv-tickv/3600*3600)/60, 2)
;    else : name = strtrim((tickv/3600) MOD 24, 2)+':'+red_nnumber((tickv-tickv/3600*3600)/60, 2)
;  endcase 
  
  return, {ticks: ticks, $
           Minor: minor, $
           Tickv: tickv, $
           Name:  name}

END

time = [8.5, 9.8]*3600
time = [8.45, 8.55]*3600
time = [8, 8]*3600 + [26.42, 26.44]*60
print, red_timestring(time)

red_timeplot, time, 0*time+1, yrange = [0, 2]


end
