; docformat = 'rst'

;+
; Convert Stonyhurst coordinates to Cartesian Helioprojective coordinates.
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, Institute for Solar Physics
; 
; 
; :Returns:
; 
; 
; :Params:
; 
;    isodate : in, type=string
; 
;      The date in ISO format. 
; 
;    stonyhurst_coordinates : in, type="fltarr(2,n) or strarr(n)"
; 
;      Stonyhurst latitude and longitude, N pairs of coordinates to be
;      converted to HPLN,HPLT coordinates. If strings, then on the
;      form "N  2.02  W  0.16", where the spaces are optional and N
;      could be S and W could be E. Stonyhurst coordinates are in
;      degrees. 
; 
;    time : in, type=double
; 
;      The time in seconds after midnight.
; 
; :Keywords:
; 
;    helio : out, optional, type="fltarr(2,n)"
;   
;      The input Stonyhurst coordinates as numbers, in arcsec.
;      Helio[0,*] = E/W = longitude; helio[1,*] = N/S = latitude.
; 
; 
; :History:
; 
;   2021-12-08 : MGL. First version.
; 
;   2021-12-15 : MGL. Use get_sun() from SolarSoft, adapted to SSTRED.
; 
;-
function red_convert_stonyhurst, stonyhurst_coordinates, isodate, time, helio = helio

  
  if size(stonyhurst_coordinates,/tname) eq 'STRING' then begin
    n_points = n_elements(stonyhurst_coordinates)  
    stonyhurst_coordinates_compressed = strcompress(stonyhurst_coordinates,/remove_all)
    helio = stregex(stonyhurst_coordinates_compressed,'[NS]([0-9.]*)[EW]([0-9.]*)',/sub,/extr,/fold)
    helio = float(helio[[2, 1], *]) ; [degrees]
    ;; Helio[0,*] = E/W = longitude = Phi
    ;; Helio[1,*] = N/S = latitude  = Theta    
    negindx = where(strmatch(stonyhurst_coordinates_compressed, '*E*'), count)    
    if count gt 0 then helio[0, negindx] *= -1    
    negindx = where(strmatch(stonyhurst_coordinates_compressed, '*S*'), count)
    if count gt 0 then helio[1, negindx] *= -1
  endif else begin
    helio = stonyhurst_coordinates
    dims = size(stonyhurst_coordinates, /dim)
    n_points = dims[1]
  endelse

  ;; helio is Stonyhurst coordinates in degrees
  
;  ;; Date handling
;  splitdate = strsplit(isodate, '[-.]', /extract)
;  year  = splitdate[0]
;  month = splitdate[1]
;  day   = splitdate[2]
;  
;  ;; Calculate b0 - From pb0r.pro in SolarSoft.
;
;  jd = julday(month, day, year) - 0.5  ; same as jd.int+jd.frac if jd=anytim2jd(yyyy-mm-dd)
;  de = jd - 2415020L + time/(3600.*24) ; Add the time as fractions of 24 h.
;  sunpos, de, ra, dec, longmed
;  lambda = longmed - (20.5D0/3600.0D0)
;  node = 73.666666D0 + (50.25D0/3600.0D0)*( (de/365.25d0) + 50.0d0 )
;  arg = lambda - node
;  b = -asin( 0.12620d0 * sin(arg*!dtor) ) * !radeg ; [degrees] Solar axis tilt angle

  ;; Constants
  r = 6.95508d8                 ; [m] wcs_rsun
  d = 1.49597870691d11          ; [m] wcs_au

  ;; Use SolarSoft's get_sun function (not the same as red_get_sun,
  ;; originally from ANA. Makes some difference!
  tmp = red_solarsoft_get_sun(isodate, dist = d, sd = sd, he_lat = b)
  ;; sd is semidiameter of disk in arc seconds. Calculate rsun in meters:
  r = d * (sd / 3600. * !dtor)  ; [m]

  
  ;; Now follow Thompson (2005) ------------------------------

  Phi   = helio[0, *] * !dtor   ; Stonyhurst latitude in radians
  Theta = helio[1, *] * !dtor   ; Stonyhurst longitude in radians
  
  Phi0 = 0.                     ; Stonyhurst heliographic longitude of the observer  
  b0 = b * !dtor                ; Tilt angle in radians
  
  ;; Heliocentric Cartesian from Stonyhurst - Eq. (11)
  x = r * cos(Theta) * sin(Phi-Phi0)                                      ; [m]
  y = r * ( sin(Theta) * cos(b0) - cos(Theta) * cos(Phi-Phi0) * sin(b0) ) ; [m]

  ;; Using eq 4 or eq 16 doesn't seem to matter.
  if 0 then begin
    ;; HPLN and HPLT in arcsec from Heliocentric Cartesian - Eq. (4)
    theta_x = x/d * !radeg * 3600. 
    theta_y = y/d * !radeg * 3600. 
  endif else begin
    ;; HPLN and HPLT in arcsec from Heliocentric Cartesian - Eq. (16)
    dd = sqrt(x*x + y*y + (d-r)^2)
    theta_x = atan(x, d-r) * !radeg * 3600. ; arg(d-r, x)
    theta_y = asin(y/d)    * !radeg * 3600.
  endelse
  
  hpc = double(helio)
  hpc[0,*] = theta_x
  hpc[1,*] = theta_y
  
  return, hpc

end

logfile = '/scratch/mats/convert_lp/2013-06-26/CRISP/downloads/sstlogs/positionLog_2013.06.26_final'
spawn, 'cat '+logfile, loglines

logitems = (strsplit(loglines,/extract)).toarray()

isodate = '2013-06-26'
time = red_time2double(strmid(logitems[*,1], 11, 8))
stony = strarr(n_elements(loglines))
for i = 0, n_elements(loglines)-1 do stony[i] = strjoin(reform(logitems[i,4:7]))

hpc = red_convert_stonyhurst(stony, isodate, time, helio = helio)

;; Carlos used SunPy to convert the same coordinates
hpc_sunpy = fltarr(2, n_elements(loglines))
openr, 1, '/home/mats/idl/sstred/cartesian-from-stonyhurst-sunpy.txt'
readf, 1, hpc_sunpy
close, 1

have_solarsoft = 0
if have_solarsoft then begin

  ;; If SolarSoft available, try lonlat2xy.pro
  ;; Gives the same results as sunpy
  hpc_lonlat2xy = lonlat2xy(helio, isodate)
  
endif

red_logdata, isodate, disktime, diskpos = diskpos


cgwindow
cgplot, /add,hpc_sunpy[0,*],hpc_sunpy[1,*], psym=16, xtitle = 'HPLN / 1"', ytitle = 'HPLT / 1"', title = 'Helioprojective Cartesian', /yno
cgplot, /add, /over, hpc[0,*],hpc[1,*], psym=16, color = 'red'
if 1 then cgplot, /add, /over, diskpos[0,*],diskpos[1,*], psym=16, color = 'purple'
if have_solarsoft then begin
  cgplot, /add, /over, hpc_lonlat2xy[0,*],hpc_lonlat2xy[1,*], psym=16, color = 'blue'
  cglegend, /add, title = ['SSTRED', 'lonlat2xy', 'SunPy'], colors = ['red', 'blue', 'black'], loc = [0.90, 0.85], align = 1,psym=16,len=0  
endif else begin
  cglegend, /add, title = ['SSTRED', 'SunPy'], colors = ['red', 'black'], loc = [0.90, 0.85], align = 1,psym=16,len=0  
endelse
cgcontrol, out = 'cartesian.pdf'



cgwindow
; Note helio[0,*] is Stonyhurst longitude, [1,*] is latitude, they are
; both in degrees.
cgplot, /add,helio[0,*], helio[1,*], psym=16, xtitle = 'HGLN / 1$\deg$ (E/W)', ytitle = 'HGLT / 1$\deg$ (N/S)', title = 'Heliographic (Stonyhurst)', /yno, color = 'red'
cgcontrol, out = 'stonyhurst.pdf'

end
