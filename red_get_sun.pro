pro red_get_sun,iy,im,id,sday,ha,dec

;RETURNS THE HOUR ANGLE AND DECLINATION OF THE SUN AT LA PALMA SSO
;CONVERTED FROM ANA ROUTINE LAPALMASUN.ANA

;  time variable for Newcomb's folmulae:
;  fraction of Julian centuries elapsed since 1900.05 (noon 1st of January)
;  = J.D.2415020.0) 
  jd=julday(im,id,iy)
  h0=(jd-2415020.0d0)/36525.0d0
  h=(jd-2415020.0d0+double(sday)/86400.D0)/36525.0d0
  hh=h*h
;  Newcomb's formulae. (page 98 Explanatory suppl. to the ephemeris)
;  mean obliquity of the ecliptic
  ehel=(0.4093198d0)-(2.2703d-4)*h-(2.86d-8)*hh
;  eccentricity of earth's orbit
  eks=(0.01675104d0)-(0.0000418d0)*h
;  mean longitude of sun
  sml=279.6967d0+36000.769d0*h
  sml=sml mod 360.D0
  sml=sml*!pi/180.0d0
;  mean anomaly
  anm=358.4758d0+35999.0498d0*h-0.00015d0*hh
  anm=anm*!pi/180.0d0
;  true longitude of sun (sl)
  cc=(1.91946d0-0.00479d0*h)*sin(anm)+0.020d0*sin(2.d0*anm)
  cc=cc*!pi/180.d0
  sl=sml+cc
;  declination of apparent sun
  dec=asin(sin(ehel)*sin(sl))
;  right ascension of apparent sun
  ra=atan(  ( cos(ehel)*sin(sl) ) ,cos(sl)  )
;  convert ra in radians to hours
  ra=ra*24.0d0/(2.d0*!pi)
  ra=ra+24.0d0*(ra lt 0.)
  dec=dec*360.0d0/(2.d0*!pi)
;  sidereal time in Greenwich at 0 UT (Newcomb)
  sidg0=6.6461D0+2400.0513D0*h0
  sidg0=sidg0 mod 24.D0
;  sidereal time in Greeenwich at any instant
  sidg=(double(sday)*1.002737908D0/3600.0D0)+sidg0  
;  longitude of observatory in degrees, negative if west
  lo=-17.880d0
;  convert longitude of obs (degrees) to time measure and find local sidereal 
;  time, longitude is positive when west
  sidl=sidg+lo*(24.0d0/360.0d0)
  ha=sidl-ra
;restrict to the -12 to +12 range
  ha=(ha+12.d0) mod 24.D0 - 12.d0
  ha=(ha-12.d0) mod 24.D0 + 12.d0
  return
end
