; docformat = 'rst'

;+
; The solar spectrum from the satlas.
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;    Jaime de la Cruz
; 
; 
; :Returns:
; 
;     The solar spectrum from the satlas.
; 
; :Params:
; 
;     xstart : in, type=float
;
;        The lower wavelength limit.
;
;     xend : in, type=float
;
;        The upper wavelength limit.
;
;     outx : out
;
;        The wavelength grid.
;
;     outy : out
;
;        The spectrum.
; 
; :Keywords:
; 
;    nm : in, optional, type=boolean
;     
;       Set this to indicate that xstart and xend are given in
;       nanometers. 
;   
;    nograv : in
;     
;   
;   
;    nocont : in
;     
;   
;   
;    cgs : in, optional, type=boolean
;     
;       Set this to get output in CGS units.
;   
;    cont : out
;     
;   
;   
;    si : in, optional, type=boolean
;     
;       Set this to get output in SI units.
;     
;    conversion_cgs : out, optional, type=fltarr
;
;       The conversion factor that transforms the spectrum (and cont)
;       to CGS units.
;
;    conversion_si : out, optional, type=fltarr
;
;       The conversion factor that transforms the spectrum (and cont)
;       to SI units.
;   
; 
; 
; :History:
; 
;    2013-07-11 : MGL. Renamed from satlas. Modified mechanism for
;                 finding the data file within the path. Should now
;                 pick the version that is in the same directory as
;                 this file, that is, in the crispred repository.
; 
;    2016-11-11: JdlCR, added the possibility to output CGS/SI units
;                for the calibration of CHROMIS
; 
;    2016-11-29 : MGL. New keywords cgs_conversion and si_conversion.
;
;    2016-12-04 : JdlCR. Removed new keywords, they break the routine.
;
;    2016-12-04 : MGL. Renamed the keywords.
;
;    2016-12-09 : JdlCR. Fixed double bug in the computation of
;                 "conversion_si". 
;
;   2017-04-18 : MGL. Remove the "conversion" keywords.
; 
;
;-
pro red_satlas, xstart, xend, outx, outy $
                , nm = nm $
                , nograv = nograv $
                , nocont = nocont $
                , cgs = cgs $
                , cont = con $
                , si = si

  ;; Find the input data
  this_dir = file_dirname( routine_filepath("red_satlas"), /mark )
  restore, this_dir+'ftsatlas.idlsave'

  if keyword_set(nm) then begin
     xstart=xstart/10.d0
     xend=xend/10.d0
  endif

  c=299792458.d0                ; light speed in m/s

  ;;pos=where(XL_FTS gt xstart AND XL_FTS lt xend)
  ;;outx=xl_fts;[pos]
  if not keyword_set(nograv) then xl_fts *= (1.d0-633.d0/c)
  pos = where(xl_fts ge xstart and xl_fts le xend)
  outx = xl_fts[pos]
  outy = yl_fts[pos]
  con = cint_fts[pos]
  
  if keyword_set(si) then begin

    clight = 2.99792458d8       ;speed of light [m/s]                                  
    aa_to_m = 1d-10                                                                        
    cm_to_m = 1d-2     
    m_to_cm = 1d2

    ;; To Watt/(s m2 Hz ster)
    conversion_si = (outx*aa_to_m)^2 / (clight * cm_to_m^2 * aa_to_m )

    outy *= conversion_si
    con  *= conversion_si

    return

  endif
  
  if keyword_set(cgs) then begin

    clight = 2.99792458d10      ;speed of light [cm/s]
    joule_2_erg = 1d7
    aa_to_cm = 1d-8

    ;; From Watt /(cm2 ster AA) to erg/(s cm2 ster cm)
    conversion_cgs = joule_2_erg / aa_to_cm      
    ;; To erg/
    conversion_cgs *= (outx*aa_to_cm)^2 / clight 

    outy *= conversion_cgs
    con  *= conversion_cgs

    return

  endif
      
  if ~keyword_set(nocont) then begin
    outy /= con
    con[*] = 1.d0
  endif

end
