function get_fpi_trans, fpi, lam, ech = ech, ecl = ecl, erh = erh, erl = erl
;
; FUNCTION get_fpi_trans
;
;    INPUT:
;           fpi: structure created with get_fpi_par
;           lam: wavelength array [Anstroms] (absolute wavelenth)
;
;   OUTPUT:
;           It outputs the transmission profile of CRISP accounting 
;           for the converging beam.
;
; KEYWORDS: 
;           erh = HRE reflectivity error
;           erl = LRE reflectivity error
;           ech = HRE cavity error [Angstroms]
;           ecl = LRE cevity error [Angstroms] 
;
;   OUTPUT:
;           it returns a structure with the nominal reflectivities
;           and cavity separations of the CRISP spectropolarimeter.
;
; Adapted from Goran Scharmer's dual_fpi2.ana by 
; Jaime de la Cruz Rodriguez (ISP-KVA 2010)
;  
; --
;
;
; MODIFICATIONS:
;
;   2010-06-30, J. de la Cruz R.: 
;               Improved quadrature for the integration with only three points.
;
; --
;
; some parameters
;
  if not keyword_set(ech) then ech = 0.d0
  if not keyword_set(ecl) then ecl = 0.d0
  if not keyword_set(erh) then erh = 0.d0
  if not keyword_set(erl) then erl = 0.d0
;
; reflectivities + error
;

  mrhr = fpi.rhr + erh
  mrlr = fpi.rlr + erl
;
; Coefficient of finesse
;
  fhr = 4. * mrhr / (1. - mrhr)^2.
  flr = 4. * mrlr / (1. - mrlr)^2.
;
; Angle averaging for converging beam at fpi.fr
;
  Tav = dblarr(n_elements(lam))
;
; Integration
;
  for n = 0, n_elements(fpi.wng) - 1 do begin
     phr = (fpi.shr) * fpi.calp[n]
     plr = (fpi.slr) * fpi.calp[n]
;
     Tav+= (1. / (1. + fhr * (sin(phr / (lam + ech)))^2.)) * $
           (1. / (1. + flr * (sin(plr / (lam + ecl)))^2.)) * $
           fpi.wng[n]
  endfor
;
  return, Tav
end
