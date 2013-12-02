function dual_fpi, fpi, lam, angle, erh = erh, erl = erl,$
                   ech = ech, ecl = ecl, thr = thr, tlr = tlr
;
; FUNCTION dual_fpi
;
;    INPUT:
;           fpi: structure created with get_fpi_par
;           lam: wavelength array [Anstroms] (absolute wavelenth)
;         angle: angle of incidence of the beam. 
;                [0 for perpendicular incidence]
;
;   OUTPUT:
;           It outputs the transmission profile of CRISP.
;           Thr = HRE individual profile (optional output)
;           Lhr = LRE individual profile (optional output)
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
; Jaime de la Cruz Rodriguez (ISP-KVA 2010)
;  
; --
;
  if not keyword_set(erh) then erh = 0.d0
  if not keyword_set(erl) then erl = 0.d0
  if not keyword_set(ech) then ech = 0.d0
  if not keyword_set(ecl) then ecl = 0.d0
;
  if(n_params() lt 3) then begin
     print, "Incorrect number of parameters"
     return, 0
  endif
; 
; Reflectivities + errors
;
  mrhr = fpi.rhr + erh
  mrlr = fpi.rlr + erl
;
; Finesse of the HRE and LRE
;
  fhr = 4.d0 * mrhr / (1.d0 - mrhr)^2.
  flr = 4.d0 * mrlr / (1.d0 - mrlr)^2.
;
; Phase difference for each etalon
;
  phr = fpi.pi2 * (fpi.shr) * cos(angle)
  plr = fpi.pi2 * (fpi.slr) * cos(angle)
;
; Transmission peaks
;
  Thr = 1.d0 / (1.d0 + fhr * (sin(phr / (lam + ech)))^2.)
  Tlr = 1.d0 / (1.d0 + flr * (sin(plr / (lam + ecl)))^2.)
;
  return, thr * tlr
end
