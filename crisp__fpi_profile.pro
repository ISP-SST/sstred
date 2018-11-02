; docformat = 'rst'

;+
; Calculates the transmission profile of the CRISP FPIs.
;
; Based on dual_fpi2.ana by Göran Scharmer (2008).
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
;    fpi : in, type=struct
;   
;      Structure created with red_get_fpi_par
;   
;    lam : in, type=fltarr
;   
;      wavelength array [Angstroms] (absolute wavelenth)
;   
; 
; :Keywords:
; 
;    ech : 
;   
;      HRE cavity error [Angstroms]
;   
;    ecl : 
;   
;      LRE cevity error [Angstroms] 
;   
;    erh : 
;   
;      HRE reflectivity error
;   
;    erl : 
;   
;      LRE reflectivity error
;   
; 
; :History:
; 
;   2008 : Scharmer. Original dual_fpi2.ana.
; 
;   2010 : JdlCR. Ported to IDL as red_get_fpi_trans,pro.
;
;   2010-06-30 : J. de la Cruz R.: Improved quadrature for the
;                integration with only three points.
;
;   2018-05-24 : MGL. Made into a class method. Incorporated the
;                relevant parameter settings from red_get_fpi_par. 
; 
;-
function crisp::fpi_profile, w, prefilter $
                             , ech = ech, ecl = ecl, erh = erh, erl = erl $
                             , offset_correction = offset_correction

  ;; Defaults
  if n_elements(ech) eq 0 then ech = 0.d0
  if n_elements(ecl) eq 0 then ecl = 0.d0
  if n_elements(erh) eq 0 then erh = 0.d0
  if n_elements(erl) eq 0 then erl = 0.d0

  ;; Parameters from get_par
  
  ;; Quadrature for the integration of Tr
  calp = cos(sqrt([0.1127017d, 0.5000000d, 0.8872983d] * (0.5d/165.d)^2.)) * 2d*!dpi 
  wng = [0.2777778, 0.4444444, 0.2777778] ; weights of the quadrature points

  ;; CRISP parameters
;  fr = 165.0d0                  ; F number for the incident beam
  shr = 787.e4                  ; HRE separation 
  slr = 295.5e4                 ; LRE separation
  w0 = double(prefilter)        ; (Default) filter center wavelength, see CASE below

  case prefilter of
    '8542': begin
      rhr  = 0.930d0            ; HRE reflectivity
      rlr  = 0.852d0            ; LRE reflectivity
    end
    '6302': begin
      rhr  = 0.935d0
      rlr  = 0.838d0
    end
    '5380': begin
      rhr  = 0.9375d0
      rlr  = 0.8510d0
    end
    '5250': begin
      w0   = 5250.5d0
      rhr  = 0.9589d0
      rlr  = 0.8561d0
    end
    '5576': begin
      w0   = 5576.5d0
      rhr  = 0.9457d0
      rlr  = 0.8903d0
    end
    '5875': begin
      w0   = 5876.0d0
      rhr  = 0.9266d0
      rlr  = 0.8567d0
    end
    '6563': begin
      rhr  = 0.9330d0
      rlr  = 0.8492d0
    end
    '6082': begin
      rhr  = 0.9381d0
      rlr  = 0.8447d0
    end
    '6173': begin
      rhr  = 0.9380d0
      rlr  = 0.8469d0
    end
    '7090': begin
      rhr  = 0.9430d0
      rlr  = 0.8565d0
    end
    '5896': begin
      rhr  = 0.9270d0
      rlr  = 0.8545d0
    end
    '5173': begin
      rhr  = 0.9585d0
      rlr  = 0.8654d0
    end
    '7772': begin
      rhr  = 0.9227d0
      rlr  = 0.8576d0
    end
  endcase
  
  ;; END Parameters from get_par

  ;; Refine nominal cavity separations, so peaks are aligned at w0.
  nhr = long(0.5d0+shr/(w0/2.d0)) 
  shr = nhr*w0/2.d0
  nlr = long(0.5d0+slr/(w0/2.d0)) 
  slr = nlr*w0/2.d0

  
  lam = w + w0

  ;; Reflectivities + error

  mrhr = rhr + erh
  mrlr = rlr + erl

  ;; Finesse
  fhr = 4. * mrhr / (1. - mrhr)^2.
  flr = 4. * mrlr / (1. - mrlr)^2.

  ;; Integration
  Tav = dblarr(n_elements(lam))
  for n = 0, n_elements(wng) - 1 do begin
    
    phr = (shr) * calp[n]
    plr = (slr) * calp[n]

    Tav += (1. / (1. + fhr * (sin(phr / (lam + ech)))^2.)) * $
           (1. / (1. + flr * (sin(plr / (lam + ecl)))^2.)) * $
           wng[n]
    
  endfor                        ; n

  if keyword_set(offset_correction) then begin

    ;; Center the transmission profile, from red::fitprefilter in the
    ;; master branch.
    
    dum = max(Tav, p)
    cc = poly_fit((lam[p-1:p+1]-w0) * 100.d0, Tav[p-1: p+1], 2)
    offset = -0.005d0 * cc[1] / cc[2]

    ;; Re-calculate with lambda offset
    Tav *= 0d
    for n = 0, n_elements(wng) - 1 do begin
      
      phr = (shr) * calp[n]
      plr = (slr) * calp[n]

      Tav += (1. / (1. + fhr * (sin(phr / (lam + offset + ech)))^2.)) * $
             (1. / (1. + flr * (sin(plr / (lam + offset + ecl)))^2.)) * $
             wng[n]
      
    endfor                      ; n


  endif

  
  return, Tav
  
end
