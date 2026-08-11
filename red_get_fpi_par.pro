; docformat = 'rst'

;+
;    Returns a structure with the nominal reflectivities and cavity
;    separations of the CRISP spectropolarimeter.
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;    Jaime de la Cruz Rodriguez, ISP-KVA, 2010
; 
; 
; :Returns:
; 
; Creates the structure 'fpi':
; 
;     w0:   nominal center of the trans peak
; 
;    rhr:   reflectivity hre
; 
;    rlh:   reflectivity lre
; 
;    shr:   separation of the hre
; 
;    slr:   separation of the lre
; 
;     fr:   f number for the incident beam (f=165 for CRISP)
; 
;    prw:   prefilter FWHM
; 
;   pfw0:   prefilter central wavelength
; 
;   calp:   quadrature for the integration of Tr
; 
;    wng:   weights of the quadrature points
; 
; :Params:
; 
; 
; :Keywords:
; 
;   line : type=string
;   
;   
;   
; 
; 
; :History:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2013-07-11 : MGL. Use red_get_fpi_sep, not get_fpi_sep.
;
;   2016-09-25 : MGL. Add CHROMIS prefilters.
; 
; 
;-
function red_get_fpi_par, line = line

  ;;pi2 = 6.28318530779d0         ; 2.d0 * PI
  pi2 = !dpi*2d
  calp = cos(sqrt([0.1127017, 0.5000000, 0.8872983] * (0.5/165.)^2.)) * pi2
  wng = [0.2777778, 0.4444444, 0.2777778]
                                ;
  fpi=create_struct('rhr',0.d0, $           ; reflectivity hre
                    'rlr',0.d0, $           ; reflectivity lre
                    'slr',0.d0, $           ; separation of the lre
                    'shr',0.d0, $           ; separation of the hre
                    'dhr',0.d0, $           ;
                    'dlr',0.d0, $           ;
                    'w0',0.d0, $            ; 
                    'fr',165.0d0, $         ; F number for the incident beam (F=165 for CRISP)
                    'thmax',0.5d0/165.d0, $ ;
                    'prw',0.d0, $           ; prefilter FWHM
                    'prw0',0.d0, $          ; prefilter central wavelength
                    'prnc',2.e0, $          ;
                    'fpi_wav',0.d0, $       ;
                    'rng', 0.0d0, $         ;
                    'calp', calp, $         ; quadrature for the integration of Tr
                    'wng', wng, $           ; weights of the quadrature points
                    'pi2', pi2 $            ; 2*pi in double precision
                   )

  lines = ['3999', '3969', '3978', '3934', '3925', '4862' $   ; CHROMIS
           , '5173', '5380', '5250', '5576', '5896', '6082' $ ; CRISP
           , '6173', '6302', '6563', '7090', '7772', '8542']  ; CRISP

  if not keyword_set(line) then begin
    print, 'Usage: fpi = red_get_fpi_par(line = <line>)'
    print, '<line> (string): ',lines
    return,0
  endif

  ;; F-ratios for the two beams:
  fr_crisp   = 165.0d0
  fr_chromis = 165.0d0                              ; FAKE DATA ******************************
  ;; We could also allow sending it in here in a keyword. The calling
  ;; method should be able to figure this out from the pinhole data.

  case line of
    '8542': begin
      fpi.w0   = 8542.d0
      fpi.rhr  = 0.930d0
      fpi.rlr  = 0.852d0
      fpi.shr  = 787.e4
      fpi.slr  = 295.5e4
      fpi.dhr  = 0.363d0
      fpi.dlr  = 0.991d0
      fpi.prw  = 10.0d0
      fpi.prw0 = 8542.0d0
      fpi.rng  = 0.7d0
      fpi.fr   = fr_crisp
    end
    '6302': begin
      fpi.w0   = 6302.d0
      fpi.rhr  = 0.935d0
      fpi.rlr  = 0.838d0
      fpi.shr  = 787.e4
      fpi.slr  = 295.5e4
      fpi.dhr  = 0.267d0
      fpi.dlr  = 0.731d0
      fpi.prw  = 4.37d0
      fpi.prw0 = 6302.d0
      fpi.rng  = 0.35d0
      fpi.fr   = fr_crisp
    end
    '5380': begin
      fpi.w0   = 5380.0d0
      fpi.rhr  = 0.9375d0
      fpi.rlr  = 0.8510d0
      fpi.shr  = 787.e4
      fpi.slr  = 295.5e4
      fpi.dhr  = 0.228d0
      fpi.dlr  = 0.624d0
      fpi.prw  = 3.4d0
      fpi.prw0 = 5382.0d0
      fpi.rng  = 0.2d0
      fpi.fr   = fr_crisp
    end
    '5250': begin
      fpi.w0   = 5250.5d0
      fpi.rhr  = 0.9589d0
      fpi.rlr  = 0.8561d0
      fpi.shr  = 787.e4
      fpi.slr  = 295.5e4
      fpi.dhr  = 0.223d0
      fpi.dlr  = 0.609d0
      fpi.prw  = 3.4d0
      fpi.prw0 = 5250.10d0
      fpi.rng  = 0.2d0
      fpi.fr   = fr_crisp
    end
    '5576': begin
      fpi.w0   = 5576.5d0
      fpi.rhr  = 0.9457d0
      fpi.rlr  = 0.8903d0
      fpi.shr  = 787.e4
      fpi.slr  = 295.5e4
      fpi.dhr  = 0.237d0
      fpi.dlr  = 0.647d0
      fpi.prw  = 3.0d0
      fpi.prw0 = 5578.2d0
      fpi.rng  = 0.2d0
      fpi.fr   = fr_crisp
    end
    '5875': begin
      fpi.w0   = 5876.0d0
      fpi.rhr  = 0.9266d0
      fpi.rlr  = 0.8567d0
      fpi.shr  = 787.e4
      fpi.slr  = 295.5e4
      fpi.dhr  = 0.237d0
      fpi.dlr  = 0.647d0
      fpi.prw  = 5.0d0
      fpi.prw0 = 5876.0d0
      fpi.rng  = 0.2d0
      fpi.fr   = fr_crisp
    end
    '6563': begin
      fpi.w0   = 6563.0d0
      fpi.rhr  = 0.9330d0
      fpi.rlr  = 0.8492d0
      fpi.shr  = 787.e4
      fpi.slr  = 295.5e4
      fpi.dhr  = 0.279d0
      fpi.dlr  = 0.761d0
      fpi.prw  = 4.9d0
      fpi.prw0 = 6563.8d0
      fpi.rng  = 0.25d0
      fpi.fr   = fr_crisp
    end
    '6082': begin
      fpi.w0   = 6173d0
      fpi.rhr  = 0.9381d0
      fpi.rlr  = 0.8447d0
      fpi.shr  = 787.e4
      fpi.slr  = 295.5e4
      fpi.dhr  = 0.
      fpi.dlr  = 0.
      fpi.prw  = 0.
      fpi.prw0 = 0.
      fpi.rng  = 0.25d0
      fpi.fr   = fr_crisp
    end
    '6173': begin
      fpi.w0   = 6173.0d0
      fpi.rhr  = 0.9380d0
      fpi.rlr  = 0.8469d0
      fpi.shr  = 787.e4
      fpi.slr  = 295.5e4
      fpi.dhr  = 0.
      fpi.dlr  = 0.
      fpi.prw  = 0.
      fpi.prw0 = 0.
      fpi.rng  = 0.25d0
      fpi.fr   = fr_crisp
    end
    '7090': begin
      fpi.w0   = 7090.0d0
      fpi.rhr  = 0.9430d0
      fpi.rlr  = 0.8565d0
      fpi.shr  = 787.e4
      fpi.slr  = 295.5e4
      fpi.dhr  = 0.
      fpi.dlr  = 0.
      fpi.prw  = 0.
      fpi.prw0 = 0.
      fpi.rng  = 0.25d0
      fpi.fr   = fr_crisp
    end
    '5896': begin
      fpi.w0   = 5896.0d0
      fpi.rhr  = 0.9270d0
      fpi.rlr  = 0.8545d0
      fpi.shr  = 787.e4
      fpi.slr  = 295.5e4
      fpi.dhr  = 0.
      fpi.dlr  = 0.
      fpi.prw  = 0.
      fpi.prw0 = 0.
      fpi.rng  = 0.25d0
      fpi.fr   = fr_crisp
    end
    '5173': begin
      fpi.w0   = 5173.0d0
      fpi.rhr  = 0.9585d0
      fpi.rlr  = 0.8654d0
      fpi.shr  = 787.e4
      fpi.slr  = 295.5e4
      fpi.dhr  = 0.
      fpi.dlr  = 0.
      fpi.prw  = 0.
      fpi.prw0 = 0.
      fpi.rng  = 0.25d0
      fpi.fr   = fr_crisp
    end
    '7772': begin
      fpi.w0   = 7772.0d0
      fpi.rhr  = 0.9227d0
      fpi.rlr  = 0.8576d0
      fpi.shr  = 787.e4
      fpi.slr  = 295.5e4
      fpi.dhr  = 0.
      fpi.dlr  = 0.
      fpi.prw  = 0.
      fpi.prw0 = 0.
      fpi.rng  = 0.50d0
      fpi.fr   = fr_crisp
    end
    ;; CHROMIS:  (Fake parameters below for now.) **************************************** NOTE!!!!
    '3999': begin              
      fpi.w0   = 3999.0d0
      fpi.rhr  = 0.9585d0
      fpi.rlr  = 0.8654d0
      fpi.shr  = 787.e4
      fpi.slr  = 295.5e4
      fpi.dhr  = 0.
      fpi.dlr  = 0.
      fpi.prw  = 0.
      fpi.prw0 = 0.
      fpi.rng  = 0.25d0
      fpi.fr   = fr_chromis
    end
    '3969': begin
      fpi.w0   = 3969.0d0
      fpi.rhr  = 0.9585d0
      fpi.rlr  = 0.8654d0
      fpi.shr  = 787.e4
      fpi.slr  = 295.5e4
      fpi.dhr  = 0.
      fpi.dlr  = 0.
      fpi.prw  = 0.
      fpi.prw0 = 0.
      fpi.rng  = 0.25d0
      fpi.fr   = fr_chromis
    end
    '3978': begin
      fpi.w0   = 3978.0d0
      fpi.rhr  = 0.9585d0
      fpi.rlr  = 0.8654d0
      fpi.shr  = 787.e4
      fpi.slr  = 295.5e4
      fpi.dhr  = 0.
      fpi.dlr  = 0.
      fpi.prw  = 0.
      fpi.prw0 = 0.
      fpi.rng  = 0.25d0
      fpi.fr   = fr_chromis
    end
    '3934': begin
      fpi.w0   = 3934.0d0
      fpi.rhr  = 0.9585d0
      fpi.rlr  = 0.8654d0
      fpi.shr  = 787.e4
      fpi.slr  = 295.5e4
      fpi.dhr  = 0.
      fpi.dlr  = 0.
      fpi.prw  = 0.
      fpi.prw0 = 0.
      fpi.rng  = 0.25d0
      fpi.fr   = fr_chromis
    end
    '3925': begin
      fpi.w0   = 3925.0d0
      fpi.rhr  = 0.9585d0
      fpi.rlr  = 0.8654d0
      fpi.shr  = 787.e4
      fpi.slr  = 295.5e4
      fpi.dhr  = 0.
      fpi.dlr  = 0.
      fpi.prw  = 0.
      fpi.prw0 = 0.
      fpi.rng  = 0.25d0
      fpi.fr   = fr_chromis
    end
    '4862': begin
      fpi.w0   = 4862.0d0
      fpi.rhr  = 0.9585d0
      fpi.rlr  = 0.8654d0
      fpi.shr  = 787.e4
      fpi.slr  = 295.5e4
      fpi.dhr  = 0.
      fpi.dlr  = 0.
      fpi.prw  = 0.
      fpi.prw0 = 0.
      fpi.rng  = 0.25d0
      fpi.fr   = fr_chromis
    end
    else: begin
      print, 'Usage: fpi = red_get_fpi_par(line = <line>)'
      print, '<line> (string): ', lines
    end
  endcase
  
  fpi.fpi_wav = fpi.w0
  
  ;; Center the cavities on fpi.w0
  red_get_fpi_sep, fpi

  return,fpi

end
