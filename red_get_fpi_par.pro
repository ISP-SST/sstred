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
; :author:
; 
;    Jaime de la Cruz Rodriguez, ISP-KVA, 2010
; 
; 
; :returns:
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
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2013-07-11 : MGL. Use red_get_fpi_sep, not get_fpi_sep.
; 
; 
;-
function red_get_fpi_par,line=line

  pi2 = 6.28318530779d0    ; 2.d0 * PI
  calp = cos(sqrt([0.1127017, 0.5000000, 0.8872983] * (0.5/165.)^2.)) * pi2
  wng = [0.2777778, 0.4444444, 0.2777778]
  ;
  fpi=create_struct('rhr',0.d0,'rlr',0.d0,'slr',0.d0,'shr',0.d0,$
                    'dhr',0.d0,'dlr',0.d0,'w0',0.d0,'fr',165.0d0,'thmax',0.5d0/165.d0,$
                    'prw',0.d0,'prw0',0.d0,'prnc',2.e0,'fpi_wav',0.d0,$
                    'rng', 0.0d0, 'calp', calp, 'wng', wng, 'pi2', pi2)
;
  lines = ['5173','5380','5250','5576', '5896','6082','6173','6302','6563','7090','7772','8542']
;
  if not keyword_set(line) then begin
     print, 'Usage: fpi = red_get_fpi_par(line = <line>)'
     print, '<line> (string): ',lines
     return,0
  endif
; 
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
     endcase
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
     endcase
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
     endcase
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
     endcase
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
     endcase
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
     endcase
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
     endcase
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
     endcase
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
     endcase
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
     endcase
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
     endcase
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
     endcase
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
     endcase
     else: begin
        print, 'Usage: fpi = red_get_fpi_par(line = <line>)'
        print, '<line> (string): ', lines
     endcase
  endcase
  fpi.fpi_wav=fpi.w0
  ;
  ; center the cavities on fpi.w0
  ;
  red_get_fpi_sep,fpi
  ;
  return,fpi
end
