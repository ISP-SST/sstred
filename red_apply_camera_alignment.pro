; docformat = 'rst'

;+
; Read the appropriate camera alignment file and apply the warping to
; an image.
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
;   im : in, type=image
; 
;     The image to be warped to the reference geometry.
;
;   cam : in, type=string
;
;     The camera from which the image comes.
; 
;   model : in, type=string
; 
;     The geometrical transform model, one of "projective" (as
;     described in the SSTRED paper) and "polywarp".
; 
; :Keywords:
; 
;   amap : in, out, optional, type="matrix or struct"
; 
;     The transform mapping as a matrix in the projective case and a
;     struct in the polywarp case. If not defined, will read from the
;     calib/alignments*.sav file and be defined on exit. If defined,
;     will be used without checking if it's appropriate for the
;     camera.
; 
;   pref : in, optional, type=string
;
;     The relevant prefilter.
;
;   preserve_size : in, optional, type=boolean
;
;      Used for the projective transform.
;
; :History:
; 
;   2025-02-18 : MGL. First version.
; 
;-
function red_apply_camera_alignment, im, model, cam $
                                     , pref = pref $
                                     , preserve_size = preserve_size $
                                     , amap = amap
 
  case model of 

    'polywarp' : begin

      if n_elements(amap) eq 0 then begin
        
        if file_test('calib/alignments_polywarp.sav') then begin
        
          ;; We have alignment data from the new polywarp pinhole
          ;; alignment procedure.
          
          restore, 'calib/alignments_polywarp.sav'
          if n_elements(pref) gt 0 then begin
            indx = where(alignments.state.camera eq cam $
                         and alignments.state.prefilter eq pref, Nalign)
          endif else begin
            indx = where(alignments.state.camera eq cam, Nalign)
          endelse
          
          case Nalign of
            0    : stop         ; Should not happen!
            1    : amap = {X: alignments[indx].map_x, $
                           Y: alignments[indx].map_y}
            else : amap = {X: mean(alignments[indx].map_x, dim=3), $
                           Y: mean(alignments[indx].map_y, dim=3) }
          endcase

        endif else begin

          red_message, 'Cannot find  calib/alignments_polywarp.sav. Did you run pinholecalib?'
          retall

        endelse

      endif else begin
        if ~isa(amap, 'STRUCT') then begin
          red_message, 'Polywarp transform specified but given amap is not a struct.'
          help, amap
          stop
        endif
      endelse
      
      return, poly_2d(im, amap.x, amap.y, 1, miss=0)
      
    end

    'projective' : begin

      if n_elements(amap) eq 0 then begin
        
        if file_test('calib/alignments.sav') then begin
      
          ;; We have data from the old projective transform matrix
          ;; pinhole calibration.
          
          restore, 'calib/alignments.sav'
          
          indx = where(alignments.state2.camera eq 'Crisp-T', Nalign)
          case Nalign of
            0    : stop         ; Should not happen!
            1    : amap = invert(      alignments[indx].map           )
            else : amap = invert( mean(alignments[indx].map, dim = 3) )
          endcase
          amap /= amap[2, 2]    
        endif else begin

          red_message, 'Cannot find  calib/alignments.sav. Did you run pinholecalib?'
          retall

        endelse

      endif else begin
        if ~isa(amap, /array, /number) then begin
          red_message, 'Projective transform specified but given amap is not an array.'
          help, amap
          stop
        endif
      endelse
      
      return, rdx_img_project(amap, im, preserve_size = preserve_size)
      
    end

    else : begin

      if keyword_set(use_old) then begin
        red_message, 'Cannot find  calib/alignments.sav. Did you run pinholecalib?'
      endif else begin
        red_message, 'Cannot find any alignment data in calib/. Did you run pinholecalib?'
      endelse

      retall
        
    endelse
    
  endcase

end




