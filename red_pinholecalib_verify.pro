; docformat = 'rst'

;+
; Verify pinhole calibration by comparing the reference (WB) image vs
; the other images/cameras after applying the calibrations.
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
; :Keywords:
; 
;   diff : in, optional, type=boolean
;   
;     Display difference image instead of blinking.
; 
;   method : in, optional, type=string
; 
;     Set to "projective" if both kinds of calibrations exist and you
;     want to verify that with the older method.
; 
;   pref : in, optional, type=string
; 
;     Verify calibration only for this prefilter.
; 
; :History:
; 
;    2025-11-02 : MGL. First version.
; 
;-
pro red_pinholecalib_verify, diff = diff, method = method, pref = pref

  if n_elements(method) eq 0 then begin
    ;; Method not specified
    case !true of
      file_test('calib/alignments_polywarp.sav') : method = 'polywarp'
      file_test('calib/alignments.sav')          : method = 'projective'
      else : begin
        red_message, 'No pinholecalib output in calib/'
        return
      end
    endcase
  endif

  case method of
    'polywarp'   : afile = 'calib/alignments_polywarp.sav'
    'projective' : afile = 'calib/alignments.sav'
    else : begin
      red_message, 'Unknown method: '+method
      return
    end
  endcase

  if ~file_test(afile) then begin
    red_message, 'File does not exist: '+afile
    return
  endif

  restore, afile
  
  case method of
    'polywarp'   : begin
      if n_elements(pref) gt 0 then begin
        indx = where(alignments.state.prefilter eq pref, Nwhere)
        if Nwhere eq 0 then return
        alignments = alignments[indx]
      endif
      fullstates = red_uniquify(alignments.state.fullstate, count = Nstates)
    end 
    'projective' : begin
      if n_elements(pref) gt 0 then begin
        indx = where(alignments.state1.prefilter eq pref or $
                     alignments.state2.prefilter eq pref, Nwhere)
        if Nwhere eq 0 then return
        alignments = alignments[indx]
      endif
      fullstates = red_uniquify(alignments.state2.fullstate, count = Nstates)
    end
  endcase
 
  print, fullstates
  
  for istate = 0, Nstates-1 do begin

    case method of
      'polywarp'   : begin
        indx = where(alignments.state.fullstate eq fullstates[istate])
        cams = alignments[indx].state.camera
        files = alignments[indx].state.filename
      end
      'projective' : begin
        indx = where(alignments.state2.fullstate eq fullstates[istate])
        cams = [alignments[indx[0]].state1.camera, alignments[indx].state2.camera]
        files = [alignments[indx[0]].state1.filename, alignments[indx].state2.filename]
      end
    endcase

    ims = red_readdata_multiframe(files, /silent)
    Ncams = n_elements(cams)
    ref_indx = where(strmatch(cams,'*-W'), complement=oth_indx, Ncomplement = Nother)
    
    dims = size(ims, /dim)
    aspect = dims[0]/float(dims[1])
    device, get_scr = sz
    if dims[0]*Nother gt sz[0] then xs = 10*(sz[0]/(10*Nother)) else xs = dims[0]
    ys = xs / aspect
    window, xs = Nother*xs, ys = ys, title = fullstates[istate]

    for icam = 0, Ncams-1 do begin
      case method of
        'polywarp'   : amap = {x:alignments[indx[icam]].map_x, $
                               y:alignments[indx[icam]].map_y }
        'projective' : if icam eq 0 then amap = diag_matrix([1.,1.,1.]) else amap = invert(alignments[indx[icam-1]].map)
      endcase
      ims[*, *, icam] = red_apply_camera_alignment(ims[*, *, icam], method $
                                                   , /preserve $
                                                   , amap = amap)
      ims[*, *, icam] /= max(ims[*, *, icam])
    endfor

    print, fullstates[istate]
    if ~keyword_set(diff) then begin
      blink_i = !true
      print, 'Hit RETURN to continue'
      while (get_kbrd(0) eq '') do begin
        if blink_i then begin
          for icam = 0, Nother-1 do tvscl, congrid(ims[*, *, ref_indx[0]], xs, ys), icam
          blink_i = !false        
        endif else begin
          for icam = 0, Nother-1 do tvscl, congrid(ims[*, *, oth_indx[icam]], xs, ys), icam
          blink_i = !true
        endelse
        wait, .5
      endwhile
    endif else begin
      for icam = 0, Nother-1 do tvscl, congrid(ims[*, *, ref_indx[0]] - ims[*, *, oth_indx[icam]], xs, ys), icam
      if istate lt Nstates-1 then begin
        s = ''
        read, 'Hit RETURN to continue', s
      endif
    endelse
    
  endfor                        ; istate
  
end

cd, '/scratch_beegfs/mats/NEW/2024-04-21/CHROMIS-test/'
cd, '/scratch_beegfs/mats/NEW/2025-08-18/CRISP-testluc/'
red_pinholecalib_verify

end
