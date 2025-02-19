; docformat = 'rst'

;+
;   Calculate the common FOV for an arbitrary
;   subset of cameras and states.
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;     Pit Sütterlin, Institute for Solar Physics, 2024
;
; :History:
;
;     2024-03-17 : PS. First version, based on old openCV routine from Tomas
;
;     2024-06-04 : PS. Replace red_warp_coord with the rdx C routine
;
;-
FUNCTION red::commonfov, align = align, $
                         cams = cams, $
                         extraclip = extraclip, $
                         output_dir = output_dir, $
                         prefilters = prefilters

  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)              
  
  if n_elements(cams) gt 0 then $
     cams = [cams] $
  else $
     if ptr_valid(self.cameras) then cams = [*self.cameras]
  
  IF n_elements(verbose) EQ 0 THEN verbose = 0
  
  if n_elements(output_dir) eq 0 then output_dir = self.out_dir+'/calib/'
  alignfile = self.out_dir+'/calib/alignments_polywarp.sav'
  
  if file_test(alignfile) then begin
    restore, alignfile
  endif else begin
    print, inam + ' : Alignment file missing.'
    print, inam + ' : Did you run a -> pinholecalib ?'
    return, 0
  endelse
  
  case n_elements(extraclip) of
    0 : extraclip = [0L, 0L, 0L, 0L]
    1 : extraclip = replicate(extraclip, 4)
    2 : extraclip = [ replicate(extraclip[0], 2), $
                      replicate(extraclip[1], 2) ]
    4 :                         ; Leave as it is.
    else : begin
      print, inam + "ERROR: Don't know how to use keyword extraclip with " $
             + strtrim(n_elements(extraclip), 2) + ' elements.'
      stop
    end
  endcase
  
  ;; filter selection
  if (Npref = n_elements(prefilters)) gt 0 then begin
    undefine, selected
    for i = 0, Npref-1 do $
       red_append, selected, where(alignments.state.prefilter EQ prefilters[i])
    if min(selected) lt 0 then begin
      print, inam, ' : Not all requested prefilters could be found!'
      return, 0
    endif else $
       alignments = alignments[selected]
  endif
  
  ;;; camera selection
  Ncams = n_elements(cams)
  if Ncams lt 2 then begin
    print, inam, ' : Need at least 2 cameras.'
    print, inam, ' : cams: ', cams
    print, inam, ' : self.cameras: ', *self.cameras
    return, 0
  endif
  
  undefine, selected
  for i = 0, Ncams-1 do $
     red_append, selected, where(alignments.state.camera EQ cams[i])
  if min(selected) lt 0 then begin
    print, inam, ' : Not all requested prefilters could be found!'
    return, 0
  endif else $
     alignments = alignments[selected]
  
  Nalign = n_elements(alignments)
  
  ;;; image dimensons.  All are identical (checked by pinholecalib)
  head = red_readhead(alignments[0].state.filename, /silent)
  dims = fxpar(head, 'NAXIS*')
  
  ref_corners = [[0, 0], [dims-1]]
  common_fov = ref_corners
  
  for ialign = 0, Nalign-1 do begin
    ;; project to respective camera plane and clip
    common_fov = rdx_point_warp(alignments[ialign].map_x, alignments[ialign].map_y, common_fov)
    common_fov[*, 0] = common_fov[*, 0] > 0 < (dims[0]-1)
    common_fov[*, 1] = common_fov[*, 1] > 0 < (dims[1]-1)
      ;;; project back, clip again
    common_fov = rdx_point_warp(alignments[ialign].revmap_x, alignments[ialign].revmap_y, common_fov)
    common_fov[*, 0] = common_fov[*, 0] > 0 < (dims[0]-1)
    common_fov[*, 1] = common_fov[*, 1] > 0 < (dims[1]-1)
  endfor                        ; ialign
  ;; round inwards and add extraclip
  common_fov[0, 0] = fix(common_fov[0, 0]) + extraclip[0] + 1
  common_fov[0, 1] = fix(common_fov[0, 1]) - extraclip[1]
  common_fov[1, 0] = fix(common_fov[1, 0]) + extraclip[2] + 1
  common_fov[1, 1] = fix(common_fov[1, 1]) - extraclip[3]
  common_fov = fix(common_fov)
  ;; clip needs to be [x0,x1,y0,y1]
  common_fov = transpose(common_fov)
  
  align = alignments
  
  return, common_fov
  
end
