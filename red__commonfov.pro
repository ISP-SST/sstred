; docformat = 'rst'

;+
;   Calculate the common FOV for an arbitrary
;   subset of cameras and states.
;
; 
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
;     2024-03-17 : First version, based on old openCV routine from Tomas
;     2024-06-04 : Pit - Replace red_warp_coord with the rdx C routine
;-

FUNCTION red::commonfov, align = align, $
                         cams = cams, $
                         extraclip = extraclip, $
                         output_dir = output_dir, $
                         prefilters = prefilters

  ;; Name of this method
inam = red_subprogram(/low, calling = inam1)              

IF n_elements(cams) GT 0 THEN $
  cams = [cams] $
ELSE $
  IF ptr_valid(self.cameras) THEN cams = [*self.cameras]

IF n_elements(verbose) EQ 0 THEN verbose = 0

IF n_elements(output_dir) EQ 0 THEN output_dir = self.out_dir+'/calib/'
alignfile = self.out_dir+'/calib/alignments_new.sav'

IF file_test(alignfile) THEN BEGIN
    restore, alignfile
ENDIF ELSE BEGIN
    print, inam + ' : Alignment file missing.'
    print, inam + ' : Did you run a -> pinholecalib ?'
    return, 0
ENDELSE

CASE n_elements(extraclip) OF
    0 : extraclip = [0L, 0L, 0L, 0L]
    1 : extraclip = replicate(extraclip, 4)
    2 : extraclip = [ replicate(extraclip[0], 2), $
                      replicate(extraclip[1], 2) ]
    4 :                         ; Leave as it is.
    ELSE : BEGIN
        print, inam + "ERROR: Don't know how to use keyword extraclip with " $
               + strtrim(n_elements(extraclip), 2) + ' elements.'
        stop
    END
ENDCASE

  ;; filter selection
IF (Npref = n_elements(prefilters)) GT 0 THEN BEGIN
    undefine, selected
    FOR i = 0, Npref-1 DO $
      red_append, selected, where(alignments.state.prefilter EQ prefilters[i])
    IF min(selected) LT 0 THEN BEGIN
        print, inam, ' : Not all requested prefilters could be found!'
        return, 0
    ENDIF ELSE $
      alignments = alignments[selected]
ENDIF

  ;;; camera selection
Ncams = n_elements(cams)
IF Ncams LT 2 THEN BEGIN
    print, inam, ' : Need at least 2 cameras.'
    print, inam, ' : cams: ', cams
    print, inam, ' : self.cameras: ', *self.cameras
    return, 0
ENDIF

undefine, selected
FOR i = 0, Ncams-1 DO $
  red_append, selected, where(alignments.state.camera EQ cams[i])
IF min(selected) LT 0 THEN BEGIN
    print, inam, ' : Not all requested prefilters could be found!'
    return, 0
ENDIF ELSE $
  alignments = alignments[selected]

Nalign = n_elements(alignments)

  ;;; image dimensons.  All are identical (checked by pinholecalib)
head = red_readhead(alignments[0].state.filename, /silent)
dims = fxpar(head, 'NAXIS*')

ref_corners = [[0, 0], [dims-1]]
common_fov = ref_corners

FOR ialign = 0, Nalign-1 DO BEGIN
      ;;; project to respective camera plane and clip
    common_fov = rdx_point_warp(alignments[ialign].map_x, alignments[ialign].map_y, common_fov)
    common_fov[*, 0] = common_fov[*, 0] > 0 < (dims[0]-1)
    common_fov[*, 1] = common_fov[*, 1] > 0 < (dims[1]-1)
      ;;; project back, clip again
    common_fov = rdx_point_warp(alignments[ialign].revmap_x, alignments[ialign].revmap_y, common_fov)
    common_fov[*, 0] = common_fov[*, 0] > 0 < (dims[0]-1)
    common_fov[*, 1] = common_fov[*, 1] > 0 < (dims[1]-1)
ENDFOR
  ;;; round inwards and add extraclip
common_fov[0, 0] = fix(common_fov[0, 0]) + extraclip[0] + 1
common_fov[0, 1] = fix(common_fov[0, 1]) - extraclip[1]
common_fov[1, 0] = fix(common_fov[1, 0]) + extraclip[2] + 1
common_fov[1, 1] = fix(common_fov[1, 1]) - extraclip[3]
common_fov = fix(common_fov)
  ;;; clip needs to be [x0,x1,y0,y1]
common_fov = transpose(common_fov)

align = alignments

return, common_fov

END
