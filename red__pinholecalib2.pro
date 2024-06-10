; docformat = 'rst'

;+
;   Compute the (forward and backward) warping matrices to map the
;   camera images to a common reference, which is defined by fitting
;   a regular grid to the selected 'reference'
;
; :Author:
; 
;     Pit Sütterlin, Institute for Solar Physics, 2024
;
; :History:
;
;   2024-03-25 : Switch warping mechanism to polynomial matrices
;   2024-04-05 : PS normalize PHs to max intensity for area computation (clouds!)
;-

PRO red::pinholecalib2, avg_fpistates = avg_fpistates $
                        , cams=cams $
                        , dir=dir $
                        , manual=manual $
                        , margin=margin $
                        , overwrite=overwrite $
                        , pref=pref $
                        , refcam=refcam $
                        , threshold=threshold $
                        , verbose=verbose $
                        , verify=verify $
                        , warp_dim=warp_dim


  ;; Name of this method
inam = red_subprogram(/low, calling = inam1)  
  
  ;; Logging
help, /obj, self, output = selfinfo 
red_writelog, selfinfo = selfinfo, logfile = logfile

IF (n_elements(threshold) eq 0) then threshold = 0.0
IF ( n_elements(dir) gt 0 ) then dir = [dir] $
ELSE IF  ptr_valid(self.pinh_dirs) THEN dir = *self.pinh_dirs
self -> getdetectors, dir = dir
IF ( n_elements(cams) GT 0 ) THEN $
  cams = [cams] $
ELSE IF ptr_valid(self.cameras) THEN $
  cams = [*self.cameras]

IF n_elements(refcam) EQ 0 THEN refcam = self.refcam
IF n_elements(verbose) EQ 0 THEN verbose = 0
IF n_elements(warp_dim) EQ 0 THEN warp_dim = 3
IF n_elements(max_shift) EQ 0 THEN max_shift = 200
IF n_elements(margin) EQ 0 THEN margin = 30
  
Ncams = n_elements(cams)
  
    ;; Search summed pinh images and camtag
ph_dir = self.out_dir+'/pinhs/'
ph_dir = red_strreplace(ph_dir, '//', '/')
  
output_dir = self.out_dir+'/calib/'
output_dir = red_strreplace(output_dir, '//', '/')
output_file = output_dir+'alignments_new.sav'
file_mkdir, output_dir

IF Ncams LT 2 THEN BEGIN
    print, inam, ' : Need at least 2 cameras.'
    print, inam, ' : dir: ',dir
    print, inam, ' : cams: ',cams
    print, inam, ' : self.cameras: ',*self.cameras
    red_writelog, /add, logfile = logfile, $
                  top_info_strings = ' : Need at least 2 cameras.'
    return
ENDIF
  
IF refcam GE Ncams THEN BEGIN
    print, inam, ' : index of reference camera out of range: ', $
           refcam, ' >= ', Ncams
    return
ENDIF
  
IF file_test(output_file) THEN IF NOT keyword_set(overwrite) THEN BEGIN
    print, inam, ' : file exists: ', output_file
    print, inam, ' : Use /OVERWTRITE to create new'
    return
ENDIF

all_files = file_search(ph_dir + '*.pinh.fits', count = nf)

   ;;; states for the reference camera
   ;;; NB: This internally takes care of the CHROMIS WB/NB combinations
self->selectfiles, files = all_files, states = states, cam = cams[refcam] $
                   , /strip_settings, selected = selection, prefilter = pref
ref_states = states[selection]
ref_states_unique = ref_states.fpi_state
  ;; We discard multiple pinhole files for the same state
ref_states_unique = ref_states[uniq(ref_states_unique, sort(ref_states_unique))]
n_ref = n_elements(ref_states_unique)
  ;;; Loop through ref states, find matches based on fpi_state, and process those

undefine, alignments
openw, unit, 'calib/refgrid.txt', /get, /append
printf, unit, '#   State         Center         Step   Angle'

FOR i_ref = 0, n_ref-1 DO BEGIN
    l_shape = [7, 3]
    idx = where((states.fpi_state EQ ref_states_unique[i_ref].fpi_state) AND $
                (states.camera NE ref_states_unique[i_ref].camera), n_idx)
    IF n_idx EQ 0 THEN CONTINUE
      ;;; process the ref image, i.e., fit the grid
    this_ref_state = ref_states_unique[i_ref]
    ref_img = red_readdata(this_ref_state.filename, /silent)
    ref_siz = size(ref_img, /dim)
      ;;; find all pinholes, locate the 'L', and fit a regular grid
    ref_m = red_separate_mask(ref_img GE max(ref_img)/10.)
    np_ref = max(ref_m)
      ;;; compute area and positions
    ref_area = lonarr(np_ref)
    ref_pos = fltarr(2, np_ref)
    FOR i = 0, np_ref-1 DO BEGIN
        ;;;ix = where_n(ref_m EQ (i+1))
        ix = where(ref_m EQ (i+1))
        ix = [[ix MOD ref_siz[0]], [ix / ref_siz[0]]]
        lu = (min(ix, dim = 1, max = ro)-3) > [0, 0]
        ro = (ro+3) < (ref_siz-1)
        bx = ref_img[lu[0]:ro[0], lu[1]:ro[1]]
        ref_pos[*, i] = lu + red_centroid(bx)
        ref_area[i] = total(bx GT max(bx)/5.)
    ENDFOR
    IF keyword_set(manual) THEN BEGIN
        print, 'Mark 3 pinholes in shape of an L'
        ref_lind = red_markpinholes(ref_m, ref_pos, title='Reference Pinholes', /keep)
        read, 'Enter length of long and short arm in gridsteps: ', l_shape
        ref_lpos = ref_pos[*, ref_lind]
        IF red_lcheck(ref_lpos, ref_lind, gridstep, ratio=l_shape) NE 1 THEN stop
    ENDIF ELSE BEGIN
          ;;; try to locate the 7/3 in the 6 largest clusters
        ix = reverse(sort(ref_area))
        FOR i=0, 3 DO FOR j=i+1, 4 DO FOR k=j+1, 5 DO BEGIN
            ref_lind = ix[[i, j, k]]
            ref_lpos = ref_pos[*, ref_lind] 
            IF red_lcheck(ref_lpos, ref_lind, gridstep, ratio=l_shape) EQ 1 THEN GOTO, found_l1
        ENDFOR
        print, inam, ' : Unable to locate the reference pinholes! Mark them manually'
        REPEAT BEGIN
            ref_lind = red_markpinholes(ref_m, ref_pos, title='Reference Pinholes')
            ref_lpos = ref_pos[*, ref_lind]
            IF red_lcheck(ref_lpos, ref_lind, gridstep, ratio=l_shape) EQ 1 THEN GOTO, found_l1
            print, inam, ' : This is not the 3 bigger pinholes!'
        ENDREP UNTIL 0
    ENDELSE
found_l1:
    print, inam, ' : Found the reference L'
      ;;; lcheck also re-orders positions to [corner, edge_l, edge_s]
      ;;; corner position == center
    X0 = ref_pos[*, ref_lind[0]]
      ;;; build the true reference grid by fitting, relative to X0
    tmp = ref_pos - X0 # replicate(1, np_ref)
    tmp_i = round(tmp/gridstep)  ;;; index
    ang = atan(ref_lpos[1, 1]-ref_lpos[1, 0], ref_lpos[0, 1]-ref_lpos[0, 0])
    IF ref_lpos[0, 1] LT ref_lpos[0, 0] THEN ang += !pi
    ;print, gridstep, ang
    par = mpfit('red_grid_func', [gridstep, ang], $
                functargs = {X:tmp_i, Y:tmp}, /quiet)
    print, inam, ' : regular grid width ', strtrim(par[0], 2), 'px, tilt ', strtrim(par[1]*!radeg, 2), 'deg'
       ;;; check that there are no double index positions
    iix = replicate(1b, np_ref)
    FOR i=0, np_ref-2 DO FOR j=i+1, np_ref-1 DO $
      IF total(tmp_i[*, i] EQ tmp_i[*, j]) EQ 2 THEN BEGIN
        iix[i] = 0
        ref_pos[*, j] = (ref_pos[*, i]+ref_pos[*, j])/2
    ENDIF
    iix = where(iix EQ 1)
    ref_pos = ref_pos[*, iix]
    reg_i = tmp_i[*, iix]
    np_ref = (size(reg_i, /dim))[1]
    reg_pos = x0#replicate(1, np_ref)-red_grid_func(par, x = reg_i, y = 0)

    printf, unit, string(this_ref_state.fpi_state, X0, par, $
                         form="(a-12,f7.2,'  ',f7.2,'   ',f7.3,f+7.3)")
    
      ;;; compute warp matrices, forward (kxy) and backward(lxy)
    map = transpose([ref_pos, reg_pos])
    red_clean_warp, map, kx, ky, warp_dim
    red_clean_warp, map[*, [2, 3, 0, 1]], lx, ly, warp_dim
    
    red_append, alignments, { state:this_ref_state, map_x: kx, map_y: ky, $
                              revmap_x: lx, revmap_y: ly }

      ;;; Now loop through dependent files and align to grid
    
    FOR i_dep = 0, n_idx-1 DO BEGIN
        state = states[idx[i_dep]]
        img = red_readdata(state.filename, /silent)
        mask = img GE max(img)/10.
        img_m = red_separate_mask(mask)
        img_np = max(img_m)
        img_area = lonarr(img_np)
        img_pos = fltarr(2, img_np)
        FOR i = 0, img_np-1 DO BEGIN
            ;;;ix = where_n(img_m EQ i+1)
            ix = where(img_m EQ (i+1))
            ix = [[ix MOD ref_siz[0]], [ix / ref_siz[0]]]
            lu = (min(ix, dim = 1, max = ro)-3) > [0, 0]
            ro = (ro+3) < (ref_siz-1)
            bx = img[lu[0]:ro[0], lu[1]:ro[1]]
            img_pos[*, i] = lu + red_centroid(bx)
            img_area[i] = total(bx GT max(bx)/5.)
        ENDFOR
        IF keyword_set(manual) THEN BEGIN
            print, 'Mark the same three pinholes as for the reference'
            img_lind = red_markpinholes(img_m, img_pos)
            img_lpos = img_pos[*, img_lind]
            IF red_lcheck(img_lpos, img_lind, gridstep, ratio=l_shape) NE 1 THEN stop
        ENDIF ELSE BEGIN
            ix = reverse(sort(img_area))
            FOR i=0, 3 DO FOR j=i+1, 4 DO FOR k=j+1, 5 DO BEGIN
                img_lind = ix[[i, j, k]]
                img_lpos = img_pos[*, img_lind] 
                IF red_lcheck(img_lpos, img_lind, gridstep, ratio=l_shape) EQ 1 THEN GOTO, found_l2
            ENDFOR
            print, inam, ' : Unable to locate the reference pinholes. Try marking them manually'
            REPEAT BEGIN
                img_lind = red_markpinholes(img_m, img_pos, title='Dependent Pinholes')
                img_lpos = img_pos[*, img_lind]
                IF red_lcheck(ref_lpos, ref_lind, gridstep, ratio=l_shape) EQ 1 THEN GOTO, found_l2
                print, inam, ' : This is not the 3 bigger pinholes!'
            ENDREP UNTIL 0
        ENDELSE

found_l2:
        XX0 = img_lpos[*, 0]
          ;;; find the matching pinholes
        tmp = img_pos-xx0#replicate(1, img_np)
        tmp_i = round(tmp/gridstep)
          ;;; this assumes large L side is horizontal.  Check/improve
        hflip = reg_i[0, ref_lind[1]] / tmp_i[0, img_lind[1]]
        vflip = reg_i[1, ref_lind[2]] / tmp_i[1, img_lind[2]]
        n_pairs = 0
        undefine, nmap
        FOR ph = 0, img_np-1 DO BEGIN
            ix = where( (reg_i[0, *] EQ hflip*tmp_i[0, ph]) AND $
                        reg_i[1, *] EQ vflip*tmp_i[1, ph], nf)
            IF nf EQ 0 THEN CONTINUE
              ;;; have a match, add the pair to map
            red_append, nmap, [tmp[*, ph]+xx0, reg_pos[*, ix]] 
            n_pairs++
        ENDFOR
        nmap = transpose(reform(nmap, 4, n_pairs))
        red_clean_warp, nmap, kx, ky, warp_dim
        red_clean_warp, nmap[*, [2, 3, 0, 1]], lx, ly, warp_dim
        
        red_append, alignments, { state:state, map_x: kx, map_y: ky, $
                                  revmap_x: lx, revmap_y: ly }
    ENDFOR    ;;; loop dependent files
    IF keyword_set(manual) THEN wdelete
ENDFOR        ;;; loop reference files

free_lun, unit

save, file = output_file, alignments

END
