; docformat = 'rst'

;+
;   Compute the (forward and backward) warping matrices to map the
;   camera images to a common reference, which is defined by fitting
;   a regular grid to the selected 'reference'.
;
; :Categories:
;
;    SST pipeline
; 
; :Author:
; 
;     Pit Sütterlin, Institute for Solar Physics, 2024
;
; :Returns:
; 
; 
; :Params:
; 
; 
; 
; 
; 
; 
; :Keywords:
; 
;   
;   
;   
; 
; 
; :History:
; 
;   2024-03-25 : PS. Switch warping mechanism to polynomial matrices.
;
;   2024-04-05 : PS. Normalize PHs to max intensity for area
;                computation (clouds!)
;
;-
pro red::pinholecalib2, avg_fpistates = avg_fpistates $ ; Keyword not used
                        , cams=cams $
                        , dir=dir $
                        , manual=manual $
                        , margin=margin $
                        , overwrite=overwrite $
                        , pref=pref $
                        , refcam=refcam $
                        , threshold=threshold $
                        , verbose=verbose $ ; Keyword not used
                        , verify=verify $   ; Keyword not used
                        , warp_dim=warp_dim

  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)  
  
  ;; Logging
  help, /obj, self, output = selfinfo 
  red_writelog, selfinfo = selfinfo, logfile = logfile

  if (n_elements(threshold) eq 0) then threshold = 0.0
  if ( n_elements(dir) gt 0 ) then dir = [dir] $
  else if  ptr_valid(self.pinh_dirs) then dir = *self.pinh_dirs
  self -> getdetectors, dir = dir
  if ( n_elements(cams) gt 0 ) then $
     cams = [cams] $
  else if ptr_valid(self.cameras) then $
     cams = [*self.cameras]
  
  if n_elements(refcam) eq 0 then refcam = self.refcam
  if n_elements(verbose) eq 0 then verbose = 0
  if n_elements(warp_dim) eq 0 then warp_dim = 3
  if n_elements(max_shift) eq 0 then max_shift = 200
  if n_elements(margin) eq 0 then margin = 30
  
  Ncams = n_elements(cams)
  
  ;; Search summed pinh images and camtag
  ph_dir = self.out_dir+'/pinhs/'
  ph_dir = red_strreplace(ph_dir, '//', '/')
  
  output_dir = self.out_dir+'/calib/'
  output_dir = red_strreplace(output_dir, '//', '/')
  output_file = output_dir+'alignments_polywarp.sav'
  file_mkdir, output_dir

  if Ncams lt 2 then begin
    print, inam, ' : Need at least 2 cameras.'
    print, inam, ' : dir: ',dir
    print, inam, ' : cams: ',cams
    print, inam, ' : self.cameras: ',*self.cameras
    red_writelog, /add, logfile = logfile, $
                  top_info_strings = ' : Need at least 2 cameras.'
    return
  endif
  
  if refcam ge Ncams then begin
    print, inam, ' : index of reference camera out of range: ', $
           refcam, ' >= ', Ncams
    return
  endif
  
  if file_test(output_file) then if not keyword_set(overwrite) then begin
    print, inam, ' : file exists: ', output_file
    print, inam, ' : Use /OVERWTRITE to create new'
    return
  endif
  
  all_files = file_search(ph_dir + '*.pinh.fits', count = nf)
  
  ;; states for the reference camera
  ;; NB: This internally takes care of the CHROMIS WB/NB combinations
  self->selectfiles, files = all_files, states = states, cam = cams[refcam] $
                             , /strip_settings, selected = selection, prefilter = pref
  ref_states = states[selection]
  ref_states_unique = ref_states.fpi_state
  ;; We discard multiple pinhole files for the same state
  ref_states_unique = ref_states[uniq(ref_states_unique, sort(ref_states_unique))]
  n_ref = n_elements(ref_states_unique)
  ;; Loop through ref states, find matches based on fpi_state, and process those
  
  undefine, alignments
  openw, unit, 'calib/refgrid.txt', /get, /append
  printf, unit, '#   State         Center         Step   Angle'
  
  for i_ref = 0, n_ref-1 do begin
    l_shape = [7, 3]
    idx = where((states.fpi_state eq ref_states_unique[i_ref].fpi_state) and $
                (states.camera NE ref_states_unique[i_ref].camera), n_idx)
    if n_idx eq 0 then continue
    ;; process the ref image, i.e., fit the grid
    this_ref_state = ref_states_unique[i_ref]
    ref_img = red_readdata(this_ref_state.filename, /silent)
    ref_siz = size(ref_img, /dim)
    ;;; For CRISP wideband mask the stronger distorted corners
    if this_ref_state.camera eq 'Crisp-W' then ref_img *= (shift(dist(ref_siz), ref_siz/2) le 0.49*ref_siz[0])
    ;; find all pinholes, locate the 'L', and fit a regular grid
    ref_m = red_separate_mask(ref_img GE max(ref_img)/10.)
    np_ref = max(ref_m)
    ;; compute area and positions
    ref_area = lonarr(np_ref)
    ref_pos = fltarr(2, np_ref)
    for i = 0, np_ref-1 do begin
      ;;ix = where_n(ref_m EQ (i+1))
      ix = where(ref_m eq (i+1))
      ix = [[ix mod ref_siz[0]], [ix / ref_siz[0]]]
      lu = (min(ix, dim = 1, max = ro)-3) > [0, 0]
      ro = (ro+3) < (ref_siz-1)
      bx = ref_img[lu[0]:ro[0], lu[1]:ro[1]]
      ref_pos[*, i] = lu + red_centroid(bx)
      ref_area[i] = total(bx gt max(bx)/5.)
    endfor
    if keyword_set(manual) then begin
      print, 'Mark 3 pinholes in shape of an L'
      ref_lind = red_markpinholes(ref_img, ref_pos, title='Reference Pinholes', /keep)
      read, 'Enter length of long and short arm in gridsteps: ', l_shape
      ref_lpos = ref_pos[*, ref_lind]
      IF red_lcheck(ref_lpos, ref_lind, gridstep, ratio=l_shape) NE 1 THEN stop
    endif else begin
      ;; try to locate the 7/3 in the 6 largest clusters
      ix = reverse(sort(ref_area))
      for i=0, 3 do for j=i+1, 4 do for k=j+1, 5 do begin
        ref_lind = ix[[i, j, k]]
        ref_lpos = ref_pos[*, ref_lind] 
        if red_lcheck(ref_lpos, ref_lind, gridstep, ratio=l_shape) eq 1 then goto, found_l1
      endfor
      print, inam, ' : Unable to locate the reference pinholes! Mark them manually'
      repeat begin
        ref_lind = red_markpinholes(ref_img, ref_pos, title='Reference Pinholes')
        ref_lpos = ref_pos[*, ref_lind]
        if red_lcheck(ref_lpos, ref_lind, gridstep, ratio=l_shape) eq 1 then goto, found_l1
        print, inam, ' : This is not the 3 bigger pinholes!'
      endrep until 0
    endelse

found_l1:

    print, inam, ' : Found the reference L for camera ', this_ref_state.camera
    ;; lcheck also re-orders positions to [corner, edge_l, edge_s]
    ;; corner position == center
    X0 = ref_pos[*, ref_lind[0]]
      ;;; build the true reference grid by fitting, relative to X0
    tmp = ref_pos - X0 # replicate(1, np_ref)
    tmp_i = round(tmp/gridstep)  ;;; index
    ang = atan(ref_lpos[1, 1]-ref_lpos[1, 0], ref_lpos[0, 1]-ref_lpos[0, 0])
    if (ref_lpos[0, 1] - ref_lpos[0, 0]) lt (-gridstep) then ang += !pi
                                ;print, gridstep, ang
    par = mpfit('red_grid_func', [gridstep, ang], $
                functargs = {X:tmp_i, Y:tmp}, /quiet)
    ref_gridstep = par[0]
    print, inam, ' : regular grid width ', strtrim(par[0], 2), 'px, tilt ', strtrim(par[1]*!radeg, 2), 'deg'
    ;; check that there are no double index positions
    iix = replicate(1b, np_ref)
    for i=0, np_ref-2 do for j=i+1, np_ref-1 do $
       if total(tmp_i[*, i] eq tmp_i[*, j]) eq 2 then begin
      iix[i] = 0
      ref_pos[*, j] = (ref_pos[*, i]+ref_pos[*, j])/2
    endif
    iix = where(iix EQ 1)
    ref_pos = ref_pos[*, iix]
    ;; ref_lind is wrong if doublets were deleted! We need correct values for flip detection
    for i=0, 2 do ref_lind[i] = $
      where((ref_pos[0, *] eq ref_lpos[0, i]) and (ref_pos[1, *] eq ref_lpos[1, i]))
    reg_i = tmp_i[*, iix]
    np_ref = (size(reg_i, /dim))[1]
    reg_pos = x0#replicate(1, np_ref)-red_grid_func(par, x = reg_i, y = 0)
    
    printf, unit, string(this_ref_state.fpi_state, X0, par, $
                         form="(a-12,f7.2,'  ',f7.2,'   ',f7.3,f+7.3)")
    
    ;; compute warp matrices, forward (kxy) and backward(lxy)
    map = transpose([ref_pos, reg_pos])
    red_clean_warp, map, kx, ky, warp_dim
    red_clean_warp, map[*, [2, 3, 0, 1]], lx, ly, warp_dim
    
    red_append, alignments, { state:this_ref_state, map_x: kx, map_y: ky, $
                              revmap_x: lx, revmap_y: ly }
    
    ;; Now loop through dependent files and align to grid
    
    for i_dep = 0, n_idx-1 do begin
      state = states[idx[i_dep]]
      img = red_readdata(state.filename, /silent)
      is_pd = strmatch(state.camera, '*-D')
      ;; mask out bad corners also for Crisp PD image
      if state.camera eq 'Crisp-D' then $
         img *= (shift(dist(ref_siz), ref_siz/2) le 0.48*ref_siz[0])
      mask = img ge max(img)/(is_pd ? 8. : 10.)
      img_m = red_separate_mask(mask)
      img_np = max(img_m)
      img_area = lonarr(img_np)
      img_pos = fltarr(2, img_np)
      for i = 0, img_np-1 do begin
        ix = where(img_m eq (i+1))
        ix = [[ix mod ref_siz[0]], [ix / ref_siz[0]]]
        lu = (min(ix, dim = 1, max = ro)-3) > [0, 0]
        ro = (ro+3) < (ref_siz-1)
        bx = img[lu[0]:ro[0], lu[1]:ro[1]]
        bx_m = img_m[lu[0]:ro[0], lu[1]:ro[1]]
        if is_pd eq 1 then begin
          img_pos[*, i] = lu + red_centroid(bx_m)
          img_area[i] = total(bx_m GT 0)
        endif else begin
          img_pos[*, i] = lu + red_centroid(bx)
          img_area[i] = total(bx GT max(bx)/5.)
        endelse
      endfor
      if keyword_set(manual) then begin
        print, 'Mark the same three pinholes as for the reference'
        img_lind = red_markpinholes(img, img_pos)
        img_lpos = img_pos[*, img_lind]
        if red_lcheck(img_lpos, img_lind, gridstep, ratio=l_shape) ne 1 then stop
      endif else begin
        ix = reverse(sort(img_area))
        for i=0, 3 do for j=i+1, 4 do for k=j+1, 5 do begin
          img_lind = ix[[i, j, k]]
          img_lpos = img_pos[*, img_lind] 
          if red_lcheck(img_lpos, img_lind, gridstep, ratio=l_shape) eq 1 then goto, found_l2
        endfor                  ; i,j
        print, inam, ' : Unable to locate the reference pinholes. Try marking them manually'
        repeat begin
          img_lind = red_markpinholes(img, img_pos, title='Dependent Pinholes')
          img_lpos = img_pos[*, img_lind]
          if red_lcheck(img_lpos, img_lind, gridstep, ratio=l_shape) eq 1 then goto, found_l2
          print, inam, ' : This is not the 3 bigger pinholes!'
        endrep until 0
      endelse

found_l2:
      print, inam, ' : Found L for camera ', state.camera, ': ', gridstep
      XX0 = img_lpos[*, 0]

      ;; find the matching pinholes
      tmp = img_pos-xx0#replicate(1, img_np)
      tmp_i = round(tmp/gridstep)
      ;; this assumes large L side is horizontal.  Check/improve
      hflip = reg_i[0, ref_lind[1]] / tmp_i[0, img_lind[1]]
      vflip = reg_i[1, ref_lind[2]] / tmp_i[1, img_lind[2]]
      n_pairs = 0
      undefine, nmap
      for ph = 0, img_np-1 do begin
        ix = where( (reg_i[0, *] eq hflip*tmp_i[0, ph]) and $
                    reg_i[1, *] eq vflip*tmp_i[1, ph], nf)
        if nf eq 0 then continue
              ;;; have a match, add the pair to map
        red_append, nmap, [tmp[*, ph]+xx0, reg_pos[*, ix]] 
        n_pairs++
      endfor
      nmap = transpose(reform(nmap, 4, n_pairs))
      red_clean_warp, nmap, kx, ky, warp_dim
      red_clean_warp, nmap[*, [2, 3, 0, 1]], lx, ly, warp_dim
      
      red_append, alignments, { state:state, map_x: kx, map_y: ky, $
                                revmap_x: lx, revmap_y: ly }
    endfor    ;; loop dependent files
    if keyword_set(manual) then wdelete
  endfor      ;; loop reference files

  free_lun, unit

  save, file = output_file, alignments

end
