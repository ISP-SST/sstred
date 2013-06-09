;+
;
;
; :history:
;
;   2013-05-?? : Keywords dx and dy added by MGL, use these if there
;                are large misalignments, like a whole row or column
;                of pinholes. Set dx to the shifts needed to
;                approximately align the NB cameras in X to the WB
;                (ref) after applying the correct flipping. Similar
;                for dy.
;
;-
PRO red::getalignclips, refrot = refrot, thres = thres, extraclip = extraclip, $
                        maxshift = maxshift, $
                        dx = dx, dy = dy  

  if n_elements(dx) eq 0 then dx = [0, 0] else dx = round(dx)
  if n_elements(dy) eq 0 then dy = [0, 0] else dy = round(dy)

  if(n_elements(thres) eq 0) then tr = 0.1 ; Not used?? MGL
  if(n_elements(maxshift) eq 0) THEN maxshift = 35
   
  ;; Procedure pinhcalib based on Pit Sutterlin's setup_ph.pro
   
  ;; Seach summed pinh images and camtag
   
  inam = 'red::getalignclips : '
  if(n_elements(extraclip) eq 0) then extraclip = [0L, 0L, 0L, 0L]
  if(n_elements(extraclip) eq 1) then extraclip = replicate(extraclip, 4)
  if(n_elements(extraclip) eq 2) then extraclip = [replicate(extraclip[0],2),replicate(extraclip[1],2)]

  ;; if(self.dopinh) then begin
  ;;    spawn, 'find ' +self.pinh_dir+'/'+self.camt+'/|grep cam',f
  ;;    camt = red_camtag(f[0])
  ;;    ;
  ;;    spawn, 'find ' +self.pinh_dir+'/'+self.camwb+'/|grep cam',f
  ;;    camw = red_camtag(f[0])
  ;;    ;
  ;;    spawn, 'find ' +self.pinh_dir+'/'+self.camr+'/|grep cam',f
  ;;    camr = red_camtag(f[0])    
  ;; endif else begin
  ;;    print, inam+' ERROR, undefined pinh_dir'
  ;;    return
  ;; endelse

  self -> getcamtags, dir = self.data_dir
  camt = self.camttag
  camr = self.camrtag
  camw = self.camwbtag
                                ;
  ft = file_search(self.out_dir+'/pinh/' + camt +'.*.pinh', count = ct)
  fr = file_search(self.out_dir+'/pinh/' + camr +'.*.pinh', count = cr)
  fw = file_search(self.out_dir+'/pinh/' + camw +'.*.pinh', count = cw)
                                ;
                                ; Get image states
  tstat = red_getstates_pinh(ft, lam = lams)
  rstat = red_getstates_pinh(fr)
  wstat = red_getstates_pinh(fw)
                                ;
                                ; Select state to align
                                ;
  allowed = [-1]
  for ii = 0L, n_elements(ft) -1 do BEGIN
     pr = where(rstat eq tstat[ii], count0)
     pw = where(wstat eq tstat[ii], count1)
     If(~(count0 gt 0) OR ~(count1 gt 0)) then continue
     print, red_stri(ii) +' '+tstat[ii]
     allowed = [temporary(allowed), ii]
  endfor
  if n_elements(allowed) gt 1 then allowed = allowed[1:*]
                                ;
  toread = 0
  read, toread, prompt = inam+'choose state to align: '
                                ;
  pos = where(allowed eq toread, count)
  if count eq 0 then begin
     print, inam + 'Error -> incorrect state number: ',toread
     return
  endif
                                ;
  print, inam+'selected state '+tstat[toread]
  pstate = tstat[toread]

  pref = (strsplit(pstate, '.',/extract))[0]
  
                                ;
                                ; load states 
                                ; ref = cam_t
                                ; slaves = camr and camw
  pos2 = where(wstat eq tstat[toread])
  ref = f0(fw[pos2])
  refs = ref
  dim = size(ref,/dim)
                                ;
  pics = fltarr(dim[0], dim[1], 2)
  pics[*,*,0]= f0(ft[toread])
  pos1 = where(rstat eq tstat[toread])
  pics[*,*,1]= f0(fr[pos1])
                                ;
  print, inam+'images to be calibrated:'
  print, ' -> '+ft[toread]
  print, ' -> '+fr[pos1]
  print, ' -> '+fw[pos2]
  tfil = ft[toread]
  rfil = fr[pos1]
  wfil = fw[pos2]
  lam = lams[toread]
                                ;
  rots = [0, 2, 5, 7]
                                ;
                                ; Rotate ref?
                                ;
  if(keyword_set(refrot)) then begin
     pp = where(rots eq refrot, count)
     if count eq 0 then begin
        print, inam+'invalid supplied refrot'
        print, inam+'valid orientations:'
        print,' -> 0 - same,   2 - flip X+Y,   5 - flip X,   7 - flip Y'
        print, inam+'ignoring refrot keyword!'
        wait, 2
     endif else begin
        print, inam+'rotating reference to position -> '+red_stri(refrot)
        ref = rotate(temporary(ref), refrot)
     endelse
  endif else refrot = 0
                                ;
                                ; define arrays to store stuff
                                ;
  np = 2
  sr = dim
  nr = n_elements(rots)
  sh = intarr(2, nr)
  cor = fltarr(nr)
  ssh = intarr(2, np+1)
  i_rot = intarr(np+1)
  i_rot[0] = refrot
                                ;
                                ; Search orientation!
                                ;
  dim-= 1
                                ;
  ;; Find pinhole grid for reference image
  print, inam + 'red_findpinholegrid ... ', format='(A,$)'
  red_findpinholegrid, ref, simx_orig, simy_orig
  print, 'done'

  ;; Remove rows and columns of pinholes that are close enough to the
  ;; FOV borders in the reference channel, that they could be outside
  ;; the FOV in some other channel.
  simx_orig = simx_orig(where(simx_orig gt maxshift and simx_orig lt (size(ref, /dim))[0]-maxshift))
  simy_orig = simy_orig(where(simy_orig gt maxshift and simy_orig lt (size(ref, /dim))[1]-maxshift))
  
  Ngridx = (size(simx_orig, /dim))[0]
  Ngridy = (size(simy_orig, /dim))[0]

  ;; Make array with ref pinhole peak intensities.
  refpeaks = fltarr(Ngridx, Ngridy)
  simx = simx_orig
  simy = simy_orig
  for igrid = 0, Ngridx-1 do begin
     for jgrid = 0, Ngridy-1 do begin
        thispeak = ref[simx[igrid]-maxshift:simx[igrid]+maxshift, simy[jgrid]-maxshift:simy[jgrid]+maxshift]
        refpeaks[igrid, jgrid] = max(thispeak)
     endfor
  endfor

  FOR im = 0, np-1 DO BEGIN

     ;; Make array with pinhole peak intensities
     peaks = fltarr(Ngridx, Ngridy)

     ;; This pic, optionally shifted.  MGL May 2013
     thispic = shift(pics[*, *, im], dx[im], dy[im])

     ;; Find matching orientation.
     ;; We could implement ABS DIFF on these matrices, this would
     ;; remove the need for the "maxshift" parameter. Which would
     ;; probably be nice for the blue tilt filter, which shifts images
     ;; quite a bit for large angles.
     FOR n = 0, nr-1 DO BEGIN
        if rots[n] eq 0 or rots[n] eq 7 then begin 
           simx = simx_orig
        endif else begin
           simx = (size(ref, /dim))[0]-reverse(simx_orig)
        endelse
        if rots[n] eq 0 or rots[n] eq 5 then begin
           simy = simy_orig
        endif else begin
           simy = (size(ref, /dim))[0]-reverse(simy_orig)
        endelse

        for igrid = 0, Ngridx-1 do begin
           for jgrid = 0, Ngridy-1 do begin
              ;;thispeak = pics[simx[igrid]-maxshift:simx[igrid]+maxshift, simy[jgrid]-maxshift:simy[jgrid]+maxshift, im]
              thispeak = thispic[simx[igrid]-maxshift:simx[igrid]+maxshift, simy[jgrid]-maxshift:simy[jgrid]+maxshift] 
              peaks[igrid, jgrid] = max(thispeak)
           endfor
        endfor

                                ;  tvscl, rebin(refpeaks,10*Ngridx,10*Ngridy,/samp) , 0
                                ;  tvscl, rebin(peaks,10*Ngridx,10*Ngridy,/samp) , 2

        rotpeaks = rotate(peaks, rots[n])
        cor[n] = correlate(refpeaks, rotpeaks)

     ENDFOR
     cm = max(cor, w)
     i_rot[im+1] = rots[w]


     ;; Use the found orientation and find shifts
     p1 = rotate(reform(pics[*, *, im]), i_rot[im+1])
     ssh[*, im+1] = shc(ref, p1, RANGE = maxshift)

     ;; Compensate for optional shifts. MGL May 2013
     ssh[*, im+1] = ssh[*, im+1] + [dx[im], dy[im]]

     if(im eq 0) then cam=camt else cam=camr
     print, inam+cam+' orientation ', strtrim(i_rot[im+1], 2), $
            ' -> shift: x,y=', strtrim(ssh[0, im+1], 2), ', ', strtrim(ssh[1, im+1], 2)


     if 1 then begin            ; This is only for testing, can be removed!
        n = w
;        p1 = reform(pics[*, *, im])
        p1 = rotate(reform(pics[*, *, im]), rots[n])
        sh[*, n] = shc(ref, p1, RANGE = 50)
        r1 = ref[sh[0, n] > 0:dim[0]+(sh[0, n] < 0), $
                 sh[1, n] > 0:dim[1]+(sh[1, n] < 0)]
        p1 = p1[(-sh[0, n]) > 0:dim[0]-(sh[0, n] > 0), $
                (-sh[1, n]) > 0:dim[1]-(sh[1, n] > 0)]
                                ;  window, 0, xsize = (size(ref, /dim))[0], ysize = (size(ref, /dim))[1]
                                ;  tvscl, r1
                                ;  window, 1, xsize = (size(p1, /dim))[0], ysize = (size(p1, /dim))[1]
                                ; tvscl, p1  

     endif

  ENDFOR
                                ;
                                ; Clip images
                                ;
  ssh= -ssh
  ssh[0, *]-= min(ssh[0, *])
  ssh[1, *]-= min(ssh[1, *])
  sx = sr[0]-max(ssh[0, *])
  sy = sr[1]-max(ssh[1, *])
  sx = 2 * (sx / 2)
  sy = 2 * (sy / 2)
  cl = rebin(ssh, 4, np+1, /sam)
  cl[1, *]+= sx-1
  cl[3, *]+= sy-1
                                ;
  FOR i = 0, np DO BEGIN
     IF i_rot[i] EQ 2 OR i_rot[i] EQ 5 THEN cl[0:1, i] = sr[0]-cl[0:1, i] $
     ELSE cl[0:1, i] += 1
     IF i_rot[i] EQ 2 OR i_rot[i] EQ 7 THEN cl[2:3, i] = sr[1]-cl[2:3, i] $
     ELSE cl[2:3, i] += 1
  ENDFOR

  ;;
  ;; Extra border clip?
  ;;
  if(total(extraclip) gt 0) then begin
     print, inam + 'extra border clip in X -> '+red_stri(extraclip[0])+' pixels'
     print, inam + 'extra border clip in Y -> '+red_stri(extraclip[1])+' pixels'

     for ii = 0, 2 do begin
        if(cl[0,ii] lt cl[1,ii]) then begin
           cl[0,ii] += extraclip[0]
           cl[1,ii] -= extraclip[1]
        endif else begin
           cl[0,ii] -= extraclip[0]
           cl[1,ii] += extraclip[1]
        endelse
        if(cl[2,ii] lt cl[3,ii]) then begin
           cl[2,ii] += extraclip[2]
           cl[3,ii] -= extraclip[3]
        endif else begin
           cl[2,ii] -= extraclip[2]
           cl[3,ii] += extraclip[3]
        endelse
        
     endfor

     sx = abs(cl[1,0] - cl[0,0]) + 1L
     sy = abs(cl[3,0] - cl[2,0]) + 1L

     dum = red_clipim(pics[*,*,0], cl[*,1])
     show, bytscl(dum, 0, 20), /nosc


     
  endif

;tighttv, red_clipim(pics[*,*,0], cl[*,1]), 0
;tighttv, red_clipim(rotate(temporary(ref), refrot), cl[*,0]), 1
;print, 555
;stop
                                ;
                                ; Align clips
  acl = 'ALIGN_CLIP='+strjoin(strtrim(cl, 2), ',')
  file_mkdir, self.out_dir +'/calib'
  openw, lun, self.out_dir+'/calib/align_clips.'+pref+'.txt', /get_lun
  for ii = 0L, 2 do begin
     printf, lun, acl[ii]
     print, '  -> Align CLIP: '+acl[ii]
  endfor
  free_lun, lun
                                ;
  save, file = self.out_dir + '/calib/align_clips.'+pref+'.sav', acl, cl, refrot, sx, sy, ssh
                                ;
  return
end
