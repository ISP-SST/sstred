; docformat = 'rst'

;+
; Calculate how images need to be clipped for a rough alignment.
;
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :author:
; 
; 
; 
; :returns:
; 
; 
; :Params:
; 
;   
; 
; :Keywords:
;    
;    thres : in, optional, type=float, default=0.1
; 
;       Threshold for identifying a strong enough pinhole.
;    
;    extraclip : 
;    
;    
;    maxshift : 
;    
; 
; 
; :history:
;
;   2013-07-24 : Use red_show rather than show.
; 
;   2013-08-30 : MGL. Added support for logging. Let the subprogram
;                find out its own name. Get camtags from self.pinh_dir
;                rather than self.data_dir, in case (the first) data
;                directory does not have data from all cameras.
;
;   2013-09-18 : MGL. Make the transmitted camera the reference and
;                remove the refroot keyword.
;
;   2013-09-25 : MGL. Display all three properly clipped and flipped
;                pinhole images in a single window to make checking of
;                alignment and coarse alignment easier to check.
;                Bugfixed the extraclip feature. Read and display also
;                properly clipped gaintables in ordet to make it
;                easier to see if extra clipping is needed. 
;
;   2013-09-25 : MGL. Now checks for existence of gaintables before
;                reading. 
;
;   2016-08-23 : THI. Rename camtag to detector and channel to camera,
;                so the names match those of the corresponding SolarNet
;                keywords.
;
;-
PRO red::getalignclips, thres = thres, extraclip = extraclip, $
                        maxshift = maxshift

  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)                

  ;; Logging
  help, /obj, self, output = selfinfo 
  red_writelog, selfinfo = selfinfo

  if(n_elements(thres) eq 0) then tr = 0.1
  if(n_elements(maxshift) eq 0) THEN maxshift = 35
   
  ;; Search summed pinh images and camtag
  if(n_elements(extraclip) eq 0) then extraclip = [0L, 0L, 0L, 0L]
  if(n_elements(extraclip) eq 1) then extraclip = replicate(extraclip, 4)
  if(n_elements(extraclip) eq 2) then extraclip = [replicate(extraclip[0],2),replicate(extraclip[1],2)]

  self -> getdetectors
  camw = self.camwbtag
  camt = self.camttag
  camr = self.camrtag
  cams = [camw, camt, camr]
  labs = ['WB', 'NBT', 'NBR']

  fw = file_search(self.out_dir+'/pinh/' + camw +'.*.pinh', count = cw)
  ft = file_search(self.out_dir+'/pinh/' + camt +'.*.pinh', count = ct)
  fr = file_search(self.out_dir+'/pinh/' + camr +'.*.pinh', count = cr)

  ;; Get image states
  wstat = red_getstates_pinh(fw)
  tstat = red_getstates_pinh(ft, lam = lams)
  rstat = red_getstates_pinh(fr)

  ;; Select state to align
  allowed = [-1]
  for ii = 0L, n_elements(ft) -1 do BEGIN
     pr = where(rstat eq tstat[ii], count0)
     pw = where(wstat eq tstat[ii], count1)
     If(~(count0 gt 0) OR ~(count1 gt 0)) then continue
     print, red_stri(ii) +' '+tstat[ii]
     allowed = [temporary(allowed), ii]
  endfor
  if n_elements(allowed) gt 1 then allowed = allowed[1:*]

  toread = 0
  read, toread, prompt = inam+' : choose state to align: '

  pos = where(allowed eq toread, count)
  if count eq 0 then begin
     print, inam + ' : Error -> incorrect state number: ',toread
     return
  endif

  print, inam+' : selected state '+tstat[toread]
  pstate = tstat[toread]
  lam = lams[toread]
  pref = (strsplit(pstate, '.',/extract))[0]
  
  ;; Prepare for putting everything in arrays
  Ncams = 3
  ;; ref = cam_t
  icamref = 1
  ;; slaves = camw and camr
  icamnonrefs = [0, 2]
  
  ;; Read ref image
  ref = f0(ft[toread])
  dim = size(ref,/dim)

  pics = fltarr(dim[0], dim[1], Ncams)
  pics[*,*,icamref] = ref

  ;; Read slave images
  pos0 = where(wstat eq tstat[toread])
  pos1 = where(rstat eq tstat[toread])

  pics[*,*,icamnonrefs[0]]= f0(fw[pos0])
  pics[*,*,icamnonrefs[1]]= f0(fr[pos1])


  print, inam+' : images to be calibrated:'
  print, ' -> '+fw[pos0]
  print, ' -> '+ft[toread] + ' (reference)'
  print, ' -> '+fr[pos1]

  ;; Define the search space, rotation, shifts in x and y.
  rots = [0, 2, 5, 7]
  xshifts = [-1, 0, 1]
  yshifts = [-1, 0, 1]

  ;; Search space dimensions
  Nsearchr = n_elements(rots)   
  Nsearchx = n_elements(xshifts)
  Nsearchy = n_elements(yshifts)

  ;; Define arrays to store stuff
  np = 2

;  sh = intarr(2, Nsearchr)
  cor = fltarr(Nsearchr, Nsearchx, Nsearchy)      ; Correlation coefficients
  ssh = intarr(2, Ncams)
  i_rot = intarr(Ncams)
  i_xshift = intarr(Ncams)
  i_yshift = intarr(Ncams)

  ;; Search orientation!
;  dim -= 1

  ;; Find pinhole grid for reference image
  print, inam + ' : red_findpinholegrid ... ', format='(A,$)'
  red_findpinholegrid, pics[*,*,icamref], simx_orig, simy_orig, thres = thres
  print, 'done'

  gridspacing = median([deriv(simx_orig),deriv(simy_orig)])

  ;; Remove rows and columns of pinholes that are close enough to the
  ;; FOV borders in the reference channel, that they could be outside
  ;; the FOV in some other channel.
  border = gridspacing/2
  simx_orig = simx_orig(where(simx_orig gt border and simx_orig lt dim[0]-border))
  simy_orig = simy_orig(where(simy_orig gt border and simy_orig lt dim[1]-border))
  simx = simx_orig
  simy = simy_orig
  
  Nsimx = (size(simx_orig, /dim))[0] 
  Nsimy = (size(simy_orig, /dim))[0]  

  ;; Make array with peak intensities for all cameras
  peaks = fltarr(Nsimx, Nsimy, Ncams)
  for icam = 0, Ncams-1 do begin

     if icam eq icamref then begin
        simx = simx_orig
        simy = simy_orig
     endif else begin
        print, inam + ' : red_findpinholegrid ... ', format='(A,$)'
        red_findpinholegrid, pics[*,*,icam], simx, simy, thres=thres
        print, 'done'
        simx = simx(where(simx gt border and simx lt dim[0]-border))
        simy = simy(where(simy gt border and simy lt dim[1]-border))
       ;; May have to reduce the grid size
        Nsimx = Nsimx < (size(simx, /dim))[0] 
        Nsimy = Nsimy < (size(simy, /dim))[0] 
     endelse

     for isim = 0, Nsimx-1 do begin
        for jsim = 0, Nsimy-1 do begin
           peaks[isim, jsim, icam] = max(pics[simx[isim]-border:simx[isim]+border $
                                              , simy[jsim]-border:simy[jsim]+border $
                                              , icam])
        endfor                  ; jsim
     endfor                     ; isim
  endfor                        ; icam

  ;; Possibly reduce peaks array size
  peaks = peaks[0:Nsimx-1, 0:Nsimy-1, *]
  
  ;; Match reference peaks with other cameras
  refpeaks = peaks[Nsearchx/2:Nsimx-Nsearchx/2-1, Nsearchy/2:Nsimy-Nsearchy/2-1, icamref]
  for im = 0, Ncams-2 do begin
     
     icam = icamnonrefs[im]

     ;; Find matching orientation and shift.
     for irot = 0, Nsearchr-1 do begin        ; Loop over orientations

        rotpeaks = rotate(peaks[*, *, icam], rots[irot])
      
        for ix = 0, Nsearchx-1 do begin     ; Loop over x shifts
           for iy = 0, Nsearchy-1 do begin ; Loop over y shifts

              shiftpeaks = rotpeaks[xshifts[ix]+Nsearchx/2:xshifts[ix]+Nsimx-Nsearchx/2-1 $
                                    , yshifts[iy]+Nsearchy/2:yshifts[iy]+Nsimy-Nsearchy/2-1]
              cor[irot, ix, iy] = correlate(refpeaks, shiftpeaks)
;              
           endfor               ; iy
        endfor                  ; ix

     endfor                     ; irot

     ;; Find best correlation
     cm = max(cor, maxloc)
     maxindx = array_indices(cor,maxloc)

     ;; Find rot and shifts that correspond to best correlation
     i_rot[icam] = rots[maxindx[0]]
     i_xshift[icam] = xshifts[maxindx[1]]
     i_yshift[icam] = yshifts[maxindx[2]]

     ;; Use the found orientation and find shifts
     p1 = rotate(reform(pics[*, *, icam]), i_rot[icam])

     mxsh = maxshift+max(abs([i_xshift[icam], i_yshift[icam]]))*gridspacing
     print, 'Max shift allowed: ', mxsh
     print, maxindx

     ssh[*, icam] = red_shc(pics[*,*,icamref], p1, RANGE = mxsh)

     print, inam+' '+cams[icam]+' : orientation ', strtrim(i_rot[icam], 2), $
            ' -> shift: x,y=', strtrim(ssh[0, icam], 2), ', ', strtrim(ssh[1, icam], 2)

  endfor                        ; im (icam)

  ;; Clip images
  ssh = -ssh
  ssh[0, *] -= min(ssh[0, *])
  ssh[1, *] -= min(ssh[1, *])
  sx = dim[0] - max(ssh[0, *])
  sy = dim[1] - max(ssh[1, *])
  sx = 2 * (sx / 2)
  sy = 2 * (sy / 2)
  cl = rebin(ssh, 4, Ncams, /sam)
  cl[1, *]+= sx-1
  cl[3, *]+= sy-1

  for i = 0, Ncams-1 do begin
     if i_rot[i] eq 2 or i_rot[i] eq 5 then begin
        cl[0:1, i] = dim[0]-cl[0:1, i] 
     endif else begin
        cl[0:1, i] += 1
     endelse
     if i_rot[i] eq 2 or i_rot[i] eq 7 then begin
        cl[2:3, i] = dim[1]-cl[2:3, i] 
     endif else begin
        cl[2:3, i] += 1
     endelse
  endfor

  ;; Extra border clip?
  if(total(extraclip) gt 0) then begin
     print, inam + ' : extra border clip in X -> '+red_stri(extraclip[0])+' pixels'
     print, inam + ' : extra border clip in Y -> '+red_stri(extraclip[1])+' pixels'

     for icam = 0, Ncams-1 do begin
        if(cl[0,icam] lt cl[1,icam]) then begin
           cl[0,icam] += extraclip[0]
           cl[1,icam] -= extraclip[1]
        endif else begin
           cl[0,icam] -= extraclip[0]
           cl[1,icam] += extraclip[1]
        endelse
        if(cl[2,icam] lt cl[3,icam]) then begin
           cl[2,icam] += extraclip[2]
           cl[3,icam] -= extraclip[3]
        endif else begin
           cl[2,icam] -= extraclip[2]
           cl[3,icam] += extraclip[3]
        endelse
     endfor

     sx = abs(cl[1,0] - cl[0,0]) + 1L
     sy = abs(cl[3,0] - cl[2,0]) + 1L

  endif

  ;; Align clips
  acl = 'ALIGN_CLIP='+strjoin(strtrim(cl, 2), ',')
  file_mkdir, self.out_dir +'/calib'
  openw, lun, self.out_dir+'/calib/align_clips.'+pref+'.txt', /get_lun
  for icam = 0L, Ncams-1 do begin
     printf, lun, acl[icam]
     print, '  -> Align CLIP: '+acl[icam]
  endfor
  free_lun, lun

  refrot = 0
  save, file = self.out_dir + '/calib/align_clips.'+pref+'.sav', acl, cl, refrot, sx, sy, ssh

  print, 'Please check that the displayed images are aligned and oriented properly.'

  szx_disp = dim[0]/2
  szy_disp = dim[1]/2
  disp_spacing = 5

  title = '        '
  for icam = 0, Ncams-1 do title = title+labs[icam]+' : '+cams[icam]+'        '
  
  window, xs = szx_disp*Ncams+(Ncams+1)*disp_spacing, ys = szy_disp+2*disp_spacing, title = title, 0
  erase, 255
  for icam = 0, Ncams-1 do begin
     dispim = red_clipim(pics[*,*,icam], cl[*,icam])
     dispim = congrid(dispim, szx_disp, szy_disp, cubic = -0.5)
     tvscl, dispim, icam*(szx_disp+disp_spacing)+disp_spacing, disp_spacing
  endfor

  ;; Display gain tables?
  window, xs = szx_disp*Ncams+(Ncams+1)*disp_spacing, ys = szy_disp+2*disp_spacing, title = title, 1
  erase, 255
  gains = fltarr(dim[0], dim[1], Ncams)

;  print, strjoin(strsplit(ft[toread], '/pinh/', /extr), '/gaintables/')
  gname = strjoin(strsplit(ft[toread], '\.pinh', /extr,/preserve,/rege), '.gain')
  gname = strjoin(strsplit(gname, '/pinh/', /extr,/preserve,/rege), '/gaintables/')
  if file_test(gname) then begin
     gains[*,*,icamref]  = f0(gname)
  endif else begin
     print, inam+' : Could not find gaintable '+gname
     print, inam+' : Need to run a->makegains'
  endelse
  gname = strjoin(strsplit(fr[pos1], '\.pinh', /extr,/preserve,/rege), '.gain')
  gname = strjoin(strsplit(gname, '/pinh/', /extr,/preserve,/rege), '/gaintables/')
  if file_test(gname) then begin
     gains[*,*,icamnonrefs[1]] = f0(gname)
  endif else begin
     print, inam+' : Could not find gaintable '+gname
     print, inam+' : Need to run a->makegains'
  endelse
  ;; The WB gain file name has no tuning info
  gname = strjoin(strsplit(fw[pos0], '\.pinh', /extr,/preserve,/rege), '.gain')
  gname = strjoin(strsplit(gname, '/pinh/', /extr,/preserve,/rege), '/gaintables/')
  gname = strsplit(gname, '.', count = nn, /extr)
  gname = strjoin([gname[0:nn-4], gname[nn-1]], '.')
  if file_test(gname) then begin
     gains[*,*,icamnonrefs[0]] = f0(gname)
  endif else begin
     print, inam+' : Could not find gaintable '+gname
     print, inam+' : Need to run a->makegains'
  endelse

  for icam = 0, Ncams-1 do begin
     dispim = red_clipim(gains[*,*,icam], cl[*,icam])
     dispim = congrid(dispim, szx_disp, szy_disp, cubic = -0.5)
     tvscl, dispim, icam*(szx_disp+disp_spacing)+disp_spacing, disp_spacing
  endfor

  ;; Print out some info about blocked rows/columns.
  gg = red_clipim(gains[*,*,icamref], cl[*,icamref]) ne 0 
  blockedrows = where(total(gg,1) lt .5, Nblock)
  if Nblock ne 0 then begin
     print, 'The following rows are significantly blocked: '
     print, blockedrows
  endif
  blockedcolumns = where(total(gg,2) lt .5, Nblock)
  if Nblock ne 0 then begin
     print, 'The following columns are significantly blocked: '
     print, blockedcolumns
  endif

  return
end
