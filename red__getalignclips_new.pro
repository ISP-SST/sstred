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
;    thres : 
;    
;    
;    extraclip : 
;    
;    
;    maxshift : 
;    
;    
;    pref : in, optional, type=string
;     
;      Indicate the prefilter you want to calculate the clips for,
;      Default is to do it for all prefilters there is data for.
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
;   2014-01-15 : MGL. Make it work with the file name conventions used
;                when pinholes are only summed for a single state per
;                camera and prefilter. New keyword: pref.
; 
;   2014-01-16 : MGL. Use gain file saved in the sumpinh step. Make a
;                prefilter loop. If more than one pinhole image per
;                prefilter, select the one that would be brightest.
;
;   2014-01-22 : MGL. Replace red_getstates with red_extractstates.
;                Implemented automatic selection of a bright tuning
;                when there are more than one. Fixed display of gain
;                tables as saved by red::sumpinh. 
;
;   2014-04-07 : THI. Use red_strreplace. Correct path for gains.
;
;-
PRO red::getalignclips_new, thres = thres $
                            , pref = pref $
                            , extraclip = extraclip $
                            , maxshift = maxshift

  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])
  
  ;; Logging
  help, /obj, self, output = selfinfo 
  red_writelog, selfinfo = selfinfo

  if(n_elements(thres) eq 0) then tr = 0.05 
  if(n_elements(maxshift) eq 0) THEN maxshift = 35
   
  ;; Prepare for putting the different cameras in arrays
  Ncams = 3
  ;; ref = camt
  icamref = 1
  ;; slaves = camw and camr
  icamnonrefs = [0, 2]

  ;; Search summed pinh images and camtag
  if(n_elements(extraclip) eq 0) then extraclip = [0L, 0L, 0L, 0L]
  if(n_elements(extraclip) eq 1) then extraclip = replicate(extraclip, 4)
  if(n_elements(extraclip) eq 2) then extraclip = [replicate(extraclip[0],2),replicate(extraclip[1],2)]

  self -> getcamtags, dir = self.pinh_dir
  camw = self.camwbtag
  camt = self.camttag
  camr = self.camrtag
  cams = [camw, camt, camr]
  labs = ['WB', 'NBT', 'NBR']

  ;; Selected prefilter or all prefilters?
  fw = file_search(self.out_dir+'/pinh/' + camw +'.*.pinh', count = cw)
  if cw eq 0 then begin
     print, inam, ' : ERROR : No wideband pinholes found in ', self.out_dir+'/pinh/'
     retall
  endif
  prefilters = strarr(cw)
  for ii = 0, cw-1 do begin
     ;; Remove directory
     prefilters[ii] = (strsplit(fw[ii], '/', /extract, count=nn))[nn-1]
     ;; Get the prefilter
     prefilters[ii] = (strsplit(prefilters[ii], '.', /extract))[1]
  endfor
  prefilters = prefilters[uniq(prefilters, sort(prefilters))]
  if n_elements(pref) ne 0 then begin
     indx = where(prefilters eq pref)
     if max(indx) eq -1 then begin
        print, inam+' : WARNING : Keyword pref does not match any pinhole file names: ', pref
        return
     endif
  endif
  Npref = n_elements(prefilters)

  for ipref = 0, Npref-1 do begin

     wfiles = file_search(self.out_dir+'/pinh/' + camw +'.'+prefilters[ipref]+'*.pinh', count = cw)
     tfiles = file_search(self.out_dir+'/pinh/' + camt +'.'+prefilters[ipref]+'*.pinh', count = ct)
     rfiles = file_search(self.out_dir+'/pinh/' + camr +'.'+prefilters[ipref]+'*.pinh', count = cr)
  
     if (cw ne ct) or (cw ne cr) then begin
        
        print, inam+' : Mismatch in available states for the different cameras.'
        print, '  Number of states for cameras WB, NBT, NBR:', cw, ct, cr
        retall

     endif

     if cw ne 1 then begin

        ;; More than a single pinhole image per camera, select the one
        ;; expected to be the brightest based on its (tuning) state.

        red_extractstates, wfiles, wav = wav, dwav = dwav, /basename
        tmp = max(abs(dwav-double(wav)), ml)
              
        wfiles = wfiles[ml]
        tfiles = tfiles[ml]
        rfiles = rfiles[ml]

     endif

     pnames = [wfiles, tfiles, rfiles]

     print, inam+' : images to be calibrated:'
;
     for icam = 0, Ncams-1 do begin
        ostring = ' -> '+pnames[icam]
        if icam eq icamref then ostring += ' (reference)'
        print, ostring
     endfor

     ref = f0(pnames[icamref])
     dim = size(ref,/dim)
     pics = fltarr(dim[0], dim[1], Ncams)
     pics[*,*,icamref] = ref

     ;; Read slave images
     pics[*,*,icamnonrefs[0]]= f0(pnames[icamnonrefs[0]])
     pics[*,*,icamnonrefs[1]]= f0(pnames[icamnonrefs[1]])

     ;; Define the search space for rotation and shifts in x and y.
     rots = [0, 2, 5, 7]
     xshifts = [-1, 0, 1]
     yshifts = [-1, 0, 1]

     ;; Search space dimensions
     Nsearchr = n_elements(rots)   
     Nsearchx = n_elements(xshifts)
     Nsearchy = n_elements(yshifts)

     ;; Define arrays to store stuff
     cor = fltarr(Nsearchr, Nsearchx, Nsearchy) ; Correlation coefficients
     ssh = intarr(2, Ncams)
     i_rot = intarr(Ncams)
     i_xshift = intarr(Ncams)
     i_yshift = intarr(Ncams)

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
           endfor               ; jsim
        endfor                  ; isim

     endfor                     ; icam

     ;; Possibly reduce peaks array size
     peaks = peaks[0:Nsimx-1, 0:Nsimy-1, *]
     
     ;; Match reference peaks with other cameras
     refpeaks = peaks[Nsearchx/2:Nsimx-Nsearchx/2-1, Nsearchy/2:Nsimy-Nsearchy/2-1, icamref]
     for im = 0, Ncams-2 do begin
        
        icam = icamnonrefs[im]

        ;; Find matching orientation and shift.
        for irot = 0, Nsearchr-1 do begin ; Loop over orientations

           rotpeaks = rotate(peaks[*, *, icam], rots[irot])
           
           for ix = 0, Nsearchx-1 do begin    ; Loop over x shifts
              for iy = 0, Nsearchy-1 do begin ; Loop over y shifts

                 shiftpeaks = rotpeaks[xshifts[ix]+Nsearchx/2:xshifts[ix]+Nsimx-Nsearchx/2-1 $
                                       , yshifts[iy]+Nsearchy/2:yshifts[iy]+Nsimy-Nsearchy/2-1]
                 cor[irot, ix, iy] = correlate(refpeaks, shiftpeaks)
              
              endfor            ; iy
           endfor               ; ix
           
        endfor                  ; irot

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

     endfor                     ; im (icam)

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
     openw, lun, self.out_dir+'/calib/align_clips.'+prefilters[ipref]+'.txt', /get_lun
     for icam = 0L, Ncams-1 do begin
        printf, lun, acl[icam]
        print, '  -> Align CLIP: '+acl[icam]
     endfor
     free_lun, lun

     refrot = 0
     save, file = self.out_dir + '/calib/align_clips.'+prefilters[ipref]+'.sav', acl, cl, refrot, sx, sy, ssh

     print, 'Please check that the displayed images are aligned and oriented properly.'

     szx_disp = dim[0]/2
     szy_disp = dim[1]/2
     disp_spacing = 5

     title = '        '+prefilters[ipref]+':      '
     for icam = 0, Ncams-1 do title = title+labs[icam]+' = '+cams[icam]+'        '
     
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
     for icam = 0, Ncams-1 do begin

        gname = red_strreplace(pnames[icam], '.pinh', '.gain')
        gname = red_strreplace(gname, '/pinh/', '/gaintables/')
        if cams[icam] EQ self.camwbtag then begin
            ;; The WB gain file name has no tuning info
            gname = strsplit(gname, '.', count = nn, /extr)
            gname = strjoin([gname[0:nn-4], gname[nn-1]], '.')
        endif
        if file_test(gname) then begin
           gains[*,*,icam]  = f0(gname)
        endif else begin
           print, inam+' : Could not find gaintable '+gname
        endelse

        dispim = red_clipim(gains[*,*,icam], cl[*,icam])
        dispim = congrid(dispim, szx_disp, szy_disp, cubic = -0.5)
        tvscl, dispim, icam*(szx_disp+disp_spacing)+disp_spacing, disp_spacing

     endfor

     ;; Print out some info about blocked rows/columns.
     ;; One could add intelligence here, suggesting a new setting for extraclip.      <-------------
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

  endfor                        ; ipref
  
  
end
