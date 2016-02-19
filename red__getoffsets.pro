; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :author:
; 
;   Jaime de la Cruz Rodriguez (Based on Pit Sutterlin's setup_ph.pro)
; 
; 
; 
; :returns:
; 
; 
; :Params:
; 
; 
; :Keywords:
; 
;    thres  : 
;   
;   
;   
;    state  : 
;   
;   
;   
;    pref : 
;   
;   
;   
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
; 
;   2013-09-09 : MGL. Added support for logging. Let the subprogram
;                find out its own name. Get camtags from self.pinh_dir
;                rather than self.data_dir, in case (the first) data
;                directory does not have data from all cameras.
;
;   2014-01-23 : MGL. Use red_extractstates instead of red_getstates
;                and local extraction of info from file names.
; 
;   2014-04-09 : TH. Use pinh_align directory.
;
;   2014-04-26 : MGL. Bugfix: lam has to be a string when written to
;                config file.
;
;   2016-02-17 : MGL. Include prefilter info when complaining about
;                not finding align clips file.
;
;-
pro red::getoffsets, thres = thres, state = state, pref=pref
  if(~keyword_set(thres)) then tr = 0.1 else tr = thres

  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])
  
  ;; Logging
  help, /obj, self, output = selfinfo 
  red_writelog, selfinfo = selfinfo

  ;; clips exist?
   
  if(~self.dopinh) then begin
     print, inam+' : ERROR, undefined pinh_dir'
     return
  endif
  if(n_elements(pref) eq 0) then pref = '*'

;  self -> getcamtags, dir = self.data_dir
  self -> getcamtags, dir = self.pinh_dir
  camt = self.camttag
  camr = self.camrtag
  camw = self.camwbtag

  tfiles = file_search(self.out_dir+'/pinh_align/' + camt + '.' + pref + '.*.pinh', count = ct)
  rfiles = file_search(self.out_dir+'/pinh_align/' + camr + '.' + pref + '.*.pinh', count = cr)
  wfiles = file_search(self.out_dir+'/pinh_align/' + camw + '.' + pref + '.*.pinh', count = cw)

  ;; Get image states
  red_extractstates, tfiles, /basename, fullstate = tstat, lambda = lams
  red_extractstates, rfiles, /basename, fullstate = rstat 
  red_extractstates, wfiles, /basename, fullstate = wstat 
;  tstat = red_getstates_pinh(tfiles, lam = lams)
;  rstat = red_getstates_pinh(rfiles)
;  wstat = red_getstates_pinh(wfiles)
   
  ;; Select state to align
   
  allowed = [-1]
  for ii = 0L, n_elements(tfiles) -1 do BEGIN
     if(keyword_set(state)) then begin
        if(tstat[ii] ne state) then begin
           print, inam + ' : Skipping state -> '+tstat[ii]
           continue
        endif
     endif
     pref = (strsplit(tstat[ii], '.',/extract))[0]
     
     IF(~file_test(self.out_dir+'/calib/align_clips.'+pref+'.sav')) THEN BEGIN
        print, inam + ' : ERROR -> align clips file not found for '+pref
        print, inam + ' :       -> you must execute red::getalignclips first!'
        continue
     ENDIF ELSE restore, self.out_dir+'/calib/align_clips.'+pref+'.sav'
     print, inam+' loading -> '+self.out_dir+'calib/align_clips.'+pref+'.sav'
     toread = ii
     pr = where(rstat eq tstat[ii], count0)
     pw = where(wstat eq tstat[ii], count1)
     If(~(count0 gt 0) OR ~(count1 gt 0)) then continue
     ;; print, red_stri(ii) +' '+tstat[ii]
     allowed = [temporary(allowed), ii]
     if n_elements(allowed) gt 1 then allowed = allowed[1:*]
      
     ;;toread = 0
     ;;read, toread, prompt = inam+' : choose state to align: '
      
     pos = where(allowed eq toread, count)
     if count eq 0 then begin
        print, inam + ' : Error -> incorrect state number: ',toread
        return
     endif
      
     ;;print, inam+' : selected state '+tstat[toread]
     pstate = tstat[toread]
                                
     ;; load states 
     ;; ref = cam_t
     ;; slaves = camr and camw
     pos2 = where(wstat eq tstat[toread])
     ref = f0(wfiles[pos2])
     dim = size(ref,/dim)
                                
     pics = fltarr(dim[0], dim[1], 2)
     pics[*,*,0]= f0(tfiles[toread])
     pos1 = where(rstat eq tstat[toread])
     pics[*,*,1]= f0(rfiles[pos1])
                                
     print, inam+' : images to be calibrated:'
     print, ' -> '+tfiles[toread]
     print, ' -> '+rfiles[pos1]
     print, ' -> '+wfiles[pos2]
     tfil = tfiles[toread]
     rfil = rfiles[pos1]
     wfil = wfiles[pos2]
     lam = lams[toread]
                                
     rots = [0, 2, 5, 7]
                                
     ;; define arrays to store stuff
      
     np = 2
     sr = dim
     nr = n_elements(rots)
     sh = intarr(2, nr)
     cor = fltarr(nr)
     i_rot = intarr(np+1)
     i_rot[0] = refrot
     
     ;; Find pinhole locations! 
      
     ref = red_clipim(temporary(ref), cl[*,0])
    
     ;; All pixels within each pinhole get the same unique number
      
     mask = red_separate_mask(ref gt tr * max(ref))
      
     ;;  Number of pinholes found
     nph = max(mask)
      
     ;; Compute PH positions (computing the CG)
     cc = fltarr(2, nph)
     FOR i = 0, nph-1 DO cc[*, i] = red_com(mask EQ i+1)
     cx = reform(cc[0, *])
     cy = reform(cc[1, *])
      
     ;; sort values - PHs should be aligned hor/vert, so this will
     ;;               give a clear step shape
     cx = cx[sort(cx)]
     cy = cy[sort(cy)]

     ;; Locate the steps and average the values of each step
     dcx = cx[1:*] - cx
     scx = [-1, where(dcx GT avg(dcx), nx), nph-1]
     simx = intarr(nx+1)
     FOR i=1, nx+1 DO simx[i-1] = round(mean(cx[scx[i-1]+1:scx[i]]))
     simx = simx[where((simx GT 32) AND (simx LT sx-32))]
      
     ;; Y axis
     dcy = cy[1:*] - cy
     scy = [-1, where(dcy GT avg(dcy), ny), nph-1]
     simy = intarr(ny+1)
     FOR i=1, ny+1 DO simy[i-1] = round(mean(cy[scy[i-1]+1:scy[i]]))
     simy = simy[where((simy GT 32) AND (simy LT sy-32))]
      
     ;; compute initial shift maps:  first all the PH positions
      
     xx = indgen(sx) # replicate(1, sy)
     yy = replicate(1, sx) # indgen(sy)
     get_lun, unit
      
     FOR im=0, np-1 DO BEGIN
        pic = red_clipim(pics[*, *, im], cl[*, im+1])
        mask = red_separate_mask(pic gt tr * max(pic))
        nph1 = max(mask)
        cc1 = fltarr(2, nph1)
      
        FOR i=0, nph1-1 DO cc1[*, i] = red_com(mask EQ i+1)
      
        cx1 = reform(cc1[0, *])
        cy1 = reform(cc1[1, *])
      
        ;; find matching PH pairs
        nf = 0
        FOR i=0, nph-1 DO BEGIN
           ix = where((abs(cx1-cc[0, i]) LT 20) AND (abs(cy1-cc[1, i]) LT 20))
           IF ix[0] GE 0 THEN BEGIN
              IF nf EQ 0 THEN BEGIN
                 refpos = cc[*, i]
                 xpos = cx1[ix[0]]-cc[0, i]
                 ypos = cy1[ix[0]]-cc[1, i]
              ENDIF ELSE BEGIN
                 refpos = [[refpos], [cc[*, i]]]
                 xpos = [[xpos], [cx1[ix[0]]-cc[0, i]]]
                 ypos = [[ypos], [cy1[ix[0]]-cc[1, i]]]
              ENDELSE
              nf++
           ENDIF
        ENDFOR
         
        ;; fit plane
                                
        if im eq 0 then fil = tfil else fil =rfil 
        fit = sfit([refpos, xpos], 1, /irr, kx=kx, /max)
        xfit = kx[0] + yy*kx[1] + xx*kx[2]
        ;;xfit = red_getplane(kx, xx, yy)
        
        xfit = fix(round(100*temporary(xfit)))
        file_mkdir, self.out_dir+'/calib/'
        fzwrite, xfit,self.out_dir+'/calib/'+ file_basename(fil,'.pinh')+'.xoffs', ' '
         
        fit = sfit([refpos, ypos], 1, /irr, kx=kx, /max)
        yfit = kx[0] + yy*kx[1] + xx*kx[2]

        ;;yfit = red_getplane(kx, xx, yy)

        yfit = fix(round(100*temporary(yfit)))
        fzwrite, yfit, self.out_dir+'/calib/'+file_basename(fil,'.pinh')+'.yoffs', ' '
         
        tfilo = file_basename(wfil,'pinh')
        ;; write out a nice file for pinholecalib.py
        opinh = file_basename(fil,'pinh')
        openw, unit, self.out_dir+'/calib/'+file_basename(fil,'.pinh')+'.cfg'
        printf, unit, 'object{'
        printf, unit, '  WAVELENGTH='+strtrim(lam, 2)
        printf, unit, '  channel{'
        printf, unit, '    FILENAME_TEMPLATE='+tfilo+'%07d'
        printf, unit, '    IMAGE_NUM='+strtrim(im, 2)
        printf, unit, '    DIVERSITY=0.00 mm'
        printf, unit, '    '+acl[0]
        printf, unit, '  }'
        printf, unit, '  channel{'
        printf, unit, '    FILENAME_TEMPLATE='+opinh+'%07d'
        printf, unit, '    IMAGE_NUM='+strtrim(im, 2)
        printf, unit, '    DIVERSITY=0.00 mm'
        printf, unit, '    '+acl[im+1]
        printf, unit, '    XOFFSET='+opinh+'xoffs'
        printf, unit, '    YOFFSET='+opinh+'yoffs'
        printf, unit, '  }'
        printf, unit, '}'
        printf, unit, 'BASIS=Zernike'
        printf, unit, 'GRADIENT=gradient_diff'
        printf, unit, 'GETSTEP=getstep_steepest_descent'
        printf, unit, 'PROG_DATA_DIR=./data/'
        printf, unit, 'MODES=2-3'
        printf, unit, 'NUM_POINTS=64'
        printf, unit, 'TELESCOPE_D=9.700000E-01'
        printf, unit, 'ARCSECPERPIX=5.920000E-02'
        printf, unit, 'PIXELSIZE=1.600000E-05'
        printf, unit, 'SIM_X='+strjoin(strtrim(simx, 2), ',')
        printf, unit, 'SIM_Y='+strjoin(strtrim(simy, 2), ',')
        printf, unit, 'CALIBRATE'
        free_lun, unit
          
        ;; File links for momfbd -> files  need a number!
        nout =  self.out_dir+'/calib/'+file_basename(wfil, 'pinh') + red_stri(im, ni='(I07)')
        if(file_test(nout)) then spawn, 'rm '+nout
        print, inam+' : creating '+file_basename(nout)
        file_link, wfil, nout
        nout =  self.out_dir+'/calib/'+file_basename(fil, 'pinh') + red_stri(im, ni='(I07)')
        if(file_test(nout)) then spawn, 'rm '+nout
        print, inam+' : creating '+file_basename(nout)
        file_link, fil, nout
     endfor

  ENDFOR
  
  return
end
