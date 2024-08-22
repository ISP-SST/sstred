; docformat = 'rst'

;+
; Make quicklook images of mosaics.
; 
; :Categories:
;
;    SST pipeline
; 
; :Author:
; 
;    Mats Löfdahl, Institute for Solar Physics
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
;   2024-07-12 : MGL. First version.
; 
;   2024-07-18 : MGL. Works for CHROMIS.
; 
;   2024-08-22 : MGL. New keyword debug.
; 
;-
pro red::quicklook_mosaic, align=align $
                           , cam=cam $
                           , compress_factor = compress_factor $
                           , core_and_wings=core_and_wings $
                           , cube=cube $
                           , dark=dark $
                           , datasets=datasets $
                           , debug = debug $
                           , destretch=destretch $
                           , maxshift=maxshift $
                           , nthreads=nthreads $
                           , prefilters = prefilters $
                           , textcolor=textcolor $
                           , use_states=use_states 

  inam = red_subprogram(/low, calling=inam1)

  if keyword_set(core_and_wings) then undefine, use_states
  
  if ~keyword_set(cam) then begin
    cindx = where(~strmatch(*self.cameras, '*-[WD]'), Ncam)
    if Ncam eq 0 then begin
      print, inam + ' : No NB camera in self.cameras?'
      print, *self.cameras
      retall
    endif
    ;; If more than one, pick the first
    cam = (*self.cameras)[cindx[0]]
  endif
  detector = self -> getdetector(cam)

  if n_elements(compress_factor) eq 0 then fac = 4 else fac = compress_factor
  
  is_wb = strmatch(cam, '*-[WB]')
  
  if ~ptr_valid(self.data_dirs) then begin
    print, inam+' : ERROR : undefined data_dir'
    return
  ENDIF

  if n_elements(textcolor) eq 0 then textcolor = 'yellow'
  if n_elements(maxshift) eq 0 then maxshift = 10
  if n_elements(nthreads) eq 0 then nthreads = 6


  ;; Limit possible datasets to mosaic directories
  mos_dirs = (*self.data_dirs)[where(strmatch(*self.data_dirs, '*mosaic*'), Nwhere)]
  IF Nwhere EQ 0 THEN stop
  
  if n_elements(datasets) eq 0 then begin
    ;; Select data sets in a menu.
    tmp = red_select_subset(file_basename(mos_dirs) $
                            , qstring = inam + ' : Select data sets' $
                            , count = Nsets, indx = ichoice)
    dirs = mos_dirs[ichoice] 
  endif else begin
    if max(datasets eq '*') eq 1 then begin
      ;; If called with datasets='*' (or, if an array, any element is
      ;; '*') then process all datasets.
      dirs = mos_dirs
    endif else begin
      ;; Select datasets that match the ones listed in the datasets
      ;; keyword.
      match2, datasets, file_basename(mos_dirs), suba, subb 
      mindx = where(subb ne -1, Nwhere)
      if Nwhere eq 0 then begin
        print, inam + ' : No data directories match the datasets keyword.'
        return
      endif
      dirs = (mos_dirs)[mindx]
    endelse
    Nsets = n_elements(dirs)
  endelse
  
  ;; Now loop over datasets (timestamps)
  for iset=0, Nsets-1 do begin

    timestamp = file_basename(dirs[iset])

    print, inam + ' : Working on '+timestamp

    files = self -> raw_search(dirs[iset] + '/' + cam, scanno=0, count=Nfiles $
                               , prefilters = prefilters )

    if files[0] eq '' then begin
      print, inam + ' : ERROR -> no frames found in '+dirs[iset]
      continue                  ; Goto next dataset
    endif
    
    umos = reform((stregex(files, 'mos([0-9][0-9])', /extract, /subexpr))[1, *])
    umos = umos[uniq(umos, sort(umos))]
    Nmos = n_elements(umos)
    
    self -> extractstates, files, states

    upref = states(uniq(states.prefilter, sort(states.prefilter))).prefilter
    Npref = n_elements(upref)

    outdir = self.out_dir +'/quicklook_mosaic/'+timestamp+'/'
    file_mkdir, outdir

    instrument = strupcase((strsplit(states[0].camera, '-', /extract))[0])
    
    if is_wb then begin

      Nstates = 1
      
    endif else begin
      
      indx = uniq(states.tun_wavelength, sort(states.tun_wavelength))
      ustat = states[indx].fullstate        
      
      states_count = 0
      undefine, pat, ustat2, sel_in
      case 1 of
        
        keyword_set(core_and_wings) : begin
          
          for ipref=0, Npref-1 do begin
            sindx = where(strmatch(ustat, '*_'+upref[ipref]+'_*'), Nmatch)
            if Nmatch eq 1 then begin
              ;; If just one state for this prefilter, then use it!
              red_append, ustat2, ustat[sindx[0]]
              red_append, sel_in, 0
            endif else begin              
              imatch = where(strmatch(ustat[sindx], '*+0*'), Nmatch)
              if Nmatch gt 1 then begin
                print, 'There is more then one spectral line in a scan with ', upref[ipref], ' prefilter.'
                print, 'You have to choose spectral points manually.'
                undefine, ustat2
                break
              endif
              if Nmatch gt 0 then begin
                red_append, ustat2, ustat[sindx[imatch]]
                red_append, sel_in, sindx[imatch]
                wv = states[indx[sindx]].tun_wavelength
                if n_elements(wv[0:imatch]) mod 2 ne 1 then ex = 1 else ex = 0
                wng = median(wv[0:imatch-ex])
                in = where(wv eq wng)
                red_append, sel_in, sindx[in]
                red_append, ustat2, ustat[sindx[in]]
                if n_elements(wv[imatch:*]) mod 2 ne 1 then ex = 1 else ex = 0
                wng = median(wv[imatch+ex:*])
                in = where(wv eq wng)
                red_append, sel_in, sindx[in]
                red_append, ustat2, ustat[sindx[in]]
              endif else begin
                np = n_elements(sindx)
                red_append, ustat2, ustat[sindx[np/4]]
                red_append, ustat2, ustat[sindx[np/2]]
                red_append, ustat2, ustat[sindx[np/2+np/4]]
                red_append, sel_in, [sindx[np/4], sindx[np/2], sindx[np/2+np/4]]
              endelse
            endelse
          endfor                ; ipref
          
          states_count = n_elements(ustat2)
          for istate=0, states_count-1 do begin
            fn = states[indx[sel_in[istate]]].filename
            prts = strsplit(fn, '[_.]', /extract)
            if strmatch(cam, '*W*') then begin
              if instrument eq 'CHROMIS' then $
                 red_append, pat, '*' + prts[-3] + '*' $
              else $
                 red_append, pat, '*' +prts[-5] + '*'
            endif else begin
              if instrument eq 'CHROMIS' then $
                 red_append, pat, '*' + prts[-3] + '*' + prts[-2] + '*' $
              else $
                 red_append, pat, '*' +prts[-5] + '*' + prts[-4] + '*' + prts[-3] + '*'
            endelse
          endfor                ; istate
          use_states = ustat2
          
        end                     ; core_and_wings
        
        n_elements(use_states) gt 0 : begin
          
          for istate=0, n_elements(use_states)-1 do begin
            imatch = where(strmatch(ustat, '*'+use_states[istate]+'*'), Nmatch)
            if Nmatch ge 1 then begin
              states_count++
              red_append, ustat2, ustat[imatch] 
              fn = states[indx[imatch]].filename
              prts = strsplit(fn, '[_.]', /extract)
              if strmatch(cam, '*W*') then begin
                if instrument eq 'CHROMIS' then $
                   red_append, pat, '*' + prts[-3] + '*' $
                else $
                   red_append, pat, '*' +prts[-5] + '*'
              endif else begin
                if instrument eq 'CHROMIS' then $
                   red_append, pat, '*' + prts[-3] + '*' + prts[-2] + '*' $
                else $
                   red_append, pat, '*' +prts[-5] + '*' + prts[-4] + '*' + prts[-3] + '*'
              endelse
            endif else print, 'There is no match for ', use_states[istate], ' state.'
          endfor                ; istate
          if states_count eq 0 then $
             print, 'There are no matches for provided states. You have to choose states manually.'
          
        end                     ; use_states
        
        else : 
        
      endcase
      
      if states_count eq 0 then begin 
        tmp = red_select_subset(states[indx].prefilter+'_'+states[indx].fullstate $
                                , qstring=inam + ' : Select states' $
                                , count=Nstates, indx=ichoice)
        ustat2 = ustat[ichoice]
        for istate=0, n_elements(ustat2)-1 do begin
          imatch = where(ustat2[istate] eq ustat)
          fn = states[indx[imatch]].filename
          prts = strsplit(fn, '[_.]', /extract)
          if strmatch(cam, '*W*') then begin
            if instrument eq 'CHROMIS' then $
               red_append, pat, '*' + prts[-3] + '*' $
            else $
               red_append, pat, '*' +prts[-5] + '*'
          endif else begin
            if instrument eq 'CHROMIS' then $
               red_append, pat, '*' + prts[-3] + '*' + prts[-2] + '*' $
            else $
               red_append, pat, '*' +prts[-5] + '*' + prts[-4] + '*' + prts[-3] + '*'
          endelse
        endfor                  ; istate
        use_states = ustat2      
      endif
      
      ustat = ustat2
      Nstates = n_elements(ustat)

    endelse                     ; is_wb

    for istate=0, Nstates-1 do begin

      if is_wb then begin

        selfiles = files
        selstates = states
        Nsel = Nfiles

        pref = selstates[0].prefilter

        namout = selstates[0].camera
        namout += '_' + 'quick'
        namout += '_' + self.isodate
        namout += '_' + timestamp
        namout += '_' + pref

        annstring = 'SST ' + self.isodate
        annstring += ' ' + timestamp
        annstring += ' ' + 'UT, '
        annstring += ' ' + selstates[0].camera
        annstring += ' ' + pref + ' $\Angstrom$'
        
      endif else begin
        
        self -> selectfiles, files=files, states=states $
                             , ustat=ustat[istate] $
                             , selected=sel, count=Nsel
        
        if Nsel eq 0 then continue
        
        selfiles = files[sel]
        selstates = states[sel]
        
        pref = selstates[0].prefilter
        
        namout = selstates[0].camera
        namout += '_' + 'quick'
        namout += '_' + self.isodate
        namout += '_' + timestamp
        namout += '_' + pref 
        namout += '_' + selstates[0].tuning

        annstring = 'SST ' + self.isodate
        annstring += ' ' + timestamp
        annstring += ' ' + 'UT, '
        annstring += ' ' + selstates[0].camera
        fpisplit = strsplit(selstates[0].tuning, '_', /extract)
        annstring += ' ' + fpisplit[0] + ' $\Angstrom$' + ' ' + fpisplit[1] + ' m$\Angstrom$'
      endelse
      
      image_scale = float(self -> imagescale(pref, /use_config))

      if ~file_test('calib/alignments.sav') then begin
        ;; Try to make sure we have alignments.

        ;; At this point, we can call sumpinh for one state per
        ;; camera, and then call pinholecalib. Hopefully this provides
        ;; the needed info in alignments.sav, so we restore it and try
        ;; again. Failure is caught in the next case statement.
        
        for icam = 0, n_elements(*self.cameras)-1 do begin
          case 1 of

            strmatch((*self.cameras)[icam], '*-D') : begin
              ;; Don't care about diversity camera D
            end
            
            strmatch((*self.cameras)[icam], '*-W') : begin
              ;; Wideband camera W
              self -> sumpinh, cams = (*self.cameras)[icam] $
                               , pref = selstates[0].prefilter
            end
            
            else : begin
              ;; Narrowband cameras NRT
              self -> sumpinh, cams = (*self.cameras)[icam] $
                               , ustat = selstates[0].fullstate $
                               , pref = selstates[0].prefilter
            end

          endcase 
        endfor                  ; icam
        self -> pinholecalib, /verify, nref=20, max_shift = 30 
      endif
      
      ;;self -> getalignment, align=align, prefilters=pref
      self -> getalignment, align=alignment
      indx = where(alignment.state2.camera eq cam and $
                   alignment.state2.prefilter eq pref $
                   , Nalign)

      help, Nalign

      if Nalign eq 0 then begin
        print, inam + ' : Pinhole calibration not done for prefilter ' + pref
        continue                ; Can't do this from within the "case" statement below
      end
      
      case Nalign of
        1    : amap =      alignment[indx].map
        else : amap = mean(alignment[indx].map, dim=3)
       endcase
      amap_inv = invert(amap)
      amap_inv /= amap_inv[2, 2] ; Normalize

      if keyword_set(debug) then begin
        print, 'Inverse amap:'
        print, amap_inv
      endif
      
      ;; Get date_avg from file headers
      date_avg = red_fitsgetkeyword_multifile(selfiles, 'DATE-AVG', count = Ndates)
      time_avg = red_time2double(strmid(date_avg, 11))

      ;; Get pointing
      red_logdata, self.isodate, time_avg, diskpos=diskpos, rsun = rsun
;      cgwindow
;      cgplot, /add, diskpos[0,*], diskpos[1,*], psym=1, color='green' $
;              , xtitle = 'HPLN / 1"', ytitle = 'HPLT / 1"' $
;              , title = 'Mosaic pointing '+self.isodate+' '+timestamp

      hdr = red_readhead(selfiles[0])
      naxis = fxpar(hdr,'NAXIS*')

      self -> get_calib, selstates[0], darkdata = dd, gaindata = gg
      mask = float(gg gt 0.)

      mmask = mask*0 + 1
      if mask[0, 0] eq 0 then begin
        indx = search2d(mask gt 0, 0,0,0,0)
        mmask[indx] = 0
      endif
      if mask[0,naxis[1]-1] eq 0 then begin
        indx = search2d(mask gt 0, 0,naxis[1]-1,0,0)
        mmask[indx] = 0
      endif
      if mask[naxis[0]-1,0] eq 0 then begin
        indx = search2d(mask gt 0, naxis[0]-1,0,0,0)
        mmask[indx] = 0
      endif
      if mask[naxis[0]-1,naxis[1]-1] eq 0 then begin
        indx = search2d(mask gt 0, naxis[0]-1,naxis[1]-1,0,0)
        mmask[indx] = 0
      endif
      ww = morph_distance(mmask,neigh=3)
      ww /= max(ww)
      weight = sin(ww*!pi/2.)^2

      if total(self.direction eq [1, 3, 4, 6]) eq 1 then naxis[[0, 1]] = naxis[[1, 0]]

      weight = rdx_img_project(amap_inv, weight, /preserve_size) >0
      mask   = rdx_img_project(amap_inv, mask,   /preserve_size)
      
      weight = red_rotate(weight, self.direction)
      mask   = red_rotate(mask,   self.direction)
      
      ;; Find out area after rotation
      ang = red_lp_angles(time_mos, self.isodate, /from_log, offset_angle = self.rotation)
      maxangle = max(abs(ang))
      ff = [maxangle,0,0,0,0,reform(ang)]
      dum = red_rotation(weight, full=ff $
                         , ang[0], 0, 0, background = bg, nthreads=nthreads)
      naxis = size(dum,/dim)
      
      t = dblarr(Nmos)
      xpos = lonarr(Nmos)
      ypos = lonarr(Nmos)
      images  = fltarr(naxis[0], naxis[1], Nmos)
      weights = fltarr(naxis[0], naxis[1], Nmos)
      masks   = fltarr(naxis[0], naxis[1], Nmos)

      for imos=0, Nmos-1 do begin

        indx = where(strpos(selfiles, 'mos'+umos[imos]) ne -1, Nwhere)
        if Nwhere eq 0 then stop

        ims = red_readdata_multiframe(selfiles[indx])
        undefine, date_mos
        for i = 0, Nwhere-1 do begin
          tmp = red_readhead(selfiles[indx[i]], date_beg = dbeg)
          red_append, date_mos, dbeg
        endfor                  ; i

        time_mos = red_time2double(strmid(date_mos, 11))
        red_logdata, self.isodate, time_mos, diskpos=diskpos_mos, rsun = rsun

        ;; Select image with highest contrast
        dims = size(ims, /dim)
        sd = fltarr(dims[2])
        for i = 0, Nwhere-1 do begin
          ims[*, *, i] = (ims[*, *, i]-dd)*gg
          sd[i] = stddev(red_centerpic(ims[*, *, i], sz = 500))
        endfor
        mx = max(sd, loc)
        im = ims[*, *, loc[0]]
        im = rdx_img_project(amap_inv, im, /preserve_size)
        im = red_rotate(im, self.direction)
        
        t[imos] = time_mos[loc[0]]
        pos = diskpos_mos[*, loc[0]]
;        red_logdata, self.isodate, t[imos], pig = pos
        
        if ~finite(pos[0]) then pos[0] = median(diskpos_mos[0,*])
        if ~finite(pos[1]) then pos[1] = median(diskpos_mos[1,*])
        
        ang = red_lp_angles(t[imos], self.isodate, /from_log, offset_angle = self.rotation)

        xpos[imos] = pos[0]/image_scale
        ypos[imos] = pos[1]/image_scale

;        cgplot, /add, /over, xpos[imos]*image_scale, ypos[imos]*image_scale, psym = 16, color = 'red'
        if keyword_set(debug) then print, 'pos:', imos, xpos[imos]*image_scale, ypos[imos]*image_scale

        images[*, *, imos]  = red_rotation(im,     ang[0], background=0.0, full=ff, nthreads=nthreads)
        weights[*, *, imos] = red_rotation(weight, ang[0], background=0.0, full=ff, nthreads=nthreads) >0
        masks[*, *, imos]   = red_rotation(mask,   ang[0], background=0.0, full=ff, nthreads=nthreads) $
                              * (weights[*, *, imos] gt 0)
        
      endfor                    ; imos

      xc = round(mean(xpos*image_scale))
      yc = round(mean(ypos*image_scale))
      annstring += ', HPLN/T: (' + red_stri(xc) + '",' + red_stri(yc) + '")'
      
      xpos = xpos - min(xpos) + naxis[0]/2 + 50
      ypos = ypos - min(ypos) + naxis[1]/2 + 50

      ;; Size of array
      Sx = max(xpos) + naxis[0]/2 + 50
      Sy = max(ypos) + naxis[1]/2 + 50

      if Sx gt Nmos*naxis[0] or Sy gt Nmos*naxis[1] then begin
        print, inam + ' : Failed to get good tile positions'
        stop
        continue
      endif
      

;      cgwindow
;      cgplot, /add, xpos, ypos, psym=16, color='red' $
;              , xrange = [0, Sx], yrange = [0, Sy] $
;              , xtitle = 'x', ytitle = 'y', title = 'Mosaic pointing '+self.isodate+' '+timestamp
      
;      tiles = fltarr(round(1.5*(max(xpos) + Sx/2. + 20.)/fac), round(1.5*(max(ypos) + Sy/2. + 20.)/fac))
      
      imap = fltarr(Sx/fac, Sy/fac, Nmos) ; image map
      wmap = fltarr(Sx/fac, Sy/fac, Nmos) ; weight map
      mmap = fltarr(Sx/fac, Sy/fac, Nmos) ; mask map
      
      for imos=0, Nmos-1 do begin

        red_progressbar, imos, Nmos, 'Mosaic tile '+red_stri(imos)

        llx = xpos[imos]/fac - naxis[0]/2/fac
        lly = ypos[imos]/fac - naxis[1]/2/fac
        urx = llx + naxis[0]/fac - 1
        ury = lly + naxis[1]/fac - 1

        llx2 = 2*(xpos[imos]/fac - naxis[0]/2/fac)
        lly2 = 2*(ypos[imos]/fac - naxis[1]/2/fac)
        urx2 = llx2 + naxis[0]/fac - 1
        ury2 = lly2 + naxis[1]/fac - 1


        if fac eq 1 then begin
;          tiles[llx2 : urx2, lly2 : ury2] = images[*, *, imos]
          imap[llx : urx, lly : ury, imos] = masks[*, *, imos] * images[*, *, imos]
          wmap[llx : urx, lly : ury, imos] = masks[*, *, imos] * weights[*, *, imos]
          mmap[llx : urx, lly : ury, imos] = masks[*, *, imos]
        endif else begin
;          tiles[llx2 : urx2, lly2 : ury2] = congrid( images[*, *, imos] $
;          , naxis[0]/fac, naxis[1]/fac, cubic = -0.5 )
          imap[llx : urx, lly : ury, imos] = congrid( masks[*, *, imos] * images[*, *, imos] $
                                                      , naxis[0]/fac, naxis[1]/fac, cubic = -0.5 )
          wmap[llx : urx, lly : ury, imos] = congrid( masks[*, *, imos] * weights[*, *, imos] $
                                                      , naxis[0]/fac, naxis[1]/fac, cubic = -0.5 ) >0
          mmap[llx : urx, lly : ury, imos] = congrid( masks[*, *, imos] $
                                                      , naxis[0]/fac, naxis[1]/fac, cubic = -0.5 ) $
                                              * (wmap[llx : urx, lly : ury, imos] gt 0)
        endelse
    
;        red_show, tiles, /re
        
      endfor                    ; imos

      cgwindow
      colors = distinct_colors(n_colors = Nmos, /num)
;      stop

      xcp = median(diskpos[0, *]) - Sx*image_scale/2.
      ycp = median(diskpos[1, *]) - Sy*image_scale/2.
      
      cgplot, /add, /nodata, [0], [0] $
              , xrange = xcp+[0, Sx*image_scale] $
              , yrange = ycp+[0, Sy*image_scale] $
              ;;, xrange = [0, Sx*image_scale], yrange = [0, Sy*image_scale] $
              , aspect = Sy/float(Sx) $
              , xtitle = 'HPLN / 1"', ytitle = 'HPLT / 1"' $
              , title = 'Mosaic FOVs '+self.isodate+' '+timestamp+' '+pref

      if keyword_set(debug) then cgplot, /add, /over, diskpos[0,*], diskpos[1,*], psym=16, color='red'
      
      if keyword_set(align) || keyword_set(destretch) then begin

        if 0 then begin

          ;; Start with center tile, then continue outward while
          ;; aligning to sum so far.

          xc = median(xpos)
          yc = median(ypos)
          dr = sqrt((xpos-xc)^2 + (ypos-yc)^2)
          ord = sort(dr)        ; The order to process

          totim = imap[*, *, ord[0]] * wmap[*, *, ord[0]]
          totw  = wmap[*, *, ord[0]]
          totm  = mmap[*, *, ord[0]]
                                ;         stop
          
          for imos = 1, Nmos-1 do begin
            
            commonmask = mmap[*, *, ord[imos]] * totm
            commonarea = total(commonmask)
            if keyword_set(debug) then print, imos, ord[imos], commonarea
            
            if 0 then begin
              im  = commonmask * (totim / (1e-5 + totw))
              imi = commonmask * imap[*, *, ord[imos]]
              red_missing, im,  /inplace, missing_type_wanted = 'median'
              red_missing, imi, /inplace, missing_type_wanted = 'median'
              shifts = red_alignoffset(im,imi)
            endif else begin
              bb = red_boundingbox(commonmask)
              if keyword_set(debug) then print, 'boundingbox : ', bb
              shifts = red_shc_mask((totim / (1e-5 + totw))[bb[0]:bb[2],bb[1]:bb[3]] $
                                     , imap[bb[0]:bb[2],bb[1]:bb[3], ord[imos]] $
                                     , commonmask[bb[0]:bb[2],bb[1]:bb[3]] $
                                     , range = 25, poly = 2)
            endelse
            if keyword_set(debug) then print, shifts
              
            ;; Accept shifts only if they are less than a few arcsec.
            dr = sqrt(total(shifts^2)) * fac * image_scale
            if dr lt maxshift then begin
              imap[*, *, ord[imos]] = red_shift_im(imap[*, *, ord[imos]], shifts[0], shifts[1], missing = 0.0)
              wmap[*, *, ord[imos]] = red_shift_im(wmap[*, *, ord[imos]], shifts[0], shifts[1], missing = 0.0) >0
              mmap[*, *, ord[imos]] = red_shift_im(mmap[*, *, ord[imos]], shifts[0], shifts[1], missing = 0.0) $
                                      * (wmap[*, *, ord[imos]] gt 0)
            endif else stop
            
            totim += imap[*, *, ord[imos]] * wmap[*, *, ord[imos]]
            totw  += wmap[*, *, ord[imos]]
            totm  OR= mmap[*, *, ord[imos]]

          endfor                ; imos

          indx = where(totm gt 0)
          mosaic = fltarr(Sx/fac, Sy/fac)
          mosaic[indx] = totim[indx]/totw[indx]
          red_missing, mosaic, /inplace, missing_type_wanted = 'nan'
          
        endif else begin

          ;; Start with the center tile, then base the order of
          ;; stitching on the largest overlapping area.

          
          xc = median(xpos)
          yc = median(ypos)
          dr = sqrt((xpos-xc)^2 + (ypos-yc)^2)
          ord = sort(dr)      
          red_progressbar, 0, Nmos, 'Add tile '+red_stri(ord[0])+' to the mosaic'
          
          edge = sobel(mmap[*, *, ord[0]])
          indx = array_indices(edge, where(edge gt 0.5*max(edge)))
          cgplot, /add, /over $
                  , xcp+indx[0,*]*fac*image_scale $
                  , ycp+indx[1,*]*fac*image_scale $
                  , psym=3, color=colors[ord[0]]
          cgtext, /add, align = 0.5, color = colors[ord[0]], charsize = 3 $
                  , xcp+median(indx[0,*]*fac*image_scale) $
                  , ycp+median(indx[1,*]*fac*image_scale) $
                  , red_stri(ord[0]) + '('+red_stri(0)+')' 
                  
          
          ;; Add center tile
          totim = imap[*, *, ord[0]] * wmap[*, *, ord[0]]
          totw  = wmap[*, *, ord[0]]
          totm  = mmap[*, *, ord[0]]

          commonarea = fltarr(Nmos)

          ;; Loop over remaining tiles
          for imos = 1, Nmos-1 do begin

            if imos ne Nmos-1 then begin
              ;; Find the tile with the most overlap with the mosaic so
              ;; far assembled
              if keyword_set(debug) then print, imos, ord[imos]
              for jmos = imos, Nmos-1 do begin
                commonmask = mmap[*, *, ord[jmos]] * totm
                commonarea[jmos] = total(commonmask)
                if keyword_set(debug) then print, jmos, ord[jmos], commonarea[jmos]
              endfor            ; jmos
              if keyword_set(debug) then print, ord[imos:Nmos-1]
              if keyword_set(debug) then print, commonarea[imos:Nmos-1]
              ;; Find max
              mxx = max(commonarea[imos:Nmos-1], maxloc)
              
              ;; Swap to get it first of remaining
              tmp = ord[imos]
              ord[imos] = ord[maxloc+imos]
              ord[maxloc+imos] = tmp
              
            endif 

            red_progressbar, imos, Nmos, 'Add tile '+red_stri(ord[imos])+' to the mosaic'
            
            ;; We should now add tile ord[imos] to the mosaic
            edge = sobel(mmap[*, *, ord[imos]])
            indx = array_indices(edge, where(edge gt 0.5*max(edge)))
            cgplot, /add, /over $
                    , xcp+indx[0,*]*fac*image_scale $
                    , ycp+indx[1,*]*fac*image_scale $
                    , psym=3, color = colors[ord[imos]]
            cgtext, /add , align = 0.5, color = colors[ord[imos]], charsize = 3 $
                    , xcp+median(indx[0,*]*fac*image_scale) $
                    , ycp+median(indx[1,*]*fac*image_scale) $
                    , red_stri(ord[imos])  + '('+red_stri(imos)+')' 
                   

            commonmask = mmap[*, *, ord[imos]] * totm

            bb = red_boundingbox(commonmask)
            if keyword_set(debug) then print, 'boundingbox : ', bb
            shifts = red_shc_mask((totim / (1e-5 + totw))[bb[0]:bb[2],bb[1]:bb[3]] $
                                  , imap[bb[0]:bb[2],bb[1]:bb[3], ord[imos]] $
                                  , commonmask[bb[0]:bb[2],bb[1]:bb[3]] $
                                  , range = 25, poly = 2)
            
            if keyword_set(debug) then print, shifts
            
            ;; Accept shifts only if they are less than a few arcsec.
            dr = sqrt(total(shifts^2)) * fac * image_scale
            if dr lt 10 then begin
              imap[*, *, ord[imos]] = red_shift_im(imap[*, *, ord[imos]], shifts[0], shifts[1], missing = 0.0)
              wmap[*, *, ord[imos]] = red_shift_im(wmap[*, *, ord[imos]], shifts[0], shifts[1], missing = 0.0) >0
              mmap[*, *, ord[imos]] = red_shift_im(mmap[*, *, ord[imos]], shifts[0], shifts[1], missing = 0.0) $
                                      * (wmap[*, *, ord[imos]] gt 0)
            endif else stop
            
            totim += imap[*, *, ord[imos]] * wmap[*, *, ord[imos]]
            totw  += wmap[*, *, ord[imos]]
            totm  OR= mmap[*, *, ord[imos]]

          endfor                ; imos

          indx = where(totm gt 0)
          mosaic = fltarr(Sx/fac, Sy/fac)
          mosaic[indx] = totim[indx]/totw[indx]
          red_missing, mosaic, /inplace, missing_type_wanted = 'nan'
          
        endelse
        
      endif else begin

        totim = total(imap*wmap, 3)
        totw = total(wmap, 3)
        totm = total(mmap, 3)
        indx = where(totm gt 0)
        mosaic = fltarr(Sx/fac, Sy/fac)
        mosaic[indx] = totim[indx]/totw[indx]
        
      endelse

      bb = red_boundingbox(mosaic gt 0)
      mosaic = mosaic[bb[0]:bb[2],bb[1]:bb[3]] >0

      red_missing, mosaic, missing_type_wanted = 'nan',/inplace
      mosaic = red_fillpix(mosaic, nthreads=nthreads)
      red_missing, mosaic, missing_value = 0,/inplace

      indx = where(mosaic gt median(mosaic)/100)
      tmp = red_histo_opt(mosaic[indx], cmin = cmin, cmax = cmax)
      
      mosaic = bytscl(mosaic >cmin <cmax)
      
      red_show, mosaic, /scroll
      
      ;; Add annotations and make the rgb cube
      thisDevice = !D.Name
      set_plot, 'Z'
      device, Set_Resolution = [bb[2]-bb[0]+1, bb[3]-bb[1]+1] $
              , Z_Buffer = 0 $
              , Set_Pixel_Depth = 24 $
              , Decomposed=0 $
              , set_font="Helvetica" $
              , /tt_font

      old_font = !P.font
      !P.font = 1
      erase
                                ;tv, mosaic
      cgimage, /keep, red_tickbox(mosaic, image_scale*fac, marg=-sx/200., /arc5, /arc10)
      
      charsize = 1.5
      charsize = (bb[3]-bb[1]+1) / 500.
      cgtext, [0.02], [0.02], annstring $
              , /normal, charsize=charsize, color=textcolor, font=1
      snap2 = tvrd(/true)
      set_plot,'X'
      !P.font = old_font
      
      write_png, outdir + namout + '.png', snap2

      print, inam + ' : Wrote mosaic to ' + outdir + namout + '.png'

    endfor                      ; istate
  endfor                        ; iset

  
end

debug = 0

;; Run quicklook_mosaic only after summing darks and flats. Pinhole
;; calibration is now automatic but do it properly with /verify if the
;; results do not look good.

;; Note: Remove automatically summed calibration data, gaintables, and
;; calib/alignments.sav after making quicklooks. Then do those steps
;; properly with commands from doit.pro.


nthreads=20
cd, '/scratch/mats/2024-07-09'
;cd, '/scratch/mats/2024-06-11'
;cd, '/scratch/mats/2023-10-17'
if 0 then begin
  cd, 'CHROMIS'
  a = chromisred("config.txt", /dev, /no)
  a -> quicklook_mosaic, debug = debug, nthreads = nthreads, cam = 'Chromis-N', /align, compress = 2, /core_and_wings ; use = '3999_+0'
;  a -> quicklook_mosaic, nthreads = nthreads, cam = 'Chromis-N', /align, /core, compress = 4
;  a -> quicklook_mosaic, nthreads = nthreads, cam = 'Chromis-W', /align
endif else begin
  cd, 'CRISP'
  a = crisp2red("config.txt", /dev, /no)
  a -> quicklook_mosaic, debug = debug, nthreads = nthreads, cam = 'Crisp-R', /align, /core_and_wings, compress = 4
;  a -> quicklook_mosaic, nthreads = nthreads, cam = 'Crisp-W', /align, pref = '8542'
;  a -> quicklook_mosaic, nthreads = nthreads, cam = 'Crisp-W', /align, pref = '6563'
;  a -> quicklook_mosaic, nthreads = nthreads, cam = 'Crisp-W', /align, pref = '6173'
endelse




end
