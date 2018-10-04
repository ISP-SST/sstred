; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
; 
; 
; 
; :Returns:
; 
; 
; :Params:
; 
; 
; :Keywords:
; 
;    mmt  : 
;   
;   
;   
;    mmr  : 
;   
;   
;   
;    filter  : 
;   
;   
;   
;    destretch  : 
;   
;   
;   
;    dir  : 
;   
;   
;   
;    square  : 
;   
;   
;   
; 
; 
; :History:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2015-03-31 : MGL. Use red_download to get turret log file.
;
;   2016-08-23 : THI. Rename camtag to detector and channel to camera,
;                so the names match those of the corresponding SolarNet
;                keywords.
; 
; 
;-
function red::polarim, dir $
                       , destretch = destretch $
                       , dir = dir $
                       , filter = filter $
                       , mmr = mmr $
                       , newflats = newflats $
                       , no_ccdtabs = no_ccdtabs $
                       , smooth = smooth $
                       , square = square $
                       , mmt = mmt
  
  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])
  
  Nelements = 16                ; Number of elements in modulation matrix

  ;; Camera/detector identification
  self->getdetectors
  wbindx     = where(strmatch(*self.cameras,'Crisp-W'))
  wbcamera   = (*self.cameras)[wbindx[0]]
  wbdetector = (*self.detectors)[wbindx[0]]
  nbtindx     = where(strmatch(*self.cameras,'Crisp-T')) 
  nbtcamera   = (*self.cameras)[nbtindx[0]]
  nbtdetector = (*self.detectors)[nbtindx[0]]
  nbrindx     = where(strmatch(*self.cameras,'Crisp-R')) 
  nbrcamera   = (*self.cameras)[nbrindx[0]]
  nbrdetector = (*self.detectors)[nbrindx[0]]


  
;  ;; Search for folders with reduced data
;  if ~keyword_set(dir) then begin
;    dir = file_search(self.out_dir + '/momfbd/*', /test_dir, count = Ndir)
;
;    if Ndir eq 0 then begin
;      print, inam + ' : No directories found in '+self.out_dir+'/momfbd'
;      return, 0
;    endif
;
;    idx = 0L
;    if Ndir gt 1 then begin
;      print, inam + ' : Found '+red_stri(Ndir)+' sub-folders:'
;      for idir = 0, Ndir - 1 do print, red_stri(idir,ni='(I2)')+' -> '+dir[idir]
;      read,idx,prompt = inam+' : Please select state ID: '
;      dir = dir[idx]
;    endif
;
;    dir = file_search(dir + '/*', /test_dir, count = Ndir)
;    idx = 0L
;    if Ndir gt 1 then begin
;      print, inam + ' : Found '+red_stri(Ndir)+' sub-folders:'
;      for idir = 0, Ndir - 1 do begin
;        print, red_stri(idir,ni='(I2)')+' -> '+file_basename(dir[idir])
;      endfor                    ; idir
;      read,idx,prompt = inam+' : Please select state ID: '
;      dir = dir[idx]
;    endif
;
;    dir+= '/cfg/results/' 
;    print, inam + ' : Processing state -> '+dir
;  endif
;  self -> getdetectors, dir = self.data_dir

  ;; Get files (right now, only momfbd is supported)
  filetype = 'momfbd'
  if self.filetype eq 'ANA' then stop ;filetype = 'f0'
  tfiles = file_search(dir+'/'+nbtdetector+'.*.'+filetype, count = Nimt)
  rfiles = file_search(dir+'/'+nbrdetector+'.*.'+filetype, count = Nimr)

  if Nimt ne Nimr then begin

    print, inam + ' : WARNING, different number of images found for each camera:'
    print, nbtcamera+' '+nbtdetector+' -> '+red_stri(Nimt)
    print, nbrcamera+' '+nbrdetector+' -> '+red_stri(Nimr)

  endif

  ;; Get states that are common to both cameras (object)
  pol = red_getstates_polarim(tfiles, rfiles, self.out_dir $
                              , camt = nbtdetector $
                              , camr = nbrdetector $
                              , camwb = wbdetector $
                              , newflats = newflats)
  Nstates = n_elements(pol)
  
  pref = pol[0]->getvar(7)
  clipfile = self.out_dir+'/calib/align_clips.'+pref+'.sav'
  if file_test(clipfile) then begin
    restore, clipfile
    tclip = cl[*,1]
    rclip = cl[*,2]
  endif
  
  ;; Modulations matrices

  ;; T-Cam
  if ~keyword_set(mmt) then begin
    
    search = self.out_dir+'/polcal/'+nbtdetector+'_'+pol[0]->getvar(7)+'_polcal.fits'
    if file_test(search) then begin
      immt = (f0(search))[0:Nelements-1,*,*]

      ;; interpolate CCD tabs?
      if ~keyword_set(no_ccdtabs) then begin
        ;; Only if Sarnoff cameras! -------------------------------------------------------
        for ii = 0, Nelements-1 do immt[ii,*,*] = red_mask_ccd_tabs(reform(immt[ii,*,*]))
      endif

      ;; Check NaNs
      for ii = 0, Nelements-1 do begin
        mask = 1B -( ~finite(reform(immt[ii,*,*])))
        idx = where(mask, count)
        if count gt 0 then immt[ii,*,*] = red_fillpix(reform(immt[ii,*,*]), mask=mask)
      endfor                    ; ii

      if n_elements(smooth) gt 0 then begin
        dpix = round(smooth)*3
        if (dpix/2)*2 eq dpix then dpix -= 1
        dpsf = double(smooth)
        psf = red_get_psf(dpix, dpix, dpsf, dpsf)
        psf /= total(psf, /double)
        for ii=0,Nelements-1 do immt[ii,*,*] = red_convolve(reform(immt[ii,*,*]), psf)
      endif

      if n_elements(tclip) eq 4 then begin
        Nx = abs(tclip[1] - tclip[0]) + 1
        Ny = abs(tclip[3] - tclip[2]) + 1
        
        print,'Clipping transmitted modulation matrix to '+red_stri(Nx)+' x '+red_stri(Ny)
        
        tmp = make_array( Nelements, Nx, Ny, type=size(immt, /type) )
        for ii=0,Nelements-1 do begin
          tmp[ii,*,*] = red_clipim( reform(immt[ii,*,*]), tclip )
        endfor
        immt = temporary(tmp)
      endif

      immt = ptr_new(red_invert_mmatrix(temporary(immt)))
      
    endif else begin
      print, inam + ' : ERROR, polcal data not found in ' + self.out_dir + '/polcal/'
      return, 0
    endelse

  endif else begin
    immt = red_invert_mmatrix(temporary(mmt))
  endelse

  ;; R-Cam
  if ~keyword_set(mmr) then begin
                                ;
    search = self.out_dir+'/polcal/'+nbrdetector+'_'+pol[0]->getvar(7)+'_polcal.fits'
    if file_test(search) then begin
      immr = (f0(search))[0:Nelements-1,*,*]

      ;; Interpolate CCD tabs?
      if ~keyword_set(no_ccdtabs) then begin
        for ii = 0, Nelements-1 do immr[ii,*,*] = red_mask_ccd_tabs(reform(immr[ii,*,*]))
      endif


      ;; Check NaNs
      for ii = 0, Nelements-1 do begin
        mask = 1B - ( ~finite(reform(immr[ii,*,*])))
        idx = where(mask, count)
        if count gt 0 then immr[ii,*,*] = red_fillpix(reform(immr[ii,*,*]), mask=mask)
      endfor                    ; ii

      ;; Smooth?
      if n_elements(smooth) gt 0 then begin
        dpix = round(smooth)*3
        if (dpix/2)*2 eq dpix then dpix -= 1
        dpsf = double(smooth)
        psf = red_get_psf(dpix, dpix, dpsf, dpsf)
        psf /= total(psf, /double)
        for ii=0,Nelements-1 do immr[ii,*,*] = red_convolve(reform(immr[ii,*,*]), psf)
      endif
      
      if n_elements(rclip) eq 4 then begin
        Nx = abs(rclip[1] - rclip[0]) + 1
        Ny = abs(rclip[3] - rclip[2]) + 1
        
        print,'Clipping reflected modulation matrix to '+red_stri(Nx)+' x '+red_stri(Ny)
        
        tmp = make_array( Nelements, Nx, Ny, type=size(immr, /type) )
        for ii=0,Nelements-1 do begin
          tmp[ii,*,*] = red_clipim( reform(immr[ii,*,*]), rclip )
        endfor                  ; ii
        immr = temporary(tmp)
      endif

      immr = ptr_new(red_invert_mmatrix(temporary(immr)))
    endif else begin
      print, inam + ' : ERROR, polcal data not found in ' + self.out_dir+'/polcal/'
      return, 0
    endelse
                                ;
  endif else begin
    immr = red_invert_mmatrix(temporary(mmr))
  endelse

  ;; Pointers to immt and immr
  for istat = 0L, Nstates - 1 do begin
    pol[istat]->setvar, 5, value = ptr_new(immt)
    pol[istat]->setvar, 6, value = ptr_new(immr)
  endfor                        ; istat

  ;; Fill border information (based on 1st image)
  print, inam + ' : Reading file -> ' + (pol[0]->getvar(9))[0]
  if self.filetype eq 'MOMFBD' then begin
    tmp = red_mozaic(momfbd_read((pol[0]->getvar(9))[0]))
  endif else begin
    tmp = f0((pol[0]->getvar(9))[0])
  endelse
  
  dimim = red_getborder(tmp, x0, x1, y0, y1, square=square)
  for istat = 0L, Nstates - 1 do pol[istat]->fillclip, x0, x1, y0, y1

  ;; Telog
  red_download, date = self.isodate, /turret, pathturret = turretfile
  if turretfile then begin
    print, inam + ' : Using SST position LOG -> ' + turretfile
    for istat = 0L, Nstates - 1 do pol[istat]->setvar, 18, value = turretfile
  endif else begin
    print, inam + ' : No Turret log file'
    stop
  endelse 

  ;; Print states
  print, inam + ' : Found '+red_stri(Nstates)+' state(s) to demodulate:'
  for istat = 0L, Nstat-1 do begin
    print, red_stri(istat, ni='(I5)') + ' -> ' + pol[istat] -> state()
  endfor                        ; istat
  
  return, pol

end
