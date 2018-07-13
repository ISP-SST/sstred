; docformat = 'rst'

;+
; Make quick-look movies and images from raw data.
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;     Mats Löfdahl, ISP
;     (Using code from an earlier version by Jaime.)
; 
; :Returns:
; 
; 
; :Params:
; 
; 
; :Keywords:
;
;    align : in, optional, type=boolean
;
;      Set this to align the cube.
;
;    cam : in, optional, type=string, default="A narrowband camera"
;   
;      Make quicklook for this camera only. 
;   
;    clip  : 
;   
;    maxshift : in, optional, type=integer, default=6
;
;      When aligning, filter out isolated shifts larger than this.
;
;    overwrite  : 
;   
;   
;   
;    x_flip  : 
;   
;   
;   
;    y_flip  : 
;   
;   
;   
; 
; 
; :History:
; 
;   2018-07-10 : MGL. First version.
; 
; 
;-
pro red::quicklook, align = align $
                    , cam = cam $
                    , clip = clip $
                    , dark = dark $
                    , gain =  gain $
                    , maxshift = maxshift $
                    , no_histo_opt = no_histo_opt $
                    , overwrite = overwrite $
                    , pattern = pattern $
                    , remote_dir = remote_dir $
                    , remote_login = remote_login $
                    , ssh_find = ssh_find $
                    , use_state = use_state $
                    , x_flip = x_flip $
                    , y_flip = y_flip 
  
  inam = red_subprogram(/low, calling = inam1)

  if ~keyword_set(cam) then begin
    cindx = where(~strmatch(*self.cameras,'*-[WD]'), Ncam)
    if Ncam eq 0 then begin
      print, inam + ' : No NB camera in self.cameras?'
      print, *self.cameras
      retall
    endif
    ;; If more than one, pick the first
    cam = (*self.cameras)[cindx[0]]
  endif
  
  if ~ptr_valid(self.data_dirs) then begin
    print, inam+' : ERROR : undefined data_dir'
    return
  endif

  if n_elements(maxshift) eq 0 then maxshift = 6

  ;; Select data sets.
  tmp = red_select_subset(file_basename(*self.data_dirs) $
                          , qstring = inam + ' : Select data sets' $
                          , count = Nsets, indx = ichoice)

  dirs = (*self.data_dirs)[ichoice] 
  
  ;; If search pattern is given, use that, otherwise just use *
  if keyword_set(pattern) then pat = "*" + pattern + "*" else pat = "*"
  

  for iset = 0, Nsets-1 do begin
                                ;stop
    searchstring = dirs[iset] + '/' + cam + '/' + pat
    files = red_file_search(searchstring, count = Nfiles)
 
    if files[0] eq '' then begin
      print, inam + ' : ERROR -> no frames found in '+dirs[iset]
      continue
    endif

    self -> extractstates, files, states
    ustat = states[uniq(states.tun_wavelength, sort(states.tun_wavelength))].fullstate

    if n_elements(ustat) gt 1 then begin
      ;; Select states.
      tmp = red_select_subset(ustat $
                              , qstring = inam + ' : Select states' $
                              , count = Nstates, indx = ichoice)
      ustat = ustat[ichoice]
    endif else begin
      Nstates = 1
    endelse
    

    for istate = 0, Nstates-1 do begin
      
      self -> selectfiles, files = files, states = states $
                           , ustat = ustat[istate] $
                           , selected = sel

      uscan = states[uniq(states.scannumber,sort(states.scannumber))].scannumber
      Nscans = n_elements(uscan)

;      sel = sel[sort(states[sel].scannumber)] ; Make sure scans are in order!
      
      ;; Load dark and flat
      self -> get_calib, states[sel[0]], darkdata = dd, gaindata = gg $
                         , darkstatus = darkstatus, gainstatus = gainstatus
      if darkstatus ne 0 then dd = 0.
      if gainstatus ne 0 then gg = 1.
    
      print, inam + ' : found scans -> '+red_stri(Nscans)

      timestamp = file_basename(dirs[iset])

      outdir = self.out_dir +'/quicklook/'+timestamp+'/'
      stop
      namout = states[sel[0]].camera
      namout += '_' + 'quick'
      namout += '_' + self.isodate
      namout += '_' + timestamp
      namout += '_' + states[sel[0]].prefilter $
                + '_' + states[sel[0]].tuning
      namout += '.mp4'

      if ~keyword_set(overwrite) and file_test(outdir+namout) then continue

      
      file_mkdir, outdir
      
;      Ntot = 100. / (Nscans - 1.0)

      print, inam + ' : saving to folder -> '+outdir 

      dim = size(dd, /dim)
      if(keyword_set(clip)) then begin
        x0 = clip[0]
        x1 = clip[1]
        y0 = clip[2]
        y1 = clip[3]
      endif else begin
        ;; Dimensions
        x0 = 0
        x1 = dim[0] - 1
        y0 = 0
        y1 = dim[1] - 1
        ;; Margin
        x0 += 50
        x1 -= 50
        y0 += 10
        y1 -= 10
      endelse

      ;; RGB cube
      cube = fltarr(dim[0], dim[1], Nscans)
      if gainstatus eq 0 then mask = red_cleanmask(gg ne 0)
      best_contrasts = fltarr(Nscans)
      
      for iscan = 0L, Nscans -1 do begin

        red_progressbar, iscan, Nscans, 'Assembling cube for movie '+namout, /predict

        self -> selectfiles, files = files, states = states $
                             , ustat = ustat[istate] $
                             , scan = uscan[iscan] $
                             , selected = sel2

        if n_elements(sel2) eq 1 then begin
          ;; Just a single file
          ims = red_readdata(states[sel[iscan]].filename)
        endif else begin
          ;; Multiple files, e.g. WB. (Ideally, you should be able to
          ;; call red_readdata with multiple file names and if would
          ;; do the right thing. Needs to produce a combined header,
          ;; too!) 
          Nframes = 0L
          for ifile = 0, n_elements(sel2)-1 do begin
            Nframes += fxpar(red_readhead(files[sel2[ifile]]),'NAXIS3')
          endfor                ; ifile
          ims = fltarr(dim[0], dim[1], Nframes)
          iframe = 0
          for ifile = 0, n_elements(sel2)-1 do begin
            ims[0, 0, iframe] = red_readdata(files[sel2[ifile]], head = head)
            iframe += fxpar(head,'NAXIS3')
          endfor                ; ifile
        endelse

        if size(ims, /n_dim) eq 3 then begin
          ;; Improvement: base this selection on wb contrast?
          Nframes = (size(ims, /dim))[2]
          contrasts = fltarr(Nframes)
          for iframe = 0, Nframes-1 do contrasts[iframe] $
             = stddev(ims[x0:x1, y0:y1, iframe])/mean(ims[x0:x1, y0:y1, iframe])
          mx = max(contrasts, ml)
          im = (ims[*, *, ml] - dd) * gg
          best_contrasts[iscan] = contrasts[ml]
        endif else begin
          im = (ims - dd) * gg
          best_contrasts[iscan] = stddev(im[x0:x1, y0:y1])/mean(im[x0:x1, y0:y1])
        endelse
  
        if gainstatus eq 0 then im = red_fillpix(im, mask=mask, nthreads=6)
        
        idx1 = where(im eq 0.0, bla, complement = idx)
        if(bla gt 0L) then im[idx1] = median(im[idx])
        
        if(keyword_set(x_flip)) then im = reverse(temporary(im), 1)
        if(keyword_set(y_flip)) then im = reverse(temporary(im), 2)
        
        cube[*, *, iscan] = im

      endfor                    ; iscan

      if keyword_set(align) then begin
        
        ;; Measure image shifts
        shifts = red_aligncube(cube, 5 $ ;, cubic = -0.5 $
                               , xbd = round(dim[0]*.9), ybd = round(dim[1]*.9), /center )

        ;; Outliers?
        indx_included = where((abs(shifts[0,*] - median(reform(shifts[0,*]),3)) le maxshift) $
                              and (abs(shifts[1,*] - median(reform(shifts[1,*]),3)) le maxshift) $
                              , complement = indx_excluded, Nincluded, Ncomplement = Nexcluded)
        if Nexcluded gt 0 then begin
          shifts[0, indx_excluded] = interpol(shifts[0, indx_included], indx_included, indx_excluded)
          shifts[1, indx_excluded] = interpol(shifts[1, indx_included], indx_included, indx_excluded)
        endif
        
        ;; Align the cube
        for iscan = 0, Nscans-1 do begin
          cube[*, *, iscan] = red_shift_im(cube[*, *, iscan] $
                                           , shifts[0, iscan], shifts[1, iscan] $
                                           , cubic = -0.5)
    
        endfor                  ; iscan
      endif
      
      if ~keyword_set(no_histo_opt) then cube = red_histo_opt(temporary(cube))

      ;; The write_video command takes an rgb cube with dimensions
      ;; [3,Nx,Ny,Nscans] 
      rgbcube = transpose(rebin(bytscl(cube[x0:x1, y0:y1, *]) $
                                , x1-x0+1, y1-y0+1, Nscans, 3, /samp) $
                          , [3,0,1,2])
      
      fps = 8
      write_video, outdir+namout, rgbcube, VIDEO_FPS=fps

      ;; Make a jpeg image of the best frame. (We could make it write
      ;; several, locally best, frames. Say, the best of every group
      ;; of 10 or so? Define keywords for this!)
      mx = max(best_contrasts, ml)
      im = bytscl(cube[x0:x1, y0:y1, ml])
      ;; Add tickmarks here?
      write_jpeg, outdir+red_strreplace(namout, '.mp4', '_scan='+strtrim(uscan[ml], 2)+'.jpg') $
                  , im, q = 100
      
    endfor                      ; istate
  endfor                        ; iset
  
end
