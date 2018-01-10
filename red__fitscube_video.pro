; docformat = 'rst'

;+
; Make a video file from data in a fitscube.
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, Institute for Solar Physics
; 
; 
; :Params:
; 
;   infile : in, type=string
; 
;     The fitscube file.
; 
; 
; :Keywords:
; 
;    clip : in, optional, type=integer, default=2
;
;      A number between 0 and 50 that indicates the percentage of
;      pixels to clip off either end of the cube histogram before
;      performing a linear stretch.
;
;    fps : in, optional, type=integer
;
;      Frames per second.
;
;    gamma : in, optional, type="float or fltarr(3)"
;
;      Log scale gamma, scalar or rgb.
;
;    golden : in, optional, type=boolean
;
;      Set to display in a golden color scheme.
;
;    iwave : in, optional, type=integer, default=0
;
;
;
;    istokes : in, optional, type=integer, default=0
;
;
;
;    outfile : in, optional, type=string, default='same as infile but with mp4 extension'
;
;       
;
;    tickmarks : in, optional, type=boolean
;   
;       Set this to draw arcsec tickmarks.
; 
; :History:
; 
;    2018-01-10 : MGL. First version.
; 
; 
;-
pro red::fitscube_video, infile $
                         , clip = clip $
                         , fps = fps $
                         , gamma = gamma $
                         , golden = golden $
                         , iwave = iwave $
                         , istokes = istokes $
                         , outfile = outfile $
                         , tickmarks = tickmarks
  
  inam = red_subprogram(/low, calling = inam1)

  ;; Check the input file
  if ~file_test(infile) then begin
    print, inam + ' : Input file does not exist: '+infile
    retall
  endif
  
  if n_elements(outfile) eq 0 then outfile = file_dirname(infile) + file_basename(infile, '.fits') + '.mp4'

  hdr = headfits(infile)
  dims = fxpar(hdr, 'NAXIS*')
  Ndims = n_elements(dims)

  if Ndims lt 3 then stop

  Nx = dims[0]
  Ny = dims[1]
  Nframes = dims[4]
  
  if n_elements(iwave) eq 0 then iwave = 0
  if n_elements(istokes) eq 0 then istokes = 0
  if n_elements(fps) eq 0 then fps = 5
  if n_elements(clip) eq 0 then clip = 2

  if keyword_set(golden) then gamma = [0.7, 1.2, 7.0]
  if n_elements(gamma) eq 1 then gamma = [gamma, gamma, gamma]
  
  vidcube = dblarr([1, Nx, Ny, Nframes])

  for iframe = 0, Nframes-1 do begin
    self -> fitscube_getframe, infile, frame, ituning = iwave, istokes = istokes, iscan = iframe
    vidcube[0, *, *, iframe] = frame
  endfor                        ; iframe

  if clip ne 0 then begin
    ;; Code to do histogram clip, stolen from cgclipscl.pro
    maxr = Max(vidcube, MIN=minr, /NAN)
    range = maxr - minr
    binsize = range / 1000.
    h = Histogram(vidcube, BINSIZE=binsize, OMIN=omin, OMAX=omax, /NAN)
    n = N_Elements(vidcube)
    cumTotal = Total(h, /CUMULATIVE)
    minIndex = Value_Locate(cumTotal, n * (clip/100.))
    if minIndex eq -1 then minIndex = 0
    while cumTotal[minIndex] eq cumTotal[minIndex + 1] do minIndex = minIndex + 1
    minThresh = minIndex * binsize + omin
    ;; Not all files can be clipped appropriately.
    maxIndex  = Value_Locate(cumTotal, n * ((100-clip)/100.))
    if (maxIndex eq -1) || (maxIndex eq N_Elements(cumTotal)) || (maxIndex eq minIndex) then stop
    ;; If you are still here, try to clip the histogram.
    while cumTotal[maxIndex] eq cumTotal[maxIndex - 1] do maxIndex = maxIndex - 1
    maxThresh = maxIndex * binsize + omin

    ;; Do the clipping
    vidcube = vidcube > minThresh < maxThresh
  endif
  
  vidcube -= min(vidcube)
  vidcube /= max(vidcube)

  ;; Do the RGB stuff
  vidcube = rebin(vidcube,[3, Nx, Ny, Nframes],/samp)
  if n_elements(gamma) ne 0 then begin
    vidcube[0, *, *, *] = vidcube[0, *, *, *]^gamma[0]
    vidcube[1, *, *, *] = vidcube[1, *, *, *]^gamma[1]
    vidcube[2, *, *, *] = vidcube[2, *, *, *]^gamma[2]
  endif 

  ;; Should be byte scaled
  vidcube *= 255d
  vidcube = byte(round(vidcube >0 <255))

  ;; Tickmarks?
  if keyword_set(tickmarks) then begin
    ;; Code inspired by Luc's tickbox_f.ana
    ma = 0B
    marg = min(dims[0:1])/20 ;30
    hmarg = 0
    sc = float(self.image_scale)
    pix = 1./sc
    Lx = fix(Nx*sc)
    Ly = fix(Ny*sc)
    thickness = 3
    ;; Draw the tickmarks, identical for all rgb colors and all frames.
    for i = 1B, Lx do begin     ; Don't draw first tick mark!
      fac = 0.4 * (2 + ((i mod 5) eq 0) + ((i mod 10) eq 0))/2.
      dd = fix(marg/2 * fac)
      x = round(hmarg+thickness-1+i*pix)
      ;; Along horizontal bottom
      vidcube[*, ((x-thickness/2) >0):((x+thickness/2) <(Nx-1)), 0:dd, *] = ma
      ;; Along horizontal top
      vidcube[*, ((x-thickness/2) >0):((x+thickness/2) <(Nx-1)), ((Ny+1-dd) >0):Ny-1, *] = ma
    endfor                      ; i
    for j = 1B, Ly do begin
      fac = 0.4 * (2 + ((j mod 5) eq 0) + ((j mod 10) eq 0))/2.
      dd = fix(marg/2 * fac)
      y = round(hmarg+thickness-1+j*pix)
      ;; Along vertical left
      vidcube[*, 0:dd, ((y-thickness/2) >0):((y+thickness/2) <(Ny-1)), *] = ma
      ;; Along vertical right
      vidcube[*, ((Nx-1-dd) >0):(Nx-1), ((y-thickness/2) >0):((y+thickness/2) <(Ny-1)), *] = ma
    endfor                      ; j
  endif

  ;; Write the video
  write_video, outfile, vidcube, video_fps = fps
  
end

