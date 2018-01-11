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
;    bit_rate : in, optional, type=integer
;
;      The target output data rate, in bits per second. See IDL
;      documentation for write_video and
;      IDLffVideoWrite::AddVideoStream.
;
;    clip : in, optional, type=integer, default=2
;
;      A number between 0 and 50 that indicates the percentage of
;      pixels to clip off either end of the cube histogram before
;      performing a linear stretch.
;
;    fps : in, optional, type=integer, default=5
;
;      Frames per second.
;
;    gamma : in, optional, type=float
;
;      Log scale gamma.
;
;    golden : in, optional, type=boolean
;
;      Set rgbgamma to produce a golden color scheme.
;
;    iwave : in, optional, type=integer, default=0
;
;      The index of the wavelength to be shown in the movie. 
;
;    istokes : in, optional, type=integer, default=0
;
;      The index of the Stokes component to be shown in the movie.
;
;    outfile : in, optional, type=string, default='same as infile but with mp4 extension'
;
;      The path to the output file.       
;                    
;    rgbgamma : in, optional, type=fltarr(3), default="[1,1,1]"
;
;      Gamma values for the RGB channels. The total gamma correction
;      is of the form data(channel)^(gamma*rgbgamma(channel)).    
;                    
;    rgbfac : in, optional, type=fltarr(3), default="[1,1,1]"
;
;      Multiplicative factors for the RGB channels.
;
;    tickcolor : in, optional, type=byte, default=0
;   
;      The color of the tickmarks, 0=black, 255=white.   
;
;    tickmarks : in, optional, type=boolean
;   
;      Set this to draw arcsec tickmarks.
; 
; :History:
; 
;    2018-01-10 : MGL. First version.
; 
; 
;-
pro red::fitscube_video, infile $
                         , bit_rate = bit_rate $
                         , clip = clip $
                         , fps = fps $
                         , gamma = gamma $
                         , golden = golden $
                         , iwave = iwave $
                         , istokes = istokes $
                         , outfile = outfile $
                         , rgbfac = rgbfac $
                         , rgbgamma = rgbgamma $
                         , tickcolor = tickcolor $
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
  if n_elements(tickcolor) eq 0 then tickcolor = 0B else tickcolor = byte(tickcolor)

  ;; Gamma and color
  if n_elements(rgbfac) eq 0 then rgbfac = [1., 1., 1.]
  if keyword_set(golden) then rgbgamma = [0.6, 1.0, 6.0]  
  if n_elements(gamma) eq 0 then gamma = 1.0
  if n_elements(rgbgamma) eq 0 then rgbgamma = [1., 1., 1.]
  totalgamma = gamma*rgbgamma
  
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
  vidcube[0, *, *, *] = rgbfac[0]*vidcube[0, *, *, *]^totalgamma[0]
  vidcube[1, *, *, *] = rgbfac[1]*vidcube[1, *, *, *]^totalgamma[1]
  vidcube[2, *, *, *] = rgbfac[2]*vidcube[2, *, *, *]^totalgamma[2]
   
  ;; Should be byte scaled
  vidcube *= 255d
  vidcube = byte(round(vidcube >0 <255))

  ;; Tickmarks?
  if keyword_set(tickmarks) then begin
    ;; Code inspired by Luc's tickbox_f.ana
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
      vidcube[*, ((x-thickness/2) >0):((x+thickness/2) <(Nx-1)), 0:dd, *] = tickcolor
      ;; Along horizontal top
      vidcube[*, ((x-thickness/2) >0):((x+thickness/2) <(Nx-1)), ((Ny+1-dd) >0):Ny-1, *] = tickcolor
    endfor                      ; i
    for j = 1B, Ly do begin
      fac = 0.4 * (2 + ((j mod 5) eq 0) + ((j mod 10) eq 0))/2.
      dd = fix(marg/2 * fac)
      y = round(hmarg+thickness-1+j*pix)
      ;; Along vertical left
      vidcube[*, 0:dd, ((y-thickness/2) >0):((y+thickness/2) <(Ny-1)), *] = tickcolor
      ;; Along vertical right
      vidcube[*, ((Nx-1-dd) >0):(Nx-1), ((y-thickness/2) >0):((y+thickness/2) <(Ny-1)), *] = tickcolor
    endfor                      ; j
  endif

  ;; Write the video
  write_video, outfile, vidcube, video_fps = fps, bit_rate = bit_rate 
  
end

