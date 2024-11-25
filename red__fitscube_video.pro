; docformat = 'rst'

;+
; Make a video file from data in a fitscube.
;
; Default is to make a temporal movie but it is also possible to make
; a spectral movie.
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
;    crop : in, optional, type=boolean
;
;      Crop by use of a GUI.
;   
;    format : in, optional, string, default='mp4'
;   
;      The container format of the movies. If 'mp4' or 'avi',
;      write_video will make such a file. If 'mov' (Mac-friendly), a
;      spawned ffmpeg command will convert write_video's mp4 output to
;      the desired format.
;
;    gamma : in, optional, type=float, default=1
;
;      Non-linear intensity scale gamma.
;
;    golden : in, optional, type=boolean
;
;      Set rgbgamma to produce a golden color scheme.
;
;    iwave : in, optional, type=integer, default=0
;
;      The index of the wavelength to be shown in the movie. The
;      string "last" is also allowed.
;
;    istokes : in, optional, type=integer, default=0
;
;      The index of the Stokes component to be shown in the movie,
;      [I,Q,U,V] correspond to indices [0,1,2,3].
;
;    outdir : in, optional, type=string, default=file_dirname(infile)
;
;      The directory in which to write the movie file.
;
;    outfile : in, optional, type=string, default='same as infile but with mp4 extension'
;
;      The name of the output file. Note that periods in this file
;      name are not allowed (except as a delimiter for the file
;      extension). The extension of this file name takes precedence
;      over the format keyword.
;                    
;    rgbgamma : in, optional, type=fltarr(3), default="[1,1,1]"
;
;      Gamma values for the RGB channels. The total gamma correction
;      is of the form data[channel]^(gamma*rgbgamma[channel]). See
;      also keyword golden.
;                    
;    rgbfac : in, optional, type=fltarr(3), default="[1,1,1]"
;
;      Multiplicative factors for the RGB channels.
;
;    shrink_fac : in, optional, type=integer
;
;      Shrink FOV dimensions by this factor.
;
;    tickcolor : in, optional, type=byte, default=0
;   
;      The color of the tickmarks, 0=black, 255=white.   
;
;    tickmarks : in, optional, type=boolean
;   
;      Set this to draw arcsec tickmarks.
;   
;    video_codec : in, optional, type=string, default='mjpeg'
;   
;      The codec to use when making the video. Used when calling
;      write_video.
;   
;    video_fps : in, optional, type=integer, default=8
;   
;      Frames per second in the movie. Used when calling write_video.
; 
; :History:
; 
;    2018-01-10 : MGL. First version.
; 
;    2018-01-12 : MGL. Use subroutine red_fitscube_getframe rather
;                 than method fitscube_getframe. New keyword: crop.
; 
;    2018-11-29 : MGL. Make also spectral movies.
; 
;    2020-12-10 : MGL. Make .mov files without conversion from .mp4.
;                 New keyword shrink_fac.
; 
;-
pro red::fitscube_video, infile $
                         , bit_rate = bit_rate $
                         , clip = clip $
                         , crop = crop $
                         , format = format $
                         , gamma = gamma $
                         , golden = golden $
                         , istokes = istokes $
                         , iscan = iscan $
                         , ituning = ituning $
                         , outdir = outdir $
                         , outfile = outfile $
                         , rgbfac = rgbfac $
                         , rgbgamma = rgbgamma $
                         , scannumber = scannumber $
                         , shrink_fac = shrink_fac $
                         , spectral = spectral $
                         , tickcolor = tickcolor $
                         , tickmarks = tickmarks $
                         , video_codec = video_codec $
                         , video_fps = video_fps
  
  inam = red_subprogram(/low, calling = inam1)

  ;; Check the input file
  if ~file_test(infile) then begin
    print, inam + ' : Input file does not exist: '+infile
    retall
  endif

  if n_elements(format) eq 0 then format = 'mov'
  case format of
    'avi' : extension = format
    'mp4' : extension = format
    else  : extension = 'mp4'
  end
  if n_elements(video_fps) eq 0 then video_fps = 5
  if n_elements(bit_rate) eq 0 then bit_rate = 40000
  if n_elements(video_codec) eq 0 then video_codec = 'mjpeg'

  ;; Gamma and color settings
  if n_elements(outdir) eq 0 then outdir = file_dirname(infile)
  if n_elements(rgbfac) eq 0 then rgbfac = [1., 1., 1.]
  if keyword_set(golden) then rgbgamma = [0.6, 1.0, 6.0]  
  if n_elements(gamma) eq 0 then gamma = 1.0
  if n_elements(rgbgamma) eq 0 then rgbgamma = [1., 1., 1.]
  totalgamma = gamma*rgbgamma

  if n_elements(clip) eq 0 then clip = 2
  if n_elements(tickcolor) eq 0 then tickcolor = 0B else tickcolor = byte(tickcolor)
  if n_elements(istokes) eq 0 then begin
    istokes = 0
  endif else begin
    if size(istokes, /tname) eq 'STRING' then begin
      case strupcase(istokes) of
        'I' : istokes = 0
        'Q' : istokes = 1
        'U' : istokes = 2
        'V' : istokes = 3 
      endcase
    endif 
  endelse
  
  hdr = headfits(infile)
  dims = fxpar(hdr, 'NAXIS*')
  Ndims = n_elements(dims)

  if Ndims lt 3 then stop

  Nx = dims[0]
  Ny = dims[1]

  pref = strtrim(red_fitsgetkeyword(hdr, 'FILTER1'), 2)
  image_scale = self -> imagescale(pref, /use_config)

  if keyword_set(spectral) then begin
    
    Nframes = dims[2]
    vidcube = dblarr([Nx, Ny, Nframes])

    if n_elements(iscan) eq 0 then iscan = 0

    ;; Read the data
    for ituning = 0, Nframes-1 do begin
      red_fitscube_getframe, infile, frame, ituning = ituning, istokes = istokes, iscan = iscan
      vidcube[*, *, ituning] = frame
    endfor                      ; ituning

    ;; Default is to scale frames individually
    if n_elements(scale_frames) eq 0 then scale_frames = 1
    
  endif else begin

    Nframes = dims[4]
    vidcube = dblarr([Nx, Ny, Nframes])
    
    if n_elements(ituning) eq 0 then begin
      ituning = 0
    endif else begin
      if size(ituning, /tname) eq 'STRING' then begin
        case strlowcase(ituning) of
          'last' : ituning = dims[2]-1
          else   : ituning = 0
        endcase
      endif  
    endelse 

    ;; Read the data
    for iscan = 0, Nframes-1 do begin
      red_fitscube_getframe, infile, frame, ituning = ituning, istokes = istokes, iscan = iscan
      vidcube[*, *, iscan] = frame
    endfor                      ; iscan

    ;; Default is to scale the cube as a whole.
    if n_elements(scale_frames) eq 0 then scale_frames = 0
    
  endelse

  ;; Shrink and crop before scaling in case the structure within the
  ;; FOV would change the scaling.

  if n_elements(shrink_fac) ne 0 then begin
    Nxx = Nx - (Nx mod shrink_fac)
    Nyy = Ny - (Ny mod shrink_fac)
    vidcube = red_crop(vidcube, size = [Nxx, Nyy], /center)
    vidcube_shrunk = dblarr([Nxx/shrink_fac, Nyy/shrink_fac, Nframes])
    for iframe = 0, Nframes-1 do begin
      vidcube_shrunk[*, *, iframe] = rebin(vidcube[*, *, iframe] $
                                           , Nxx/shrink_fac $
                                           , Nyy/shrink_fac $
                                           , /samp)
    endfor                      ; iframe
    vidcube = vidcube_shrunk
  endif
  
  if keyword_set(crop) then vidcube = red_crop(vidcube)
  newdims = size(vidcube, /dim)
  Nx = newdims[0]
  Ny = newdims[1]  

  if keyword_set(scale_frames) then begin

    n = Nx*Ny

    for ituning = 0, Nframes-1 do begin
      
      if clip ne 0 then begin
        ;; Code to do histogram clip, stolen from cgclipscl.pro
        maxr = max(vidcube[*, *, ituning], MIN=minr, /NAN)
        range = maxr - minr
        binsize = range / 1000.
        h = histogram(vidcube[*, *, ituning], BINSIZE=binsize, OMIN=omin, OMAX=omax, /NAN)
        cumTotal = total(h, /CUMULATIVE)
        minIndex = value_locate(cumTotal, n * (clip/100.))
        if minIndex eq -1 then minIndex = 0
        while cumTotal[minIndex] eq cumTotal[minIndex + 1] do minIndex = minIndex + 1
        minThresh = minIndex * binsize + omin
        ;; Not all files can be clipped appropriately.
        maxIndex  = value_locate(cumTotal, n * ((100-clip)/100.))
        if (maxIndex eq -1) || (maxIndex eq N_Elements(cumTotal)) || (maxIndex eq minIndex) then stop
        ;; If you are still here, try to clip the histogram.
        while cumTotal[maxIndex] eq cumTotal[maxIndex - 1] do maxIndex = maxIndex - 1
        maxThresh = maxIndex * binsize + omin

        ;; Do the clipping
        vidcube[*, *, ituning] = vidcube[*, *, ituning] > minThresh < maxThresh
      endif 

      ;; Scale to [0,1] range
      vidcube[*, *, ituning] -= min(vidcube[*, *, ituning])
      vidcube[*, *, ituning] /= max(vidcube[*, *, ituning])

    endfor                      ; ituning
    
  end else begin

    if clip ne 0 then begin
      ;; Code to do histogram clip, stolen from cgclipscl.pro
      maxr = max(vidcube, MIN=minr, /NAN)
      range = maxr - minr
      binsize = range / 1000.
      h = histogram(vidcube, BINSIZE=binsize, OMIN=omin, OMAX=omax, /NAN)
      n = n_elements(vidcube)
      cumTotal = Total(h, /CUMULATIVE)
      minIndex = Value_Locate(cumTotal, n * (clip/100.))
      if minIndex eq -1 then minIndex = 0
      while cumTotal[minIndex] eq cumTotal[minIndex + 1] do minIndex = minIndex + 1
      minThresh = minIndex * binsize + omin
      ;; Not all files can be clipped appropriately.
      maxIndex  = value_locate(cumTotal, n * ((100-clip)/100.))
      if (maxIndex eq -1) || (maxIndex eq N_Elements(cumTotal)) || (maxIndex eq minIndex) then stop
      ;; If you are still here, try to clip the histogram.
      while cumTotal[maxIndex] eq cumTotal[maxIndex - 1] do maxIndex = maxIndex - 1
      maxThresh = maxIndex * binsize + omin

      ;; Do the clipping
      vidcube = vidcube > minThresh < maxThresh
    endif
  
    vidcube -= min(vidcube)
    vidcube /= max(vidcube)

  endelse

  ;; Do the RGB stuff
  rgbcube = fltarr(3, Nx, Ny, Nframes)
  rgbcube[0, *, *, *] = rgbfac[0]*vidcube^totalgamma[0]
  rgbcube[1, *, *, *] = rgbfac[1]*vidcube^totalgamma[1]
  rgbcube[2, *, *, *] = rgbfac[2]*vidcube^totalgamma[2]

  ;; Should be byte scaled
  rgbcube *= 255d
  rgbcube = byte(round(rgbcube >0 <255))

  ;; Tickmarks?
  if keyword_set(tickmarks) then begin
    ;; Code inspired by Luc's tickbox_f.ana
    marg = min([Nx, Ny])/20 >30 ;30
    hmarg = 0
    sc = float(image_scale)
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
      rgbcube[*, ((x-thickness/2) >0):((x+thickness/2) <(Nx-1)), 0:dd, *] = tickcolor
      ;; Along horizontal top
      rgbcube[*, ((x-thickness/2) >0):((x+thickness/2) <(Nx-1)), ((Ny+1-dd) >0):Ny-1, *] = tickcolor
    endfor                      ; i
    for j = 1B, Ly do begin
      fac = 0.4 * (2 + ((j mod 5) eq 0) + ((j mod 10) eq 0))/2.
      dd = fix(marg/2 * fac)
      y = round(hmarg+thickness-1+j*pix)
      ;; Along vertical left
      rgbcube[*, 0:dd, ((y-thickness/2) >0):((y+thickness/2) <(Ny-1)), *] = tickcolor
      ;; Along vertical right
      rgbcube[*, ((Nx-1-dd) >0):(Nx-1), ((y-thickness/2) >0):((y+thickness/2) <(Ny-1)), *] = tickcolor
    endfor                      ; j
  endif


  
  ;; Construct the output file name
  if n_elements(outfile) eq 0 then begin
    
    tag = '_'+(['I', 'Q', 'U', 'V'])[istokes]

    if keyword_set(spectral) then begin
      dum = red_fitsgetkeyword(infile, 'SCANNUM',variable_values=scannum)
      scannumbers = reform(scannum.values)
      if n_elements(scannumber) ne 0 then begin
        iscan = (where(scannumbers eq scannumber))[0]
      endif
      tag += '_scan'+strtrim(scannumbers[iscan], 2)
    endif else begin
      red_fitscube_getwcs, infile, coordinates = coordinates
      tag += '_' + string(coordinates[ituning,0].wave[0,0], format = '(f7.3)')+'nm'
    endelse
    
    outfile = file_basename(infile, '.fits') + tag 

    ;; Note that periods (.) are not allowed in the file name of the movie.
    outfile = red_strreplace(outfile, '.', '_', n = 100)

    outfile += '.' + extension
    
  endif

  
  ;; Write the video
  file_delete, outdir + '/' + outfile, /allow
  write_video, outdir + '/' + outfile $
               , rgbcube $
               , video_fps = video_fps $
               , video_codec = video_codec $
               , bit_rate = bit_rate 

  if format eq 'mov' then begin
    ;; Convert to Mac-friendly (and smaller) .mov file using recipe from Tiago     
    mname = outdir + '/' + red_strreplace(outfile, '.'+extension,'.'+format)
    file_delete,mname, /allow
    spawn, 'ffmpeg -n -i "' + outdir + '/' + outfile $
           + '" -c:v libx264 -preset slow -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -crf 26 -tune grain "' $
           + mname + '"'
    file_delete, outdir + '/' + outfile
  endif
  
end

