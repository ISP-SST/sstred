; docformat = 'rst'

;+
; Return the average (sum/number) of the image frames defined by a
; list of file names. 
;
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, ISP
; 
; 
; :Returns:
; 
;    The average image. 
; 
; 
; :Params:
; 
;    files_list : in, type=strarr
;   
;      The list of image files, the contents of which to sum.
;   
; 
; :Keywords:
; 
;    time_ave  : out, optional, type=string
;   
;       The average timestamp of the summed frames.
;   
;    time_beg  : out, optional, type=string
;   
;       The beginning of the first exposure in the summed frames.
;   
;    time_end  : out, optional, type=string
;
;       The end of the last exposure in the summed frames.
;   
;    summed  : out, optional, type=dblarr
;   
;       The summed frames, without division with the number of summed frames.
;   
;    nsum : out, optional
;
;       Number of frames actually summed
;
;    old :  in, optional, type=boolean 
;   
;       Set this for data with the "old" header format.
;   
;    check  : in, optional, type=boolean
;   
;       Check for bad frames before summing.
;   
;    lun : in, optional, type=integer 
;   
;       Logical unit number for where to store checking results.
;   
;    gain : in, optional
; 
;       Gain frame with bad pixels zeroed. If backscatter correction
;       is supposed to be performed, the gain must have been
;       calculated from a backscatter corrected flat.
; 
;    dark : in, optional
; 
;       Dark frame.
; 
;    select_align : in, optional, type=boolean
; 
;       Align before summing if set. User gets to select feature to
;       align on. 
;
;    pinhole_align : in, optional, type=boolean
; 
;       Align before summing if set. Use brightest feature (assumed to
;       be a pinhole) to align on. If both select_align and
;       pinhole_align are set, pinhole_align is ignored.
;
;    nthreads  : in, optional, type=integer
;   
;       The number of threads to use for backscatter correction.
;   
;    backscatter_gain : in, optional
;   
;       The gain needed to do backscatter correction. Iff
;       backscatter_gain and backscatter_psf are given, do the
;       correction. 
;   
;    backscatter_psf : in, optional
;   
;       The psf needed to do backscatter correction. Iff
;       backscatter_gain and backscatter_psf are given, do the
;       correction.
;   
;   xyc : out, optional, type="lonarr(2)" 
;   
;       The (x,y) coordinates of the feature used for alignment.
;   
;    filter  : in, optional, type=int, default=3
;       Size of the medianfilter to use when checking data.
;   
; 
; :History:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2013-08-15 : MGL. Added documentation. Code from last year that
;                allows subpixel alignment before summing.
;
;   2013-08-16 : MGL. New code to optionally do backscatter
;                correction, as well as gain and dark correction, and
;                filling of bad pixels. Removed unused keyword
;                "notime". Use red_time2double rather than
;                time2double.
; 
;   2013-08-27 : MGL. Let the subprogram find out its own name.
; 
;   2013-09-09 : MGL. Move per-frame correction for dark, gain, and
;                back scatter to inside of conditional for pinhole
;                alignment. 
; 
;   2013-09-13 : MGL. Lower the limit for least number of frames
;                needed to do checking from 10 to 3.
; 
;   2013-12-10 : PS. Keyword lim
; 
;   2013-12-11 : PS. Also pass back the number of summed frames
; 
;   2014-01-27 : MGL. Use red_com rather than red_centroid.
; 
;   2014-04-11 : MGL. New keywords select_align and xyc. Limit number
;                of repeats in the alignment so bad data can
;                terminate.
; 
;   2016-05-20 : MGL. Rewrite the setting up and reading part to
;                support files with more than a single frame. Use
;                red_readdata and red_redhead.
; 
;   2016-05-22 : MGL. Rewrite the summing part, take care of the easy
;                cases first.
; 
;   2016-05-22 : MGL. Remove "old" keyword. If we have "old" ana
;                headers, they should be taken care of in
;                red_readdata. 
; 
;   2016-05-23 : MGL. Added time_beg and time_end keywords. Changed
;                time keyword to time_ave. Do summing in double
;                precision.  
; 
;   2016-05-24 : MGL. Removed filtering of FITS headers.  
; 
;   2016-05-25 : MGL. Use red_progressbar.  
; 
;   2016-05-30 : MGL. Report how many files are being summed.   
; 
;   2016-06-09 : MGL. Bugfix in printout of what frames were rejected.
;                Bugfix in summing without checking.
; 
;   2016-09-21 : THI. Make the size of the medianfilter a parameter.   
; 
;   2016-09-22 : MGL. Make time_beg etc. with 6 decimals.
;
; 
;-
function red_sumfiles, files_list $
                       , time_ave = time_ave $
                       , time_beg = time_beg $
                       , time_end = time_end $
                       , summed = summed $
                       , nsum = nsum $
                       , check = check $ 
                       , lim = lim $
                       , lun = lun $ 
                       , pinhole_align = pinhole_align $ 
                       , select_align = select_align $ 
                       , gain = gain $   
                       , dark = dark $   
                       , xyc = xyc $   
                       , backscatter_gain = backscatter_gain $
                       , backscatter_psf = backscatter_psf $
                       , nthreads = nthreads $
                       , filter = filter

  ;; Name of this subprogram
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  if n_elements(nthreads) eq 0 then nthreads = 2 else nthreads = nthreads

  DoDescatter = n_elements(backscatter_gain) gt 0 and n_elements(backscatter_psf) gt 0

  Nfiles = n_elements(files_list)

  ;; Find out dimensions of the data
  Nframes_per_file = lonarr(Nfiles)
  Nframes = 0
  for ifile = 0, Nfiles-1 do begin
    head = red_readhead(files_list[ifile], /silent)
    case sxpar(head, 'NAXIS') of
      2 : Nframes_per_file[ifile] = 1
      3 : Nframes_per_file[ifile] = sxpar(head, 'NAXIS3')
      else : begin
        print, inam + ' : No support for files other than single 2D frames or 3D data cubes' 
        print, files_list[ifile]
        print, head
        stop
      end
    endcase
    if ifile eq 0 then begin
      dim = [sxpar(head, 'NAXIS1'), sxpar(head, 'NAXIS2')]
    endif else begin
      if dim[0] ne sxpar(head, 'NAXIS1') or dim[1] ne sxpar(head, 'NAXIS2') then begin
        print, inam + ' : [X,Y] dimensions do not match earlier files.'
        print, files_list[ifile]
        print, 'This file:',  [sxpar(head, 'NAXIS1'), sxpar(head, 'NAXIS2')]
        print, 'Earlier file(s):', dim 
        stop
      endif
    endelse
  endfor                        ; ifile
  Nframes = total(Nframes_per_file)


  if keyword_set(check) then docheck = 1B else docheck = 0B
  if docheck and Nframes lt 3 then begin
    print, inam+" : Not enough statistics, don't do checking."
    docheck = 0B
  endif

  if keyword_set(pinhole_align) or keyword_set(select_align) then begin
    if n_elements(gain) eq 0 or n_elements(dark) eq 0 then begin
      print, inam+' : Need to specify gain and dark when doing /pinhole_align or /select_align'
      help, gain, dark
      stop
    endif
  endif

  if n_elements(gain) eq 0 then begin
    ;; Default gain is unity
    gain = dblarr(dim)+1D
  endif else begin
    ;; Make a mask for fillpixing. We want to exclude really large
    ;; areas along the border of the image and therefore irrelevant
    ;; for pinhole calibration.
    mask = red_cleanmask(gain)
  endelse 

  if n_elements(dark) eq 0 then begin
    ;; Default zero dark correction
    dark = dblarr(dim)
  endif 

  ;; If just a single frame, return it now. 
  if Nframes eq 1 then begin
    head = red_readhead(files_list[0], /silent)
    time_beg = red_time2double((strsplit(fxpar(head, 'DATE-BEG'), 'T', /extract))[1])
    time_end = red_time2double((strsplit(fxpar(head, 'DATE-END'), 'T', /extract))[1])
    time_ave = red_timestring((time_beg + time_end) / 2d, Nsecdecimals = 6)
    time_beg = red_timestring(time_beg, Nsecdecimals = 6)
    time_end = red_timestring(time_end, Nsecdecimals = 6)
    ;; Include gain, dark, fillpix, backscatter here? This case
    ;; should really never happen...
    return, red_readdata(files_list[0], /silent)
  endif

  ;; Needed for warning messages and progress bars.
  ntot = 100. / (Nframes - 1.0)
  bb = string(13B)

  times = dblarr(Nframes)
  times_beg = dblarr(Nframes)
  times_end = dblarr(Nframes)

  if docheck then begin
    ;; Set up for checking
    cub = intarr(dim[0], dim[1], Nframes)
  endif

  ;; Loop over files
  iframe = 0
  for ifile = 0, Nfiles-1 do begin
    
    if DoCheck then begin
      cub[0, 0, iframe] = red_readdata(files_list[ifile], header = head, /silent)
      red_progressbar, iframe, Nframes, inam+' : loading files in memory', clock = clock
    endif else begin
      head = red_readhead(files_list[ifile], /silent)
    endelse

    cadence = sxpar(head, 'CADENCE')
    time_beg = red_time2double((strsplit(fxpar(head, 'DATE-BEG'), 'T', /extract))[1])
    time_end = red_time2double((strsplit(fxpar(head, 'DATE-END'), 'T', /extract))[1])

    times[iframe] = time_beg + cadence*findgen(Nframes_per_file[ifile]) + fxpar(head, 'XPOSURE')/2

    times_beg[iframe] = time_beg + cadence*findgen(Nframes_per_file[ifile])
    times_end[iframe] = (times_beg[iframe:iframe+Nframes_per_file[ifile]-1] + fxpar(head, 'XPOSURE')) <time_end

    iframe += Nframes_per_file[ifile]
    
  endfor                        ; ifile
  print

  if DoCheck then begin

    mval = total(total(cub, 1), 1)/(dim[0]*dim[1])

    ;; Find bad frames
    if n_elements(lim) eq 0 then lim = 0.0175   ; Allow 2% deviation from the mean value
    if n_elements(filter) eq 0 then filter = 3  ; default is to use a medianfilter of width 3
    if filter mod 2 eq 0 then filter += 1       ; force odd filtersize
    tmean = median(mval)
    mmval = median( mval, filter )
    ;; Set edges to neighbouring values since the median filter does
    ;; not modify endpoints.
    mmval[0:filter/2-1] = mmval[filter/2]        
    mmval[Nfiles-filter/2:Nfiles-1] = mmval[Nfiles-filter/2-1] 

    goodones = abs(mval - mmval) LE lim * tmean ; Unity for frames that are OK.
    idx = where(goodones, Nsum, complement = idx1, Ncomplement = Nrejected)     

    if Nrejected gt 0 then begin
      print, inam+' : rejected frames :'
      for irejected = 0, Nrejected-1 do begin 
        ;; Find file name and frame number within file
        ifile = 0
        while idx1[irejected] gt total(Nframes_per_file[0:ifile]) do ifile += 1
        if ifile gt 0 then Nsub = round(total(Nframes_per_file[0:ifile-1])) else Nsub = 0
        msg = files_list[ifile]+', frame # '+strtrim(idx1[irejected]-Nsub, 2)
        print, msg
        if(keyword_set(lun)) then printf, lun, msg
      endfor                    ; irejected
    endif else print, inam+' : all files seem to be fine'

  endif else begin              ; docheck
    
    ;; If no checking, all frames are considered OK.
    goodones = replicate(1, Nframes) 
    idx = where(goodones, Nsum, complement = idx1)
    
  endelse                       ; docheck

  ;; Set time stamps to potentially be returned as keywords.
  time_beg = red_timestring(min(times_beg[idx]), Nsecdecimals = 6)
  time_ave = red_timestring(mean(times[idx]),    Nsecdecimals = 6)
  time_end = red_timestring(max(times_end[idx]), Nsecdecimals = 6)

  
  ;; Do the summing
  
  case 1 of

    ~DoCheck $
       and ~DoDescatter $
       and ~keyword_set(pinhole_align) $
       and ~keyword_set(select_align) : begin
      
      ;; Easy case first, no checking, no alignment and no
      ;; descattering. No checking means we have not read the data
      ;; yet. 

      ;; Tested with dark frames on 2016-05-23.

      summed = dblarr(dim)
      iframe = 0

      for ifile = 0, Nfiles-1 do begin
        
        red_progressbar, clock = clock, iframe, Nframes $
                         , inam+' : reading and summing '+strtrim(Nsum, 2)+' files'
        cub = red_readdata(files_list[ifile], /silent)

        summed += total(cub, 3, /double)
        iframe += Nframes_per_file[ifile]

      endfor                    ; ifile

      average = summed / Nsum

      ;; Dark and gain correction 
      average = (average - dark)*gain

      ;; Fill the bad pixels 
      if n_elements(mask) gt 0 then average = red_fillpix(average, mask = mask)

    end

    ~DoDescatter $
       and ~keyword_set(pinhole_align) $
       and ~keyword_set(select_align) : begin

      ;; Another easy case, checked data but no alignment and no
      ;; descattering. The data are checked means they are already
      ;; read in.

      ;; Tested with dark frames on 2016-05-23.

      print, bb, inam+' : summing '+strtrim(Nsum, 2)+' good files in memory.'

      summed = total(cub[*, *, idx], 3, /double)

      average = summed / Nsum

      ;; Dark and gain correction 
      average = (average - dark)*gain

      ;; Fill the bad pixels 
      if n_elements(mask) gt 0 then average = red_fillpix(average, mask = mask)

    end
    
    else : begin
      
      ;; And now for the more time consuming cases, where all frames
      ;; have to be dealt with individually. Either for alignment,
      ;; or for descattering, or both.
      
      
      firstframe = 1B           ; When doing alignment we use the first (good) frame as reference.
      dc_sum = [0.0, 0.0]
      summed = dblarr(dim)

      ifile = 0
      ii = 0                    ; Within-file-counter

      ;; Loop over all frames
      for iframe = 0L, Nframes-1 do begin
        
        if goodones[ii] then begin ; Sum only OK frames

          if docheck then begin

            ;; If checked, we already have the frames in memory.
            thisframe = double(cub[*, *, iframe])
            red_progressbar, clock = clock, iframe, Nframes $
                             , inam+' : summing '+strtrim(Nsum, 2)+' checked frames'

          endif else begin

            ;; If not checked, we (sometimes) have to read the frames in.
            if ii eq 0 or ii ge Nframes_per_file[ifile] then begin
              cub = red_readdata(files_list[ifile], /silent)
              ii = 0
            endif
            thisframe = double(cub[*, *, ii])
            ii += 1
            
            red_progressbar, clock = clock, iframe, Nframes $
                             , inam+' : summing '+strtrim(Nsum, 2)+' frames'

          endelse
          
          if keyword_set(pinhole_align) or keyword_set(select_align) then begin
            
            ;; If we are doing sub-pixel alignment, then we need to
            ;; correct each frame for dark and gain.
            thisframe = (thisframe - dark)*gain

            ;; We also need to do any descattering correction of each frame.
            if DoDescatter then begin
              thisframe = red_cdescatter(thisframe $
                                         , backscatter_gain, backscatter_psf $
                                         , /verbose, nthreads = nthreads)
            endif

            ;; And fill the bad pixels 
            thisframe = red_fillpix(thisframe, mask = mask)

            if firstframe then begin
              
              marg = 100
              if keyword_set(select_align) then begin
                ;; Select feature to align on with mouse.
                if max(dim) gt 1000 then fac = max(dim)/1000. else fac = 1
                window, 0, xs = 1000, ys = 1000
                tvscl, congrid(thisframe, dim[0]/fac, dim[1]/fac, cubic = -0.5)
                print, 'Use the mouse to click on feature to align on.'
                cursor, xc, yc, /device
                xyc = round([xc, yc]*fac) >marg <(dim-marg-1)
                subsz = 300
                subim = thisframe[xyc[0]-subsz/2:xyc[0]+subsz/2-1, xyc[1]-subsz/2:xyc[1]+subsz/2-1]
                tvscl, subim
              end else begin
                ;; Find brightest pinhole spot that is reasonably centered in
                ;; the FOV.
                subim = thisframe[marg:dim[0]-marg, marg:dim[1]-marg]
                mx = max(subim, maxloc)
                ncol = dim[1]-2*marg+1
                xyc = lonarr(2)
                xyc[0] = maxloc / ncol
                xyc[1] = maxloc MOD ncol
                xyc += marg     ; Position of brightest spot in original image
              endelse

              ;; Establish subfield size sz, shrunk to fit.
              sz = 99              
              subim = thisframe[xyc[0]-sz/2:xyc[0]+sz/2, xyc[1]-sz/2:xyc[1]+sz/2]
              tot_init = total(subim/max(subim) gt 0.2)
              repeat begin
                sz -= 2
                subim = thisframe[xyc[0]-sz/2:xyc[0]+sz/2, xyc[1]-sz/2:xyc[1]+sz/2]
                tot = total(subim/max(subim) gt 0.2)
              endrep until tot lt tot_init

              sz += 4           ; A bit of margin
              
            endif               ; firstframe

            ;; Iteratively calculate centroid of thisframe, and then
            ;; shift the frame so it aligns with the first frame.
            
            ;; Iteratively shift the image to get the spot centered in
            ;; the subfield, because that's where the centroiding is
            ;; most accurate.
            dc1 = [0.0,0.0]
            Nrep = 0
            repeat begin
              ;;print, 'Shift:', dc1
              im_shifted = red_shift_im(thisframe, dc1[0], dc1[1])                        ; Shift the image
              subim0 = im_shifted[xyc[0]-sz/2:xyc[0]+sz/2, xyc[1]-sz/2:xyc[1]+sz/2]       ; New subimage

              subim0 = subim0 / stdev(subim0)
              cnt = red_com(subim0) ; Centroid after shift
              dcold = dc1           ; Old shift
              dc1 = dc1 + (sz/2.0 - cnt)
              ;;print, 'Shift change vector length:',
              ;;sqrt(total((dc1-dcold)^2))
              Nrep += 1
            endrep until sqrt(total((dc1-dcold)^2)) lt 0.01 or Nrep eq 100
            ;; Iterate until shift changes less than 0.01 pixel

            if firstframe then begin
              ;; Keep as reference
              dc0 = dc1
            endif else begin
              dc = dc1-dc0       ; This is the shift!
              dc_sum += dc       ; ...add it to the total
              thisframe = red_shift_im(thisframe, dc[0], dc[1]) 
              firstframe = 0B
            endelse             ; firstframe

          endif                 ; keyword_set(pinhole_align)

          summed += thisframe

        endif                   ; goodones[ii] 

      endfor                    ; iframe

      average = summed / Nsum

      if total(abs(dc_sum)) gt 0 then begin ; shift average image back, to ensure an average shift of 0.
        average = red_shift_im( average, -dc_sum[0]/Nsum, -dc_sum[1]/Nsum )
      endif

      if ~keyword_set(pinhole_align) and ~keyword_set(select_align) then begin

        ;; Some actions already done for each frame in the case of
        ;; pinhole alignment.

        ;; Dark and gain correction 
        average = (average - dark)*gain

        if DoDescatter then begin
          ;; Backscatter correction (done already for pinholes)
          average = red_cdescatter(average, backscatter_gain, backscatter_psf $
                                   , /verbose, nthreads = nthreads)
        endif

        ;; Fill the bad pixels 
        if n_elements(mask) gt 0 then  average = red_fillpix(average, mask = mask)

      endif

    end
    
  endcase

  return, average
  
end                             ; red_sumfiles

dnames = file_search('/mnt/sand15n/Incoming/2016.05.09/Dark/12:43:00/Chromis-W/cam*fits', count = Ndarks)

print, Ndarks

dsum_check = red_sumfiles(dnames[0:2] $
                          , time_ave = time_ave_check $
                          , time_beg = time_beg_check $
                          , time_end = time_end_check $
                          , summed = summed_check $
                          , nsum = nsum_check $  
                          , /check $
                         )

dsum = red_sumfiles(dnames[0:2] $
                    , time_ave = time_ave $
                    , time_beg = time_beg $
                    , time_end = time_end $
                    , summed = summed $
                    , nsum = nsum $
                   )

stats, dsum
stats, dsum_check
stats, dsum_check - dsum


end
