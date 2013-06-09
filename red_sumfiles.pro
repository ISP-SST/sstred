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
; 
; 
; 
; :returns:
; 
; 
; :Params:
; 
;    files_list : 
;   
;   
;   
; 
; :Keywords:
; 
;    time  : 
;   
;   
;   
;    summed  : 
;   
;   
;   
;    old  : 
;   
;   
;   
;    check  : 
;   
;   
;   
;    lun  : 
;   
;   
;   
;    notime : 
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
;-
function red_sumfiles, files_list, time = time, summed = summed, old = old, check = check, lun = lun, notime=notime
                                ;
  new = 1B
  If(keyword_set(old)) then new = 0B
                                ;
  nt = n_elements(files_list)
  fzread, tmp, files_list[0], h
  tmp = long(tmp)
                                ;
  dum = fzhead(files_list[0])
  dum = strsplit(dum, ' =', /extract)
  if(new) then begin
     t1 = dum[18]
     t2 = dum[21]
     time = (red_time2double(t1) + red_time2double(t2)) * 0.5d0
  endif else time = red_time2double(dum[2])
                                ;
  if(nt gt 1) then begin
     ntot = 100. / (nt - 1.0)
     bb = string(13B)
                                ;
     if(~keyword_set(check)) then begin
        for ii = 1L, nt -1 do begin
           tmp+= f0(files_list[ii])
           dum = fzhead(files_list[ii])
           dum = strsplit(dum, ' =', /extract)
           if(new) then begin
              t1 = dum[18]
              t2 = dum[21]
              time+= (red_time2double(t1) + red_time2double(t2)) * 0.5d0
           endif else time += red_time2double(dum[2])
                                ;
           print, bb, 'red_sumfiles : adding files -> ', ntot * ii, '%', $
                  FORMAT = '(A,A,F5.1,A,$)'
        endfor
        print, ' '
        summed = tmp
        tmp = float(double(tmp) / double(nt))
        time = red_time2double(time / double(nt), /dir)

     endif else begin
        dim = size(tmp, /dim)
        cub = intarr(dim[0], dim[1], nt)
        times = dblarr(nt)
        mval = fltarr(nt)
                                ;
        for ii = 0L, nt -1 do begin
           cub[*,*,ii] = f0(files_list[ii])
           mval[ii] = mean(cub[*,*,ii])
                                ;
           dum = fzhead(files_list[ii])
           dum = strsplit(dum, ' =', /extract)
           
           if(new) then begin
              t1 = dum[18]
              t2 = dum[21]
              times[ii] = (red_time2double(t1) + red_time2double(t2)) * 0.5d0
           endif else times[ii] = red_time2double(dum[2])
           print, bb, 'red_sumfiles : loading files in memory -> ', ntot * ii, '%', $
                  FORMAT = '(A,A,F5.1,A,$)'
        endfor
        print, ' '
        
        if(n_elements(lim) eq 0) then lim = 0.0175 ; Allow 2% deviation from the mean value
        tmean = median(mval)
        mmval = median(mval, 3)
        idx = where(abs(mval - mmval) LE lim * tmean, count, complement = idx1)
                                ;tt = times - times[0]
                                ; plot, tt, abs(mval-mmval), psym=3, ystyle = 3, xstyle = 3
                                ; oplot, tt,replicate(tmean + tmean*lim, nt), /line
                                ; loadct,3,/silent
                                ; if(count ne nt) then oplot, tt,(abs(mval-mmval))[idx1], psym=7, color = 165
                                ; loadct,0,/silent
        
        if(count ne nt) then begin
           print, 'red_sumfiles : rejected frames :'
           print, transpose((files_list[idx1]))
           IF(keyword_set(lun)) THEN BEGIN
              printf, lun, files_list[idx1]
           endif
        endif else print, 'red_sumfiles : all files seem to be fine'
                                ;stop
                                ;summed = total((temporary(cub))[*,*,idx], 3, /preserve) 
        summed = lonarr(dim[0], dim[1])
        for ii = 0, n_elements(idx)-1 do begin
           summed += cub[*,*,idx[ii]]
           print, string(13b), 'red_sumfiles : summing checked frames -> ', 100./(count-1.0) * ii, '%',FORMAT='(A,A,F5.1,A,$)'
        endfor
        print, ' '
        tmp = summed / float(count)
        time = red_time2double(mean(times[idx]), /dir)
     endelse
  endif

  return, tmp
END


;; Make (inverse) gain table from flat field (or sum thereof). MLö 2008
;; function red_flat2gain, flat $
;;                         , gain_nozero = gain_nozero $
;;                         , smoothsize = smoothsize $
;;                         , badthreshold = badthreshold $
;;                         , maxgain = maxgain $
;;                         , mingain = mingain $
;;                         , preserve_sarnoff_taps = preserve_sarnoff_taps 

;;   ;; If preserve_sarnoff_taps is set, don't zero borders
;;   ;; between Sarnoff taps.

;;   ;; Unsharp masking smoothing kernel width
;;   if n_elements(smoothsize) eq 0 then smoothsize = 7.    

;;   ;; Unsharp masking threshold for bad pixels 
;;   if n_elements(badthreshold) eq 0 then badthreshold = .8 

;;   ;; Another pair of thresholds, this one on the gain itself
;;   if n_elements(maxgain) eq 0 then maxgain = 5.
;;   if n_elements(mingain) eq 0 then mingain = 0.

;;   ;; First make the basic gain table, which is a normalized inverse
;;   ;; flat. 

;;   indx = where(flat)
;;   mq = median(flat(indx))
;;   gain = 0.0*flat
;;   gain(indx) = mq/flat(indx)    

;;   ;; Then find bad pixels.

;;   ;; Unsharp masking
;;   ;psf = get_psf(smoothsize*7, smoothsize*7,smoothsize,smoothsize)
;;   ;psf /= total(psf, /double)
;;   ;smgain = red_convolve(gain, psf)
;;   smgain = smooth(gain, smoothsize, /edge, /nan)
;;   dgain = abs(gain - smgain)

;;   ;; Find bad pixels
;;   if keyword_set(preserve_sarnoff_taps) then begin
;;      badindex = where((dgain gt badthreshold) $ ; array with badpixels
;;                       or (gain le mingain) $ 
;;                       or (gain gt maxgain) $
;;                       , complement = goodindx) ; and the good ones
;;   endif else begin
;;      ;; Zero borders between Sarnoff taps
;;      xx = indgen(size(gain,/dim), /long) mod (size(gain,/dim))[0]  
;;      yy = indgen(size(gain,/dim), /long) / (size(gain,/dim))[0]  

;;      badindex = where((dgain gt badthreshold) $ ; array with badpixels
;;                       or (gain le mingain) $ 
;;                       or (gain gt maxgain) $
;;                       or ((xx mod 128) eq 0) or ((yy mod 512) eq 0) $ ; taps
;;                       , complement = goodindx) ; and the good ones
;;   endelse
;;   print, 'There are ', (size(badindex, /dim))[0], ' bad pixels.'

;;   ;; Normalize with respect to the good pixels
;;   gain = gain/mean(gain[goodindx])

;;   ;; If requested:
;;   IF arg_present(gain_nozero) then gain_nozero = gain

;;   ;; Set the bad pixels to zero.
;;   if (size(badindex, /dim))[0] gt 0 then gain[badindex] = 0.0 

;;   return, gain

;; end                             ; red_flat2gain
;
function red_get_psf,nx,ny,cx,cy
  psf=dblarr(nx,ny)
  nx2=double(nx/2)
  ny2=double(ny/2)
  cyy=cy/(2.d0*sqrt(2.d0*alog(2.d0)))
  cxx=cx/(2.d0*sqrt(2.d0*alog(2.d0)))
  for yy=0.,ny-1. do for xx=0.,nx-1 do begin
     u=((xx-nx2)/cxx)^2.+((yy-ny2)/cyy)^2.
     psf[xx,yy]=exp(-u*0.5d0)
  endfor
  return,psf
end
