; docformat = 'rst'

;+
; Read data frames from multiple files. Assume all frames have the
; same dimensions.
;
; Keywords will be used when calling red_readdata for the individual
; files. 
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
; :Returns:
;
;    A datacube with the data from the files.
; 
; :Params:
; 
;   files : The file names.
; 
; :Keywords:
;
;    date_beg : out, optional, type=strarr
;
;       DATE-BEG keywords for all frames.
; 
;   framenums : in, out, optional, array
; 
;      If given as input, frames with these frame numbers are
;      returned. As output, the frame numbers of the returned frames
;      or -1 for frames that were not read.
; 
;   status : out, optional
; 
;      Status numbers for each file in fnames, 0 for success, -1 for
;      failure. .
; 
; :History:
; 
;   2022-08-22 : MGL. First version.
; 
;   2022-12-22 : MGL. Allow the files to have multiple frames.
; 
;   2024-09-06 : MGL. New keyword framenums.
; 
;-
function red_readdata_multiframe, files $
                                  , date_beg = date_beg $
                                  , framenums = framenums $
                                  , status = status $
                                  , _ref_extra = extra

  Nfiles = n_elements(files)

  status = bytarr(Nfiles)

  ;; Get all frame numbers from the files
  framecnt = red_fitsgetkeyword_multifile(files, 'NAXIS3') >1
  framenum = red_fitsgetkeyword_multifile(files, 'FRAMENUM')
  frameinc = red_fitsgetkeyword_multifile(files, 'FRAMEINC') >1
  for ifile = 0, Nfiles-1 do red_append, all_framenums, framenum[ifile] + indgen(framecnt[ifile])*frameinc[ifile]

  ;; Find out spatial dimensions of datacube
  hdr = red_readhead(files[0])
  Nx = red_fitsgetkeyword(hdr, 'NAXIS1', count=cnt)
  Ny = red_fitsgetkeyword(hdr, 'NAXIS2', count=cnt)

  if n_elements(framenums) eq 0 then begin

    ;; Simplified reading if all frames are returned.
    
    framenums = all_framenums   ; To be returned

    Nframes = n_elements(framenums)
    datacube = fltarr(Nx, Ny, Nframes)
    date_beg = strarr(Nframes)
    
    for ifile =  0, Nfiles-1 do begin
      
      red_progressbar, ifile, Nfiles, 'Reading data frames from multiple files'
      
      if ifile eq 0 then iframe = 0 else iframe = total(framecnt[0:ifile-1])

      datacube[0, 0, iframe] = red_readdata(files[ifile] $
                                            , date_beg = these_date_beg $
                                            , _strict_extra = extra $  
                                            , status = thisstatus)        

      ;; Status = -1 if failed, 0 if success.
      status[ifile] = thisstatus

      ;; Also set relevant framenums to -1 if reading failed!
      if thisstatus eq -1 then framenums[iframe] = replicate(-1, framecnt[ifile])

      if thisstatus eq 0 then date_beg[iframe] = these_date_beg
      
    endfor                      ; ifile
    
    return, datacube

  endif else begin

    ;; Called with framenums, so we need to select those frame numbers
    ;; only. We do not want to read all frames in all files and keep
    ;; them in memory so we'll read one file at a time and
    ;; select based on the frame numbers.

    Nframes = n_elements(framenums)
    datacube = fltarr(Nx, Ny, Nframes)
    date_beg = strarr(Nframes)

    for ifile =  0, Nfiles-1 do begin
      
      red_progressbar, ifile, Nfiles, 'Reading data frames from multiple files'

      these_framenums = framenum[ifile] + indgen(framecnt[ifile])*frameinc[ifile]

      match2, framenums, these_framenums, suba, subb

      Nmatch = round(total(suba ge 0))
      if Nmatch ne round(total(subb ge 0)) then stop ; Should not happen?

      if Nmatch eq 0 then begin
        status[ifile] = 0       ; If not success, at least not a failure!
        continue
      endif
      
      ;; There is at least one match, so we need to read this file
      data = red_readdata(files[ifile] $
                          , date_beg = these_date_beg $
                          , _strict_extra = extra $  
                          , status = thisstatus)     

      ;; Status = -1 if failed, 0 if success.
      status[ifile] = thisstatus

      if thisstatus eq -1 then begin
        ;; Also set matching framenums to -1 if reading failed!
        framenums[ifile] = -1
        continue                ; Nothing more to do with this file
      endif
        
      ;; Put frames corresponding to matching framnumbers in datacube
      for iframe = 0, framecnt[ifile]-1 do begin
        indx = where(subb ne -1) ; These frames from this file match
        datacube[*, *, subb[indx]] = data[*, *, indx]
;        if subb[iframe] eq -1 then continue ; This is not a matching frame
;        ;;  print, framenums[subb[iframe]], these_framenums[iframe]
;        datacube[*, *, subb[iframe]] = data[*, *, iframe]
        date_beg[subb[indx]] = these_date_beg[indx]
      endfor                    ; iframe
      
    endfor                      ; ifile
    
    return, datacube

  endelse
  
end

files = file_search('/scratch/mats/2023-10-17/CRISP/data/08:39:58/Crisp-W/*_00000_*6173*', count = Nfiles)

help, files
ims1 = red_readdata_multiframe(files, framenums = framenums1, status = status1, date_beg = date_beg1)

framenums2 = [1024, 1036, 1048, 1200, 1300, 1400]
ims2 = red_readdata_multiframe(files, framenums = framenums2, status = status2, date_beg = date_beg2)

stats,ims1[*,*,6] - ims2[*,*,1]
stats,ims1[*,*,3] - ims2[*,*,0]
stats,ims1[*,*,282] - ims2[*,*,4]

time_beg1 = red_time2double(strmid(date_beg1, 11))
time_beg2 = red_time2double(strmid(date_beg2, 11))

red_timeplot, time_beg1, framenums1, psym = 16, /yno
cgplot, /over, time_beg2, framenums2, psym = 16, color = 'cyan'

end
