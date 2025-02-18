; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CHROMIS pipeline
; 
; 
; :Author:
; 
; 
; 
; 
; :Keywords:
;
;    files : in, optional, type=strarr
;
;       Flat files to make gains out of.
;
;    nthreads  : 
;   
;   
;   
;    cam  : 
;   
;   
;   
;    pref  : 
;   
;   
;   
;    min  : 
;   
;   
;   
;    max  : 
;   
;   
;   
;    bad : 
;   
;   
;   
;    smoothsize  : 
;   
;   
;   
; 
; 
; :History:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2013-08-27 : MGL. Added support for logging. Let the subprogram
;                find out its own name.
; 
;   2016-02-15 : MGL. Use loadbackscatter. Remove keyword descatter,
;                new keyword no_descatter.
; 
;   2017-04-06 : MGL. Make version for CHROMIS without backscatter or
;                zeroing of detector tap borders.
; 
;   2017-12-20 : MGL. Store gains with metadata headers. 
; 
;   2018-04-19 : MGL. Rewrote using filenames() method, separate
;                ordinary flats and cavity-free flats.
; 
;   2018-08-21 : MGL. New keyword files.
; 
;   2022-08-29 : MGL. CHROMIS --> RED.
;
;   2024-11-02 : JdlCR. Modifications for new
;                demodulation/flat-fielding scheme.
;
;-
pro red::makegains, bad=bad $
                    , cam = cam $
                    , files = files $
                    , max = max $
                    , min = min $
                    , nthreads = nthreads $
                    , pref = pref $
                    , smoothsize = smoothsize 

  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)    

  red_make_prpara, prpara, bad
  red_make_prpara, prpara, cam 
  red_make_prpara, prpara, max 
  red_make_prpara, prpara, min 
  red_make_prpara, prpara, pref 
  red_make_prpara, prpara, smoothsize

  if n_elements(files) eq 0 then begin
    tosearch = self.out_dir+'/flats/*flat.fits'
    files = file_search(tosearch, count = Nfiles)
    if Nfiles eq 0 then begin
      print, inam+' : No flats found : ' + tosearch
    endif
  endif

  ;; We want to run separately for regular and cavity-free flats 
  cindx = where(strmatch(files, '*cavityfree*'), complement = findx, Ncav, ncomplement = Nnocav)
  lccindx = where(strmatch(files[cindx], '*lc?_cavityfree*'), complement = oindx, Ncav, ncomplement = Noldcav)

  ;; split cavityfree types in two separate lists
  lccindx = cindx[lccindx]  
  cindx = cindx[oindx]
  

  for iscav = 0, 2 do begin

    if(iscav eq 0) then begin
      if(Nnocav gt 0) then begin
        flatname = files[findx]
      endif else begin
        continue
      endelse
    endif else if(iscav eq 1) then begin
      if(Noldcav gt 0) then begin
        flatname = files[cindx]
      endif else begin
        continue
      endelse
    endif else if(iscav eq 2) then begin
      if(Ncav gt 0) then begin
        flatname = files[lccindx]
      endif else begin
        continue
      endelse
    endif
   
    
    Nfiles = n_elements(flatname)
    if Nfiles eq 0 then continue
    
    self -> extractstates, flatname, states

    case iscav of 
      2: begin
        gainname = self -> filenames('cavityfree_lc_gain', states)
      end
      1: begin
        gainname = self -> filenames('cavityfree_gain', states)
      end
      0: begin
        gainname = self -> filenames('gain', states)
      end
    endcase
    
    for ifile = 0L, Nfiles -1 do begin
      
      tmp = strsplit(file_basename(flatname[ifile]), '._', /extract)
      if(keyword_set(pref)) then begin
        if(tmp[5] ne pref) then begin
          print, inam+' : skipping prefilter -> '+tmp[5]
          continue
        endif
      endif

      flat = red_readdata(flatname[ifile], head = hdr)
      
      ;; Only one camera?
      if n_elements(cam) ne 0 then if tmp[0] NE cam then continue
      
      gain = self -> flat2gain(flat, ma=max, mi=min, bad=bad $
                               , /preserve, smoothsize=smoothsize)
      
;    namout = file_basename(flatname[ifile], '.flat.fits')+'.gain'
;    outdir = self.out_dir+'/gaintables/'
      
      ;; Edit the header
      red_fitsaddkeyword, hdr, 'FILENAME', file_basename(gainname[ifile])
      self -> headerinfo_addstep, hdr, prstep = 'RECIPROCAL' $
                                  , prproc = inam, prpara = prpara
    
      ;; Output gaintable
      file_mkdir, file_dirname(gainname[ifile])
      print, inam+' : Saving '+gainname[ifile]
      ;;fzwrite, float(gain), outdir+namout, ' '
      overwrite = 1
      red_writedata, gainname[ifile], float(gain), header = hdr,$
                     filetype='fits', overwrite = overwrite
      
    endfor                      ; ifile
  endfor                        ; iscav
  
end
