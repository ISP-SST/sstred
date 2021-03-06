; docformat = 'rst'

;+
; Download backgain descatter data if needed, then return the gain and
; psf in the parameters.
; 
; :Categories:
;
;    SST observations
; 
; 
; :Author:
; 
;    Mats Löfdahl, 2016-02-15
; 
; 
; :Params:
; 
;    cam : in, type=string
;
;       The camera tag.
;
;    pref : in, type=string
;
;       The prefilter.
;
;    bgain : out, type=fltarr
;
;       The backscatter gain.
;
;    bpsf : out, type=fltarr
;
;       The backscatter psf.
;
;
; :Keywords: 
;
;    bgfile : out, optional, type=string
;
;       The name of the backscatter gain file, if existing.
;
;
;    bpfile : out, optional, type=string
;
;       The name of the backscatter psf file, if existing.
;
;    no_load : in, optional, type=boolean
;
;       Do not load the backscatter data.
;
;    write : in, optional, type=boolean
;
;       Normally, loadbackscatter reads the files. With /write, it
;       writes instead the backgain and psf provided in bgain and
;       bpsf.
;
;
; :History:
; 
;     2016-02-15 : MGL. First version.
;
;     2016-02-17 : MGL. New keyword "write_gain". Only read files if
;                  parameters bgain (and bpsf) are present, otherwise
;                  just construct the file names.
;
;     2016-02-19 : MGL. Allow to write both gain and psf.
;
;     2019-03-28 : MGL. New keyword no_load.
;
;-
pro red_loadbackscatter, cam, date, dir, pref, bgain, bpsf $
                         , bgfile = bgfile $
                         , bpfile = bpfile $
                         , no_load = no_load $
                         , write = write

  year = (strsplit(date, '.-', /extract))[0]
  
  bgfile = dir + path_sep() + cam + '.backgain.' + pref + '_' + year + '.f0'
  bpfile = dir + path_sep() + cam + '.psf.' + pref + '_' + year + '.f0'
  
  if keyword_set(write) then begin

    ;; Only write the files if possible.

    if size(bgain, /n_dim) lt 1 then begin
      print, 'red_loadbackscatter : Cannot write the provided bgain.'
      help, bgain
      stop
    endif else begin
      fzwrite, bgain, bgfile, ' '
    endelse

    if size(bpsf, /n_dim) lt 1 then begin
      print, 'red_loadbackscatter : Cannot write the provided bpsf.'
      help, bpsf
      stop
    endif else begin
      fzwrite, bpsf, bpfile, ' '
    endelse
    
  endif else begin
    
    ;; Construct file names

    if arg_present(bgain) or arg_present(bpsf) then begin
      if ~file_test(bgfile) or ~file_test(bpfile) then red_download, date=date, backscatter = pref
    endif
    
    if arg_present(bgain) then begin
      
      ;; Read the gain if wanted
      
      if file_test(bgfile) then begin
        
        print, 'red_loadbackscatter : Loading backscatter gain for ' + cam + ', ' + pref
        bgain = f0(bgfile)
        
      endif else begin
        
        print, 'red_loadbackscatter : Backscatter gain not available for ' + cam + ', ' + pref
        print, bgfile
        stop
        
      endelse                   ; exists
    endif                       ; wants gain

    
    if arg_present(bpsf) then begin

      ;; Read the psf if wanted

      if file_test(bpfile) then begin
        
        print, 'red_loadbackscatter : Loading backscatter psf for ' + cam + ', ' + pref
        bpsf  = f0(bpfile)
        
      endif else begin
        
        print, 'red_loadbackscatter : Backscatter psf not available for ' + cam + ', ' + pref
        print, bpfile
        stop
        
      endelse                   ; exists
    endif                       ; wants psf

  endelse                  

end
