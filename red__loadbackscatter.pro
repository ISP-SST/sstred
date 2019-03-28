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
;     2019-03-28 : MGL. Use the no_load keyword when calling
;                  red_loadbackscatter. 
;
;-
pro red::loadbackscatter, cam, pref, bgain, bpsf $
                          , bgfile = bgfile $
                          , bpfile = bpfile $
                          , write = write

  if arg_present(bgain) or arg_present(bpsf) then begin
    red_loadbackscatter, cam, self.isodate, self.descatter_dir, pref, bgain, bpsf, $
                         bgfile = bgfile , bpfile = bpfile, write = write
  endif else begin
    red_loadbackscatter, cam, self.isodate, self.descatter_dir, pref, /no_load, $
                         bgfile = bgfile , bpfile = bpfile
  endelse 

end
