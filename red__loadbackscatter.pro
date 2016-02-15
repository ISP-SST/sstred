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
; :History:
; 
;     2016-02-15 : MGL. First version.
;
;
;
;-
pro red::loadbackscatter, cam, pref, bgain, bpsf

  year = (strsplit(self.isodate, '-', /extract))[0]
  
  bgfile = self.descatter_dir + '/' + cam + '.backgain.'+pref+'_'+year+'.f0'
  bpfile = self.descatter_dir + '/' + cam + '.psf.'+pref+'.f0'
  
  if ~file_test(bgfile) or ~file_test(bpfile) then self -> download, backscatter = pref

  if file_test(bgfile) and file_test(bpfile) then begin
  
     print, 'Loading backscatter data for ' + cam + ', ' + pref
     bgain = f0(bgfile)
     bpsf  = f0(bpfile)

  endif else begin

     print, 'Backscatter data not available for ' + cam + ', ' + pref
     stop
     
  endelse
end
