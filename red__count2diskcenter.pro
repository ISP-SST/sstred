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
;    pref : 
;   
;   
;   
; 
; :Keywords:
; 
;    cam  : 
;   
;   
;   
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2016-02-15 : MGL. Use loadbackscatter.
;
;   2016-08-23 : THI. Rename camtag to detector and channel to camera,
;                so the names match those of the corresponding SolarNet
;                keywords.
; 
; 
;-
function red::count2diskcenter, pref, cam = cam
  inam = ' red::count2diskcenter :'
  self->getdetectors, dir = self.data_dir
                                
  fc1 = self.out_dir + '/prefilter_fits/'+self.camttag+'.'+pref+'.prefilter_pars.f0'
  fc2 = self.out_dir + '/prefilter_fits/'+self.camrtag+'.'+pref+'.prefilter_pars.f0'
                                
  if(file_test(fc1) AND file_test(fc2)) then begin
     i0 = (f0(fc1))[0]
     i1 = (f0(fc2))[0]
     aver = 0.5 * (i0 + i1)
  endif else begin
     print, inam + 'ERROR, Could not find prefilter fit parameters'
     stop
  endelse
  print, inam + 'Transmitted -> ', i0
  print, inam + 'Reflected   -> ', i1
  print, inam + 'Average     -> ', aver

  ;; load WB flat?

  if(~keyword_set(cam)) then cam = self.camwbtag

  fzread,ff, self.out_dir+'/flats/'+cam+'.'+pref+'.flat', h

  if((pref EQ '8542') AND (self.descatter_dir NE '')) then begin
    self -> loadbackscatter, cam, pref, bg, psf
    ff = rdx_descatter(ff, bg, psf)
  endif
  
  ;; Now, we know the variation of the WB images on time and the ratio
  ;; between WB and NW at the time of the flats: We can use ratios to
  ;; obtain the correction factor for the NB.

  var = aver / median(ff)

  return, var
end
