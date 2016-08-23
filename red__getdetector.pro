; docformat = 'rst'

;+
; Return the camera-tag for a given (camera-)channel.
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :author:
; 
;     Tomas Hillberg, ISP
; 
; 
; :returns:
;
;    camera-tag : type=string
; 
; 
; :Params:
; 
;    camera : in, type=string
;   
;      Name of the camera. Ex. 'Crisp-W'
;   
; 
; :Keywords:
;
; 
; :history:
; 
;   2016-05-19 : First version.
;
;   2016-08-23 : THI. Rename camtag to detector and channel to camera,
;                so the names match those of the corresponding SolarNet
;                keywords.
; 
;-
function red::getdetector, camera
    
    ; No tags, try to load or parse filenames
    if ~ptr_valid(self.detectors) then begin
        self->getdetectors
    endif
    
    ; Still no camera- or detectorlists, return empty result
    if( ~ptr_valid(self.cameras) || ~ptr_valid(self.detectors) ) then begin
        return, ''
    endif
    
    ; Sizes does not match, return empty result
    if( n_elements(*self.cameras) ne n_elements(*self.detectors) ) then begin
        return, ''
    endif
    
    pos = where(*self.cameras eq camera)
    
    ; No match or multiple matches.
    if( n_elements(pos) ne 1 || max(pos) lt 0 ) then begin
        return, ''
    endif
    
    return, (*self.detectors)[pos[0]]

end
