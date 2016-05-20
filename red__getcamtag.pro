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
;    channel : in, type=string
;   
;      Name of the channel. Ex. 'Crisp-W'
;   
; 
; :Keywords:
;
; 
; :history:
; 
;   2016-05-19 : First version.
; 
;-
function red::getcamtag, channel
    
    ; No tags, try to load or parse filenames
    if ~ptr_valid(self.cam_tags) then begin
        self->getcamtags
    endif
    
    ; Still no tags, or no channel-list.
    if( ~ptr_valid(self.cam_channels) || ~ptr_valid(self.cam_tags) ) then begin
        return, ''
    endif
    
    ; Sizes does not match.
    if( n_elements(*self.cam_channels) ne n_elements(*self.cam_tags) ) then begin
        return, ''
    endif
    
    pos = where(*self.cam_channels eq channel)
    
    ; No match or multiple matches.
    if( n_elements(pos) ne 1 || max(pos) lt 0 ) then begin
        return, ''
    endif
    
    return, (*self.cam_tags)[pos[0]]

end
