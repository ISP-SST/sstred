; docformat = 'rst'

;+
; Return the filename for the flatfield image for a given (camera-)channel.
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
;     flatname : type=string
; 
; 
; :Params:
; 
;     cam : in, type=string
;   
;         Name of the camera channel. Ex. 'Crisp-W'
;   
; 
; :Keywords:
; 
;     state : in, optional, type=string
;   
;         State string.
;   
;     tag : in, optional, type=boolean
;
;         Interpret cam as a "tag" (Ex. 'camXX'), i.e. do not
;         bother looking up the camtag.
; 
;     data : out, optional
;
;         Return the flat image.
; 
;     header : out, optional
;
;         Return the flat header.
; 
;     summed_name : out, optional, type=string
;
;         Return the filename for the summed (un-normalized) flat.
; 
;     summed_data : out, optional
;
;         Return the summed flat image.
; 
;     summed_header : out, optional
;
;         Return the summed flat header.
; 
; 
; :history:
; 
;   2016-05-19 : First version.
; 
; 
;-
function red::getflat, cam, $
                         state = state, $
                         tag = tag, $
                         data = data, $
                         header = header, $
                         summed_name = summed_name, $
                         summed_data = summed_data, $
                         summed_header = summed_header
    
    if ~keyword_set(tag) then begin
        camtag = self->RED::getcamtag(cam)
    endif else camtag = cam
    
    ; append state
    if( n_elements(state) eq 1 ) then  camtag += '.' + state
    
    filename = self.out_dir+'/flats/'+camtag + '.flat'
    summed_name = self.out_dir+'/flats/summed/'+camtag + '.summed.flat'

    if file_test(filename) then begin
        if arg_present(data) then begin
            data = red_readdata(filename, header=header)
        endif else if arg_present(header) then begin
            header = red_readhead(filename)
        endif
    endif

    if file_test(summed_name) then begin
        if arg_present(summed_data) then begin
            summed_data = red_readdata(summed_name, header=summed_header)
        endif else if arg_present(summed_header) then begin
            summed_header = red_readhead(summed_name)
        endif
    endif

    return, filename
  
end
