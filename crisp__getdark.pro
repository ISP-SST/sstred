; docformat = 'rst'

;+
; Return the filename for the darkfield image for a given (camera-)channel.
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
;     darkname : type=string
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
;     tag : in, optional, type=boolean
;
;         Interpret cam as a "tag" (Ex. 'camXX'), i.e. do not
;         bother looking up the camtag.
; 
;     data : out, optional
;
;         Return the dark image.
; 
;     header : out, optional
;
;         Return the dark header.
; 
;     summed_name : out, optional, type=string
;
;         Return the filename for the summed (un-normalized) dark.
; 
;     summed_data : out, optional
;
;         Return the summed dark image.
; 
;     summed_header : out, optional
;
;         Return the summed dark header.
; 
; 
; :history:
; 
;   2016-05-19 : First version.
; 
; 
;-
function crisp::getdark, cam, $
                         tag = tag, $
                         data = data, $
                         header = header, $
                         summed_name = summed_name, $
                         summed_data = summed_data, $
                         summed_header = summed_header
    
    if ~keyword_set(tag) then begin
        camtag = self->RED::getcamtag(cam)
    endif else camtag = cam
    
    filename = self.out_dir+'/darks/'+camtag+'.dark'
    summed_name = self.out_dir+'/darks/'+camtag+'.summed.0000001'

    if file_test(filename) then begin
        if arg_present(data) then begin
            fzread, data, filename, header
        endif else if arg_present(header) then begin
            header = fzhead(filename)
        endif
    endif

    if file_test(summed_name) then begin
        if arg_present(summed_data) then begin
            fzread, summed_data, summed_name, summed_header
        endif else if arg_present(summed_header) then begin
            summed_header = fzhead(summed_name)
        endif
    endif

    return, filename
  
end
