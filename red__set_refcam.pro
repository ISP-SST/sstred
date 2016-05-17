PRO red::set_refcam, cam

    self->getcamtags        ; load tags
    IF size(cam,/type) EQ 7 THEN BEGIN
        idx = -1
        IF ptr_valid(self.cam_tags) THEN idx = where(*self.cam_tags EQ cam)
        IF idx lt 0 THEN BEGIN
            print,'set_refcam: Could not find ' + cam + ' in camera list.'
            idx = 0
        ENDIF
    ENDIF ELSE BEGIN
        idx = byte(cam)
        IF idx ge n_elements(*self.cam_tags) THEN BEGIN
            print,'set_refcam: Index ', idx, ' out of bounds.'
            idx = 0
        ENDIF
    ENDELSE
    self.refcam = idx

END 
