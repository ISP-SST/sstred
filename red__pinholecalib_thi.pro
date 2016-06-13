; docformat = 'rst'

;+
;   Find the transformation matrices that maps the reference channel
;   onto the other ones, use the transform to calculated the offset
;   files needed for MOMFBD alignment.
;
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;     Tomas Hillberg, Institute for Solar Physics, 2015
;
; 
; 
; :Params:
; 
;   
; 
; :Keywords:
;    
;    threshold : in, optional, type=float, default=0.1
; 
;       Threshold for identifying a strong enough pinhole.
;
;    
;    extraclip : 
;
;
;    refcam : in, optional, type=integer, default=0
;    
;      Select reference-channel.
;
;    
;    verbose : in, optional, type=integer, default=0
;    
;      Provide more output
;
;    
;    pref : in, optional, type=string
;     
;      Indicate the prefilter you want to calculate the clips for,
;      Default is to do it for all prefilters there is data for.
;
; 
; :History:
;
;   2015-12-01 : New implementation that uses OpenCV functionality for finding
;                and aligning the pinholes. The offset files are now directly
;                computed instead of fitted.
;
;   2016-06-13 : MGL. Various cosmetic edits. CHange bunch of if
;                statements to a case statement. Move make_corners to
;                the red_ namespace and its own file.
;
;-
pro red::pinholecalib_thi, threshold = threshold $
                         , nref = nref $
                         , pref = pref $
                         , extraclip = extraclip $
                         , dir = dir $
                         , cams = cams $
                         , refcam = refcam $
                         , verbose = verbose

    ;; Name of this method
    inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])
  
    ;; Logging
    help, /obj, self, output = selfinfo 
    red_writelog, selfinfo = selfinfo, logfile = logfile

    if(n_elements(threshold) eq 0) then threshold = 0.3
    if(n_elements(nref) eq 0) then nref = 4
    if( n_elements(dir) gt 0 ) then dir = [dir] $
    else if ptr_valid(self.pinh_dirs) then dir = *self.pinh_dirs
    self -> getcamtags, dir = dir
    if( n_elements(cams) gt 0 ) then cams = [cams] $
    else if ptr_valid(self.cam_tags) then cams = [*self.cam_tags]

    if(n_elements(refcam) eq 0) THEN refcam = self.refcam
    if(n_elements(verbose) eq 0) THEN verbose = 0
   
    Ncams = n_elements(cams)

    case n_elements(extraclip) of
       0 : extraclip = [0L, 0L, 0L, 0L]
       1 : extraclip = replicate(extraclip, 4)
       2 : extraclip = [ replicate(extraclip[0], 2), $
                         replicate(extraclip[1], 2) ]
       4 :                      ; Leave as it is.
       else : begin
          inam + "ERROR: Don't know how to use keyword extraclip with " $
             + strtrim(n_elements(extraclip), 2) + ' elements.'
          stop
       end
    endcase

    ;; Search summed pinh images and camtag
    ph_dir = self.out_dir+'/pinhs/'
    ph_dir = red_strreplace(ph_dir,'//','/')
    
    output_dir = self.out_dir+'/calib/'
    output_dir = red_strreplace(output_dir,'//','/')
    file_mkdir, output_dir
    
    if Ncams lt 2 then begin
        print, inam, ' : Need at least 2 cameras.'
        print, inam, ' : dir: ',dir
        print, inam, ' : cams: ',cams
        print, inam, ' : self.cam_tags: ',*self.cam_tags
        red_writelog, /add, logfile = logfile, top_info_strings = ' : Need at least 2 cameras.'
        return
    endif

    if refcam ge Ncams then begin
        print, inam, ' : index of reference camera out of range: ', refcam, ' >= ', Ncams
        return
    endif

    ;; Selected prefilter or all prefilters?
    files = file_search( ph_dir + '*.pinh' )
    self->selectfiles, files=files, prefilter=pref, states=states, /strip_settings
    nf = n_elements( files )
    if nf eq 0 then begin
        print, inam, ' : ERROR : No pinhole files found in ', ph_dir
        red_writelog, /add, logfile = logfile $
                      , top_info_strings = ' : ERROR : No pinhole files found in ' + ph_dir
        retall
    endif
    
    self->selectfiles, files=files, states=states, cam = cams[refcam], selected=ref_sel
    fullstate = [states[ref_sel].fullstate]
    state_list = fullstate[ uniq(fullstate, sort(fullstate)) ]
    ns = n_elements(state_list)
    head = red_readhead( files[0], /silent )
    dims = fxpar(head, 'NAXIS*')

    if( n_elements(dims) ne 2 ) then begin
        print, inam, ' : pinhole data is not 2-dimensional.'
        return
    endif

    corners = red_make_corners( [ extraclip(0), dims[0]-extraclip(1)-1 , $
                                  extraclip(2), dims[1]-extraclip(3)-1 ] )
 
    h_init = fltarr(3, 3, Ncams)

    aligns = fltarr(3, 3, ns, Ncams)
    clips = fltarr(ns, 4)
     
    ;; Loop over states 
    for ss = 0L, ns - 1 do begin

       self->selectfiles, files=files, states=states, cam = cams[refcam] $
                          , ustat=state_list[ss], selected=ref_sel
        if( n_elements(ref_sel) ne 1 || ref_sel lt 0 ) then begin
            print, inam, " : multiple reference images for state='", state_list[ss], "'"
            continue
        endif
   
        common_fov = corners
        ref_img = red_readdata( files[ref_sel], /silent )

        for icam = 0, Ncams-1 do begin
            if icam EQ refcam then begin
                print, ' -> ' + files[ref_sel] + ' (reference)'
            endif else begin
               self->selectfiles, files=files, states=states, cam = cams[icam] $
                                  , ustat=state_list[ss], selected=sel
                if( n_elements(sel) ne 1 || sel lt 0 ) then begin
                   print, inam, " : multiple or missing images for cam='" $
                          , cams[icam],"' state='", state_list[ss], "'"
                    continue
                endif
                print, ' -> ' + files[sel]

                img = red_readdata(files[sel], /silent)
                if max(h_init(*,*,icam)) gt 0 then begin
                    this_init = h_init(*,*,icam)
                endif

                this_transform = rdx_img_align( ref_img, img, nref=nref, h_init=this_init $
                                                , threshold=threshold, verbose=verbose )

                aligns(*,*,ss,icam) = this_transform
                h_init(*,*,icam) = temporary(this_init)

                ; transform the corners of the FOV to the non-reference camera
                common_fov = common_fov # this_transform
                idx = where(common_fov(*,2) ne 0, COMPLEMENT=idx_c)
                if max(idx) ne -1 then begin
                    common_fov(idx,0) /= common_fov(idx,2)
                    common_fov(idx,1) /= common_fov(idx,2)
                endif
                if max(idx_c) ne -1 then common_fov(idx_c,*) = 0
                common_fov(*,2) = 1
                
                ; if a mapped corner falls outside the range, clip it
                idx = where(common_fov(*,0) lt corners(0,0) )
                if max(idx) ne -1 then  common_fov(idx,0) = corners(0,0)
                idx = where(common_fov(*,0) gt corners(1,0) )
                if max(idx) ne -1 then  common_fov(idx,0) = corners(1,0)
                idx = where(common_fov(*,1) lt corners(0,1) )
                if max(idx) ne -1 then  common_fov(idx,1) = corners(0,1)
                idx = where(common_fov(*,1) gt  corners(2,1) )
                if max(idx) ne -1 then  common_fov(idx,1) = corners(2,1)
                
                ; transform back to reference camera
                common_fov = common_fov # invert(this_transform)
                idx = where(common_fov(*,2) ne 0, COMPLEMENT=idx_c)
                if max(idx) ne -1 then begin
                    common_fov(idx,0) /= common_fov(idx,2)
                    common_fov(idx,1) /= common_fov(idx,2)
                endif
                if max(idx_c) ne -1 then common_fov(idx_c,*) = 0
                common_fov(*,2) = 1

                ; clip again
                idx = where( common_fov(*,0) lt corners(0,0) )
                if max(idx) ne -1 then  common_fov(idx,0) = corners(0,0)
                idx = where( common_fov(*,0) gt corners(1,0) )
                if max(idx) ne -1 then  common_fov(idx,0) = corners(1,0)
                idx = where( common_fov(*,1) lt corners(0,1) )
                if max(idx) ne -1 then  common_fov(idx,1) = corners(0,1)
                idx = where( common_fov(*,1) gt corners(2,1) )
                if max(idx) ne -1 then  common_fov(idx,1) = corners(2,1)
                
            endelse
        endfor      ; icam

        clips(ss,0) = max(common_fov([0,2],0))
        clips(ss,1) = min(common_fov([1,3],0))
        clips(ss,2) = max(common_fov([0,1],1))
        clips(ss,3) = min(common_fov([2,3],1))
        
        cl = intarr(4,Ncams)
        cl(*,refcam) = [ ceil( max(clips(*,0))), $
                         floor(min(clips(*,1))), $
                         ceil( max(clips(*,2))), $
                         floor(min(clips(*,3))) ] + 1 ;   N.B. align_clip index is 1-based

        sx = long(max(cl(0:1,refcam)) - min(cl(0:1,refcam)) + 1) 
        sy = long(max(cl(2:3,refcam)) - min(cl(2:3,refcam)) + 1)

        ref_origin = [ min(cl(0:1,refcam))-1, min(cl(2:3,refcam))-1 ]

        indices = [ [[dindgen(sx)#replicate(1.d0, sy) + ref_origin(0)]], $
                    [[replicate(1.d0, sx)#dindgen(sy) + ref_origin(1)]], $
                    [[replicate(1.d0, sx, sy)]]]

        print, inam+' : generating offset files:'
        center = transpose([ (cl(0,refcam)+cl(1,refcam))/2, (cl(2,refcam)+cl(3,refcam))/2, 1 ])

        for icam = 0, Ncams-1 do begin
            if icam NE refcam then begin
            
                cam_origin = [ min(cl(0:1,icam))-1, min(cl(2:3,icam))-1 ]
                self->selectfiles, files=files, states=states, cam = cams[icam] $
                                   , ustat=state_list[ss], selected=sel
                if( n_elements(sel) ne 1 || sel lt 0 ) then begin
                    continue
                endif
                this_file = output_dir + file_basename(files[sel],'.pinh')

                ;; For the non-refernce channels, map the reference
                ;; center and cutout a region of the right size
                
                swap = transpose([0,0,1]) # reform(aligns(*,*,ss,icam)) gt center
                tmp_center = fix(center # reform(aligns(*,*,ss,icam)))
                cl(*,icam) = [ tmp_center(0)-(1-2*swap(0))*abs(center(0)-cl(  swap(0),refcam)), $
                               tmp_center(0)+(1-2*swap(0))*abs(center(0)-cl(1-swap(0),refcam)), $
                               tmp_center(1)-(1-2*swap(1))*abs(center(1)-cl(2+swap(1),refcam)), $
                               tmp_center(1)+(1-2*swap(1))*abs(center(1)-cl(3-swap(1),refcam)) ]

                ; generate offset files
                offs = reform(indices, sx*sy, 3) # aligns(*,*,ss,icam)
                idx = where(offs(*,2) ne 0, COMPLEMENT=idx_c)
                if max(idx) ne -1 then begin
                    offs(idx,0) /= offs(idx,2)
                    offs(idx,1) /= offs(idx,2)
                endif
                if max(idx_c) ne -1 then offs(idx_c,0:1) = 0
                offs = reform(offs, sx, sy, 3) 

                if aligns(0,0,ss,icam) lt 0 then offs = reverse(offs,1)
                if aligns(1,1,ss,icam) lt 0 then offs = reverse(offs,2)

                offs -= indices
                
                offs(*,*,0) -= (cam_origin[0] - ref_origin[0])
                offs(*,*,1) -= (cam_origin[1] - ref_origin[1])
                
                if (cl(1,icam)-cl(0,icam)) lt 0 then begin
                    offs = reverse(offs,1)
                    offs(*,*,0) *= -1
                endif
                
                if (cl(3,icam)-cl(2,icam)) lt 0 then begin
                    offs = reverse(offs,2)
                    offs(*,*,1) *= -1
                endif

                if aligns(0,0,ss,icam)*(cl(1,icam)-cl(0,icam)) lt 0 OR $
                   aligns(1,1,ss,icam)*(cl(3,icam)-cl(2,icam)) lt 0 then begin
                    print, inam + ' : sanity check failed - the sign of the transform does not match the align-clip.'
                    print, aligns(*,*,ss,icam)
                    print, cl(*,icam)
                endif

                print, ' -> ' + this_file +'.(x|y)offs'
                red_writedata, this_file + '.xoffs', fix(round(100*offs(*,*,0))) $
                               , filetype='ANA', /overwrite
                red_writedata, this_file + '.yoffs', fix(round(100*offs(*,*,1))) $
                               , filetype='ANA', /overwrite
                    
            endif
         endfor                 ; icam


        ; Save align clips etc.
        acl = 'ALIGN_CLIP='+strjoin(strtrim(cl, 2), ',')
        clipfile = output_dir+'align_clips.'
        if state_list[ss] ne '' then clipfile += state_list[ss] + '.'
        print, ' -> ' + clipfile + '(txt|sav)'
        openw, lun, clipfile+'txt', /get_lun
        for icam = 0L, Ncams-1 do begin
            printf, lun, acl[icam]
            print, '  -> '+cams[icam]+': '+acl[icam]
        endfor
        free_lun, lun
        
        ssh = round(reform(aligns(2,0:1,*))) ;  NB: the shifts calulated with the old routines have different meaning
        refrot = 0                           ;  TODO: extract rotation from transformation matrix
        align = reform(aligns[*,*,ss,*])

        save, file = clipfile+'sav', align, acl, cl, refrot, sx, sy, ssh


     endfor                     ; ipref
  
end
