; docformat = 'rst'

;+
;   Find the transformation matrices that maps the reference channel onto the other ones,
;   use the transform to calculated the offset files needed for MOMFBD alignment.
;
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :author:
; 
;     Tomas Hillberg, Institute for Solar Physics, 2015
;
; 
; :returns:
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
; :history:
;
;   2015-12-01 : New implementation that uses OpenCV functionality for finding
;                and aligning the pinholes. The offset files are now directly
;                computed instead of fitted.
;
;-

function make_corners, clip
    dim = size(clip,/dim)
    corners = intarr(4,3)
    if dim(0) gt 3 then begin
        corners(*,2) = 1
        corners([0,2], 0) = clip(0)
        corners([1,3], 0) = clip(1)
        corners([0,1], 1) = clip(2)
        corners([2,3], 1) = clip(3)
    endif
    return, corners
end


PRO red::pinholecalib_thi, threshold = threshold $
                         , pref = pref $
                         , extraclip = extraclip $
                         , refcam = refcam $
                         , verbose = verbose

    ;; Name of this method
    inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])
  
    ;; Logging
    help, /obj, self, output = selfinfo 
    red_writelog, selfinfo = selfinfo, logfile = logfile

    if(n_elements(threshold) eq 0) then threshold = 0.1
    if(n_elements(refcam) eq 0) THEN refcam = 0
    if(n_elements(verbose) eq 0) THEN verbose = 0
   
    Ncams = 3

    ;; Search summed pinh images and camtag
    if(n_elements(extraclip) eq 0) then extraclip = [0L, 0L, 0L, 0L]
    if(n_elements(extraclip) eq 1) then extraclip = replicate(extraclip, 4)
    if(n_elements(extraclip) eq 2) then extraclip = [replicate(extraclip[0],2),replicate(extraclip[1],2)]

    ph_dir = self.out_dir+'/pinh_align/'
    ph_dir = red_strreplace(ph_dir,'//','/')
    
    output_dir = self.out_dir+'/calib_thi/'
    output_dir = red_strreplace(output_dir,'//','/')
    file_mkdir, output_dir
    
    self -> getcamtags, dir = self.pinh_dir
    
    cams = [ self.camwbtag, self.camttag, self.camrtag ]    ;  TODO autodetect cameras
    Ncams = n_elements(cams)
    if Ncams lt 2 then begin
        print, inam, ' : Need at least 2 cameras.'
        red_writelog, /add, logfile = logfile, top_info_strings = ' : Need at least 2 cameras.'
        return
    endif

    ;; Selected prefilter or all prefilters?
    fw = file_search( ph_dir + self.camwbtag +'.*.pinh', count = cw)
    if cw eq 0 then begin
        print, inam, ' : ERROR : No wideband pinholes found in ', ph_dir
        red_writelog, /add, logfile = logfile, top_info_strings = ' : ERROR : No wideband pinholes found in ' + ph_dir
        retall
    endif
    
    prefilters = strarr(cw)
    for ii = 0, cw-1 do begin
        ;; Remove directory
        prefilters[ii] = (strsplit(fw[ii], '/', /extract, count=nn))[nn-1]
        ;; Get the prefilter
        prefilters[ii] = (strsplit(prefilters[ii], '.', /extract))[1]
    endfor
    prefilters = prefilters[uniq(prefilters, sort(prefilters))]

    if n_elements(pref) ne 0 then begin
        indx = where(prefilters eq pref)
        if max(indx) eq -1 then begin
            print, inam+' : WARNING : Keyword pref does not match any pinhole file names: ', pref
            red_writelog, /add, logfile = logfile, top_info_strings = 'WARNING : Keyword pref does not match any pinhole file names: ' + pref
            return
        endif
        prefilters = prefilters(indx)
    endif
    Npref = n_elements(prefilters)
  
    tmp = f0(fw[0])
    dim = size(tmp,/dim)

    corners = make_corners( [ extraclip(0), dim[0]-extraclip(1)-1, extraclip(2), dim[1]-extraclip(3)-1] )
  
    h_init = fltarr(3, 3, Ncams)

    for ipref = 0, Npref-1 do begin

        files = [ [file_search(ph_dir + cams[0] +'.'+prefilters[ipref]+'*.pinh', count = cw)], $
                  [file_search(ph_dir + cams[1] +'.'+prefilters[ipref]+'*.pinh', count = ct)], $
                  [file_search(ph_dir + cams[2] +'.'+prefilters[ipref]+'*.pinh', count = cr)]]


        if (cw ne ct) or (cw ne cr) then begin
            print, inam+' : Mismatch in available states for the different cameras.'
            print, '  Number of states for cameras WB, NBT, NBR:', cw, ct, cr
            continue
        endif

        align = fltarr(3, 3, cw, Ncams)
        clips = fltarr(cw, 4)
     
        for istate = 0, cw-1 do begin
         
            print, inam+' : calibrating channels:'
            common_fov = corners
            ref_img = f0(files[istate,refcam])
            for icam = 0, Ncams-1 do begin
                if icam EQ refcam then begin
                    print, ' -> '+files[istate,refcam] + ' (reference)'
                endif else begin
                    print, ' -> '+files[istate,icam]
                    img = f0(files[istate,icam])
                    if max(h_init(*,*,icam)) gt 0 then begin
                        this_init = h_init(*,*,icam)
                    endif
                    this_transform = rdx_img_align( ref_img, img, nref=4, h_init=this_init, threshold=threshold, verbose=verbose )

                    align(*,*,istate,icam) = this_transform
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
            endfor

            clips(istate,0) = max(common_fov([0,2],0))
            clips(istate,1) = min(common_fov([1,3],0))
            clips(istate,2) = max(common_fov([0,1],1))
            clips(istate,3) = min(common_fov([2,3],1))

        endfor                         ; istate

        ; average over tunings
        align_avg = total(align,3)/cw
        

        cl = intarr(4,Ncams)
        cl(*,refcam) = [ ceil( max(clips(*,0))), $
                          floor(min(clips(*,1))), $
                          ceil( max(clips(*,2))), $
                          floor(min(clips(*,3)))] + 1   ;   N.B. align_clip index is 1-based

        sx = long(max(cl(0:1,refcam)) - min(cl(0:1,refcam)) + 1) 
        sy = long(max(cl(2:3,refcam)) - min(cl(2:3,refcam)) + 1)

        ref_origin = [ min(cl(0:1,refcam))-1, min(cl(2:3,refcam))-1, 1 ]

       
        indices = [ [[dindgen(sx)#replicate(1.d0, sy) + ref_origin(0)]], $
                    [[replicate(1.d0, sx)#dindgen(sy) + ref_origin(1)]], $
                    [[replicate(1.d0, sx, sy)]]]

        print, inam+' : generating offset files:'
        center = transpose([ (cl(0,refcam)+cl(1,refcam))/2, (cl(2,refcam)+cl(3,refcam))/2, 1 ])
        for icam = 0, Ncams-1 do begin
            if icam NE refcam then begin
            
                ; for the non-refernce channels, map the reference center and cutout a region of the right size
                swap = transpose([0,0,1]) # reform(align_avg(*,*,icam)) gt center
                tmp_center = fix(center # reform(align_avg(*,*,icam)))
                cl(*,icam) = [ tmp_center(0)-(1-2*swap(0))*abs(center(0)-cl(swap(0),refcam)), $
                               tmp_center(0)+(1-2*swap(0))*abs(center(0)-cl(1-swap(0),refcam)), $
                               tmp_center(1)-(1-2*swap(1))*abs(center(1)-cl(2+swap(1),refcam)), $
                               tmp_center(1)+(1-2*swap(1))*abs(center(1)-cl(3-swap(1),refcam)) ]
                               
                ; generate offset files
                for istate = 0, cw-1 do begin
                
                    offs = reform(indices, sx*sy, 3) # align(*,*,istate,icam)
                    idx = where(offs(*,2) ne 0, COMPLEMENT=idx_c)
                    if max(idx) ne -1 then begin
                        offs(idx,0) /= offs(idx,2)
                        offs(idx,1) /= offs(idx,2)
                    endif
                    if max(idx_c) ne -1 then offs(idx_c,0:1) = 0
                    offs = reform(offs, sx, sy, 3) 

                    if align(0,0,istate,icam) lt 0 then offs = reverse(offs,1)
                    if align(1,1,istate,icam) lt 0 then offs = reverse(offs,2)

                    offs -= indices
                    
                    offs(*,*,0) -= (min(cl(0:1,icam)) - ref_origin(0) - 1)
                    offs(*,*,1) -= (min(cl(2:3,icam)) - ref_origin(1) - 1)
                    
                    if (cl(1,icam)-cl(0,icam)) lt 0 then begin
                        offs = reverse(offs,1)
                        offs(*,*,0) *= -1
                    endif
                    
                    if (cl(3,icam)-cl(2,icam)) lt 0 then begin
                        offs = reverse(offs,2)
                        offs(*,*,1) *= -1
                    endif

                    if align(0,0,istate,icam)*(cl(1,icam)-cl(0,icam)) lt 0 OR $
                       align(1,1,istate,icam)*(cl(3,icam)-cl(2,icam)) lt 0 then begin
                        print, inam + ' : sanity check failed - the sign of the transform does not match the align-clip.'
                        print, align(*,*,istate,icam)
                        print, cl(*,icam)
                    endif
    
                    this_file = output_dir + file_basename(files[istate,icam],'.pinh')
                    print, ' -> ' + this_file +'.(x|y)offs'
                    fzwrite, fix(round(100*offs(*,*,0))), this_file + '.xoffs'
                    fzwrite, fix(round(100*offs(*,*,1))), this_file + '.yoffs'
                    
                endfor
            endif
        endfor


        ; Save align clips etc.
        acl = 'ALIGN_CLIP='+strjoin(strtrim(cl, 2), ',')
        openw, lun, output_dir+'align_clips.'+prefilters[ipref]+'.txt', /get_lun
        for icam = 0L, Ncams-1 do begin
            printf, lun, acl[icam]
            print, '  -> '+cams[icam]+': '+acl[icam]
        endfor
        free_lun, lun
        
        ssh = round(reform(align(2,0:1,*)))         ;  NB: the shifts calulated with the old routines have different meaning
        refrot = 0                                  ;  TODO: extract rotation from transformation matrix

        save, file = output_dir + 'align_clips.'+prefilters[ipref]+'.sav', align, acl, cl, refrot, sx, sy, ssh


    endfor                        ; ipref
  
end
