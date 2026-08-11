pro red_make_offs, mapping, xoffs, yoffs, align_clip, ref_clip=ref_clip

    
    if(n_elements(align_clip) ne 4) then align_clip = [1L, 1024L, 1L, 1024L]
    if(n_elements(ref_clip) ne 4) then ref_clip = align_clip
    
    if mapping[0,0]*(ref_clip[1]-ref_clip[0])*(align_clip[1]-align_clip[0]) lt 0 then begin
        align_clip[0:1] = reverse(align_clip[0:1])
    endif
    if mapping[1,1]*(ref_clip[3]-ref_clip[2])*(align_clip[3]-align_clip[2]) lt 0 then begin
        align_clip[2:3] = reverse(align_clip[2:3])
    endif
    
    nx = long(max(align_clip[0:1]) - min(align_clip[0:1]) + 1)
    ny = long(max(align_clip[2:3]) - min(align_clip[2:3]) + 1)

    ref_origin = [ min(ref_clip[0:1])-1, min(ref_clip[2:3])-1 ]

    xoffs = fltarr(nx,ny)
    yoffs = fltarr(nx,ny)

    indices = [ [[dindgen(nx)#replicate(1.d0, ny) + ref_origin[0]]], $
                [[replicate(1.d0, nx)#dindgen(ny) + ref_origin[1]]], $
                [[replicate(1.d0, nx, ny)]]]

    offs = reform(indices, nx*ny, 3) # mapping

    idx = where(offs(*,2) ne 0, COMPLEMENT=idx_c)
    if max(idx) ne -1 then begin
        offs[idx,0] /= offs[idx,2]
        offs[idx,1] /= offs[idx,2]
    endif
    if max(idx_c) ne -1 then offs[idx_c,0:1] = 0
    offs = reform(offs, nx, ny, 3)

    if mapping[0,0] lt 0 then begin
        offs = reverse(offs,1)
    endif
    if mapping[1,1] lt 0 then begin
        offs = reverse(offs,2)
    endif

    offs -= indices

    offs[*,*,0] -= (min(align_clip[0:1]) - ref_origin[0] - 1)
    offs[*,*,1] -= (min(align_clip[2:3]) - ref_origin[1] - 1)
    
    if (align_clip[1]-align_clip[0]) lt 0 then begin
        offs = reverse(offs,1)
        offs[*,*,0] *= -1
    endif
    
    if (align_clip[3]-align_clip[2]) lt 0 then begin
        offs = reverse(offs,2)
        offs[*,*,1] *= -1
    endif
   
    xoffs = fix(round(100*offs[*,*,0]))
    yoffs = fix(round(100*offs[*,*,1]))
        
end
