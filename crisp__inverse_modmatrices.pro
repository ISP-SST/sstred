; docformat = 'rst'

;+
; Read inverse modulation matrix maps from disk, make and store them
; if necessary.
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, Institute for Solar Physics
; 
; 
; :Returns:
; 
; 
; :Params:
; 
; 
; 
; 
; 
; 
; :Keywords:
; 
;    smooth : in, optional, type=integer, default=9
; 
;       Width in pixels of smoothing kernel applied to modulation
;       matrix maps.
; 
;   x01y01 : in, optional, type=intarr[4]
;   
;      Define clip area [x0,x1,y0,y1] -> [x0:x1,y0:y1].
; 
; 
; :History:
; 
; 
; 
; 
; 
; 
;-
pro crisp::inverse_modmatrices, prefilter, dir $
                                , camr = camr $
                                , camt = camt $
                                , immr = immr $
                                , immt = immt $
                                , no_ccdtabs = no_ccdtabs $
                                , overwrite = overwrite $
                                , smooth = smooth $
                                , x01y01 = x01y01
  
  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  if n_elements(smooth) eq 0 then smooth = 9

  if n_elements(camr) gt 0 then begin
    red_append, cams, camr
    red_append, prefixes, 'r'
  endif
  
  if n_elements(camt) gt 0 then begin
    red_append, cams, camt
    red_append, prefixes, 't'
  endif

  Ncams = n_elements(cams)

  for icam = 0, Ncams-1 do begin

    detector = self->getdetector( cams[icam] )

    fname = dir + '/imm'+prefixes[icam]+'.fits'

    if file_test(fname) and ~keyword_set(overwrite) then begin

      imm = red_readdata(fname, head = hdr)
      ;; Get prpara info from hdr and check if they match current
      ;; parameters?
      
    endif else begin

      ;; Geometrical distortion map
      self -> getalignment, align=align, prefilters=prefilter
      indx = where(align.state2.camera eq cams[icam], Nalign)
      case Nalign of
        0    : stop             ; Should not happen!
        1    : amap = invert(      align[indx].map           )
        else : amap = invert( mean(align[indx].map, dim = 3) )
      endcase

      ;; Polcal output
      pname = self.out_dir+'/polcal/'+detector+'_'+prefilter+'_polcal.fits'
      
      if file_test(pname) then begin

        mm = red_readdata(pname) ; Modulation matrix
        dims = size(mm, /dim)
        Nelements = dims[0]
        if n_elements(x01y01) eq 4 then begin
          ;; The size of the momfbd-restored images
          Nx = x01y01[1]-x01y01[0]+1
          Ny = x01y01[3]-x01y01[2]+1
        endif else begin
          Nx = dims[1]
          Ny = dims[2]
        endelse
        
        ;; Interpolate Sarnoff CCD tabs?
        if strmatch((red_camerainfo(detector)).model,'Sarnoff*') then begin
          if ~keyword_set(no_ccdtabs) then begin
            for ii = 0, Nelements-1 do mm[ii,*,*] = red_mask_ccd_tabs(reform(mm[ii,*,*]))
          endif
        endif

        ;; Check NaNs
        for ii = 0, Nelements-1 do begin
          mask = finite(reform(mm[ii,*,*]))
          idx = where(mask, count)
          if count gt 0 and count lt Nx*Ny then $
             mm[ii,*,*] = red_fillpix(reform(mm[ii,*,*]), mask=mask)
        endfor                  ; ii

        ;; Smooth?
        if smooth gt 0 then begin
          dpix = round(smooth)*3
          if (dpix/2)*2 eq dpix then dpix -= 1
          dpsf = double(smooth)
          psf = red_get_psf(dpix, dpix, dpsf, dpsf)
          psf /= total(psf, /double)
          for ii=0,Nelements-1 do mm[ii,*,*] = red_convolve(reform(mm[ii,*,*]), psf)
        endif

        print,'Transforming and clipping transmitted modulation matrix to ' $
              + red_stri(Nx) + 'x' + red_stri(Ny)

        tmp = make_array( Nelements, Nx, Ny, type=size(mm, /type) )
        for ii=0,Nelements-1 do begin
          tmpp = rdx_img_project(amap, reform(mm[ii,*,*])) ; Apply the geometrical mapping
          if n_elements(x01y01) eq 4 then begin
            tmp[ii,*,*] = tmpp[x01y01[0]:x01y01[1],x01y01[2]:x01y01[3]] ; Clip to the selected FOV
          endif else begin
            dims = size(tmpp, /dim)
            case 1 of           ; Adjust array size after rdx_img_project. (But does this position it correctly?)
              dims[0] ge Nx and dims[1] ge Ny : tmp[ii, *, *] = tmpp[0:Nx-1, 0:Ny-1]
              dims[0] ge Nx : tmp[ii, *, 0:dims[1]-1] = tmpp[0:Nx-1, *]
              dims[1] ge Ny : tmp[ii, 0:dims[0]-1, *] = tmpp[*, 0:Ny-1]
              else : tmp[ii, 0:dims[0]-1, 0:dims[1]-1] = tmpp
            endcase
          endelse
        endfor                  ; ii
        mm = temporary(tmp)
        imm = red_invert_mmatrix(temporary(mm)) ; Inverse modulation matrix
        
      endif else begin
        print, inam + ' : ERROR, polcal data not found in ' + self.out_dir + '/polcal/'
        return
      endelse

      ;; Store some prpara info in the header? 
      red_writedata, fname, imm
      
    endelse

    case prefixes[icam] of
      'r'  : immr = imm
      't'  : immt = imm
      else : stop               ; Should not happen!
    endcase
    
  endfor                        ; icam

end
