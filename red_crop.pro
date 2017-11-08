; docformat = 'rst'

;+
; Crop an image or image cube, out-of-bounds pixels.
;
; A 3D array is interpreted as an image cube where the first two
; dimensions are spatial.
;
; If no or only partial FOV information is given with keywords, use
; the XROI GUI to either modify an initial ROI/FOV or define a new one
; from scratch. The last ROI is used.
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
;    A cropped and possibly padded version of the input array.
; 
; :Params:
; 
;    ims : in, type=array
; 
;      An image or an image cube (up to 8 dimensions, where the first
;      two are the image dimensions) to be cropped.
; 
; 
; :Keywords:
;
;    centered : in, optional, type=boolean
;
;      Set this to get a centered FOV.
;
;    corners : in, out, optional, type=array
;
;      Corner coordinates [llx,lly,urx,ury]. If the FOV is selected
;      specified with other keywords the corners are returned in this
;      keyword.
;
;    nocropping : in, optional, type=boolean
;
;      Don't actually do any cropping, return after establishing
;      the "corners" keyword. 
;
;    pad : in, optional, type="number or string", default=0
;
;      Either a number used to fill in out-of-bounds pixels or a
;      string with a method to calculate such a number from the
;      cropped image ("mean" or "median").
;
;    size : in, optional, type="integer or array"
; 
;      Set this to a FOV size to only select center of FOV
;      interactively. If an array, it is interpreted as sx=size[0],
;      sy=size[1]. If a single integer, sx=sy=size.
; 
;    xc : in, out, optional, type=integer
;   
;      X pixel coordinate of center of FOV. 
;
;    yc : in, out, optional, type=integer
;   
;      Y pixel coordinate of center of FOV. 
; 
; 
; :History:
; 
;    2017-11-06 : MGL. First version.
; 
;    2017-11-07 : MGL. New keyword nocropping. Use the XROI GUI to
;                 define the corners if keywords do not provide enough
;                 information.
; 
;-
function red_crop, ims $
                   , centered = centered $
                   , corners = corners $
                   , nocropping = nocropping $
                   , pad = pad $
                   , size = size $
                   , xc = xc $
                   , yc = yc 

  inam = red_subprogram(/low, calling = inam1)

  ;; Original image dimensions
  dims = size(ims, /dim)

  ;; Interpret the size keyword if given
  case n_elements(size) of
    0 : 
    1 : begin
      Sx = size
      Sy = size          
    end
    else : begin
      Sx = size[0]
      Sy = size[1]
    end
  endcase

  if keyword_set(centered) then begin
    xc = dims[0]/2
    yc = dims[1]/2
  endif

  ;; Make sure we have the corners of the FOV
  case 1 of
    n_elements(size) gt 0 and n_elements(xc) gt 0 and n_elements(yc) gt 0 : begin
      ;; The corners are completely specified by size and center
      ;; coordinates.
      corners = lonarr(4)
      corners[0:1] = [xc-Sx/2, yc-Sy/2]          ; Lower
      corners[2:3] = corners[0:1] + [Sx, Sy] - 1 ; Upper
    end
    n_elements(corners) eq 4 : begin
      ;; Carry on, I guess we can use the corner keyword as input!
    end
    else : begin

      ;; We are lacking some information, use mouse to select ROI from
      ;; first image in the cube. Make a predefined ROI that the user
      ;; can change, using whatever partial info that was given.

      ;; Make it centered as default.
      if n_elements(xc) gt 0 then X_in = xc else X_in = dims[0]/2
      if n_elements(yc) gt 0 then Y_in = yc else Y_in = dims[1]/2
      if n_elements(size) gt 0 then begin
        ;; We have the size keyword, use Sx and Sy derived from it.
        X_in += [-1,  1, 1, -1]*Sx/2
        Y_in += [-1, -1, 1,  1]*Sy/2
      endif else begin
        ;; Just pick a size, why not half the original dimensions?
        X_in += [-1,  1, 1, -1]*dims[0]/4
        Y_in += [-1, -1, 1,  1]*dims[1]/4
      endelse
      roi_in = OBJ_NEW('IDLgrROI', X_in, Y_in)

      ;; Define the display image in an unusual way to allow for
      ;; multiple dimensions.
      dispim = bytscl(reform(ims[0:dims[0]*dims[1]-1], dims[0:1]))

      print
      print, 'Use the XROI GUI to either modify an initial ROI/FOV or define a new one from scratch.'
      print, 'Select Quit in the File menu. The last ROI is used.'
      print
      
      ;; Fire up the XROI GUI.
      xroi, dispim, regions_in = [roi_in], regions_out = roi, /block $
            , tools = ['Translate-Scale', 'Rectangle'] $
            , title = 'Modify or define ROI'
      roi[-1] -> getproperty, roi_xrange = roi_x
      roi[-1] -> getproperty, roi_yrange = roi_y

      obj_destroy, roi_in
      obj_destroy, roi

      xc = round(mean(roi_x))
      yc = round(mean(roi_y))
      Sx = round(roi_x[1])-round(roi_x[0])
      Sy = round(roi_y[1])-round(roi_y[0])
      corners = lonarr(4)
      corners[0:1] = [xc-Sx/2, yc-Sy/2]          ; Lower
      corners[2:3] = corners[0:1] + [Sx, Sy] - 1 ; Upper

      size = [Sx, Sy]
      
    end
  endcase

  if keyword_set(nocropping) then begin
    print, inam + ' : Returning after establishing the corners.'
    return, 0
  endif
  
  Sx = long(corners[2] - corners[0] + 1)
  Sy = long(corners[3] - corners[1] + 1)

  ;; Check that at least some part of the cropped FOV is within
  ;; bounds. 
  if corners[0] ge dims[0] or $
     corners[1] ge dims[1] or $
     corners[2] lt 0 or $
     corners[3] lt 0 then begin
    print, inam + ' : Cropped FOV is completely outside of original FOV.'
    ;; The best we can do here is to return an array with the right
    ;; size filled with padding.
    if n_elements(pad) eq 0 || size(pad, /tname) eq 'STRING' then begin
      ;; We'll use the default zero for string-valued pad because
      ;; "mean" and "median" have no well defined meaning in this case
      ;; (since we normally apply the operation to the cropped FOV).
      padding = 0 
    endif else begin
      ;; Numerical pad is OK.
      padding = pad
    endelse
    if n_elements(dims) eq 2 then newdims = [Sx, Sy] else newdims = [Sx, Sy, dims[2:*]]
    return, padding + make_array(dimension = newdims, type = size(ims, /type))
  endif
  
  if size(ims, /n_dim) lt 2 then begin
    print, inam + ' : Cannot handle dimensions < 2 or > 8.'
    help, ims
    stop
  endif

  ;; Actual spatial corners in the input array, safeguarded for
  ;; out-of-bounds indices.
  in_corners = lonarr(4)
  in_corners[0:1] = corners[0:1] >0 <(dims-1)   ; Lower
  in_corners[2:3] = corners[2:3] >0 <(dims-1)   ; Upper

  ;; Spatial corners in the output array. Normally [0,Sx-1,0,Sy-1] but
  ;; needs to match the possibly undersized (because of out-of-bounds)
  ;; cropped image.
  out_corners = [0, 0, Sx-1, Sy-1] ; If everything is within bounds
  if corners[0] lt 0       then out_corners[0] = -corners[0] >0 <(Sx-1)
  if corners[1] lt 0       then out_corners[1] = -corners[1] >0 <(Sy-1)
  if corners[2] ge dims[0] then out_corners[2] = (Sx-1 + (dims[0]-corners[2]-1)) >0 <(Sx-1)
  if corners[3] ge dims[1] then out_corners[3] = (Sy-1 + (dims[1]-corners[3]-1)) >0 <(Sy-1)
  
  ;; Do we need padding?
  need_padding = out_corners[2] - out_corners[0] ne Sx or $
                 out_corners[3] - out_corners[1] ne Sy 

  ;; Define padding if not to be calculated later
  case 1 of
    n_elements(pad) eq 0 : begin
      ;; Keyword pad undefined, use default.
      padding = 0    
    end
    size(pad, /tname) eq 'STRING' : begin
      ;; Keyword pad is a string, leave padding undefined if
      ;; calculated later, else use default.
      if pad ne 'median' && pad ne 'mean' then padding = 0 
    end
    else : begin
      ;; Use numerical value of keyword pad
      padding = pad
    end
  endcase
  
  case size(ims, /n_dim) of

    2 : begin
      ;; This is a 2D image, return a 2D image
      newim = make_array(type = size(ims, /type), dimension = [Sx, Sy])

      if need_padding then begin
        ;; Pad the image array before doing the cropping.
        if n_elements(padding) ne 0 then begin
          ;; Use padding defined above
          newim += padding
        endif else begin
          ;; Keyword pad is the name of an operation to be applied to
          ;; the cropped image.
          case pad of
            'median' : newim += median( ims[in_corners[0]:in_corners[2], in_corners[1]:in_corners[3]] )
            'mean'   : newim += mean(   ims[in_corners[0]:in_corners[2], in_corners[1]:in_corners[3]] )
            else     : stop
          endcase
        endelse
      endif

      ;; Do the cropping
      newim[out_corners[0]:out_corners[2], out_corners[1]:out_corners[3]] $
         = ims[in_corners[0]:in_corners[2], in_corners[1]:in_corners[3]]

      return, newim

    end

    3 : begin

      ;; This is a 3D image cube, return a 3D image cube
      newims = make_array(type = size(ims, /type), dimension = [Sx, Sy, dims[2]])

      if need_padding then begin
        ;; Pad the image array before doing the cropping.
        if n_elements(padding) ne 0 then begin
          newims += padding
        endif else begin
          ;; Keyword pad is the name of an operation to be applied to
          ;; each cropped image.
          for i2 = 0, dims[2]-1 do begin
            case pad of
              'median' : newims[*, *, i2] $
                  = median(ims[in_corners[0]:in_corners[2], in_corners[1]:in_corners[3], i2])
              'mean'   : newims[*, *, i2] $
                 = mean(ims[in_corners[0]:in_corners[2], in_corners[1]:in_corners[3], i2])
              else     : stop
            endcase
          endfor                ; i2
        endelse
      endif

      ;; Do the cropping
      newims[out_corners[0]:out_corners[2], out_corners[1]:out_corners[3], *] $
         = ims[in_corners[0]:in_corners[2], in_corners[1]:in_corners[3], *]

      return, newims

    end
    
    4 : begin
      
      ;; This is a 4D image cube, return a 4D image cube
      newims = make_array(type = size(ims, /type), dimension = [Sx, Sy, dims[2:3]])

      if need_padding then begin
        ;; Pad the image array before doing the cropping.
        if n_elements(padding) ne 0 then begin
          newims += padding
        endif else begin
          ;; Keyword pad is the name of an operation to be applied to
          ;; each cropped image.
          for i2 = 0, dims[2]-1 do $
             for i3 = 0, dims[3]-1 do begin
            case pad of
              'median' : newims[*, *, i2, i3] $
                  = median(ims[in_corners[0]:in_corners[2], in_corners[1]:in_corners[3], i2, i3])
              'mean'   : newims[*, *, i2, i3] $
                 = mean(ims[in_corners[0]:in_corners[2], in_corners[1]:in_corners[3], i2, i3])
              else     : stop
            endcase
          endfor                ; i2, i3
        endelse
      endif

      newims[out_corners[0]:out_corners[2], out_corners[1]:out_corners[3], *, *] $
         = ims[in_corners[0]:in_corners[2], in_corners[1]:in_corners[3], *, *]
      
      return, newims

    end
    
    5 : begin
      
      ;; This is a 5D image cube, return a 5D image cube
      newims = make_array(type = size(ims, /type), dimension = [Sx, Sy, dims[2:4]])

      if need_padding then begin
        ;; Pad the image array before doing the cropping.
        if n_elements(padding) ne 0 then begin
          newims += padding
        endif else begin
          ;; Keyword pad is the name of an operation to be applied to
          ;; each cropped image.
          for i2 = 0, dims[2]-1 do $
             for i3 = 0, dims[3]-1 do $
                for i4 = 0, dims[3]-1 do begin
            case pad of
              'median' : newims[*, *, i2, i3, i4] $
                  = median(ims[in_corners[0]:in_corners[2], in_corners[1]:in_corners[3], i2, i3, i4])
              'mean'   : newims[*, *, i2, i3, i4] $
                 = mean(ims[in_corners[0]:in_corners[2], in_corners[1]:in_corners[3], i2, i3, i4])
              else     : stop
            endcase
          endfor                ; i2, i3, i4
        endelse
      endif

      newims[out_corners[0]:out_corners[2], out_corners[1]:out_corners[3], *, *, *] $
         = ims[in_corners[0]:in_corners[2], in_corners[1]:in_corners[3], *, *, *]

      return, newims

    end     

    6 : begin

      ;; This is a 6D image cube, return a 6D image cube
      newims = make_array(type = size(ims, /type), dimension = [Sx, Sy, dims[2:5]])

      if need_padding then begin
        ;; Pad the image array before doing the cropping.
        if n_elements(padding) ne 0 then begin
          newims += padding
        endif else begin
          ;; Keyword pad is the name of an operation to be applied to
          ;; each cropped image.
          for i2 = 0, dims[2]-1 do $
             for i3 = 0, dims[3]-1 do $
                for i4 = 0, dims[3]-1 do $
                   for i5 = 0, dims[3]-1 do begin
            case pad of
              'median' : newims[*, *, i2, i3, i4, i5] $
                  = median(ims[in_corners[0]:in_corners[2], in_corners[1]:in_corners[3], i2, i3, i4, i5])
              'mean'   : newims[*, *, i2, i3, i4, i5] $
                 = mean(ims[in_corners[0]:in_corners[2], in_corners[1]:in_corners[3], i2, i3, i4, i5])
              else     : stop
            endcase
          endfor                ; i2, i3, i4, i5
        endelse
      endif

      newims[out_corners[0]:out_corners[2], out_corners[1]:out_corners[3], *, *, *, *] $
         = ims[in_corners[0]:in_corners[2], in_corners[1]:in_corners[3], *, *, *, *]

      return, newims
      
    end
    
    7 : begin

      ;; This is a 7D image cube, return a 7D image cube
      newims = make_array(type = size(ims, /type), dimension = [Sx, Sy, dims[2:6]])

      if need_padding then begin
        ;; Pad the image array before doing the cropping.
        if n_elements(padding) ne 0 then begin
          newims += padding
        endif else begin
          ;; Keyword pad is the name of an operation to be applied to
          ;; each cropped image.
          for i2 = 0, dims[2]-1 do $
             for i3 = 0, dims[3]-1 do $
                for i4 = 0, dims[3]-1 do $
                   for i5 = 0, dims[3]-1 do $
                      for i6 = 0, dims[3]-1 do begin
            case pad of
              'median' : newims[*, *, i2, i3, i4, i5, i6] $
                  = median(ims[in_corners[0]:in_corners[2], in_corners[1]:in_corners[3], i2, i3, i4, i5, i6])
              'mean'   : newims[*, *, i2, i3, i4, i5, i6] $
                 = mean(ims[in_corners[0]:in_corners[2], in_corners[1]:in_corners[3], i2, i3, i4, i5, i6])
              else     : stop
            endcase
          endfor                ; i2, i3, i4, i5, i6
        endelse
      endif

      newims[out_corners[0]:out_corners[2], out_corners[1]:out_corners[3], *, *, *, *, *] $
         = ims[in_corners[0]:in_corners[2], in_corners[1]:in_corners[3], *, *, *, *, *]

      return, newims
      
    end
    
    8 : begin

      ;; This is a 8D image cube, return a 8D image cube
      newims = padding + make_array(type = size(ims, /type), dimension = [Sx, Sy, dims[2:7]])

      if need_padding then begin
        ;; Pad the image array before doing the cropping.
        if n_elements(padding) ne 0 then begin
          newims += padding
        endif else begin
          ;; Keyword pad is the name of an operation to be applied to
          ;; each cropped image.
          for i2 = 0, dims[2]-1 do $
             for i3 = 0, dims[3]-1 do $
                for i4 = 0, dims[3]-1 do $
                   for i5 = 0, dims[3]-1 do $
                      for i6 = 0, dims[3]-1 do $
                         for i7 = 0, dims[3]-1 do begin
            case pad of
              'median' : newims[*, *, i2, i3, i4, i5, i6, i7] $
                  = median(ims[in_corners[0]:in_corners[2], in_corners[1]:in_corners[3], i2, i3, i4, i5, i6, i7])
              'mean'   : newims[*, *, i2, i3, i4, i5, i6, i7] $
                 = mean(ims[in_corners[0]:in_corners[2], in_corners[1]:in_corners[3], i2, i3, i4, i5, i6, i7])
              else     : stop
            endcase
          endfor                ; i2, i3, i4, i5, i6, i7
        endelse
      endif

      newims[out_corners[0]:out_corners[2], out_corners[1]:out_corners[3], *, *, *, *, *, *] $
         = ims[in_corners[0]:in_corners[2], in_corners[1]:in_corners[3], *, *, *, *, *, *]

      return, newims
      
    end
    
    else : begin
      print, inam + ' : Cannot handle dimensions > 8.'
      help, ims
      stop
    end

  endcase

end

oldimage = findgen(100, 100, 5, 3, 4)

; Test out-of-bounds
;newimage = red_crop(oldimage, corners = corners, xc = 99, yc = 2, size = 12, pad = 'median')

; Test oversized image
;newimage = red_crop(oldimage, corners = corners, /center, size = 120, pad = 'median')


; Test ROI selection
newimage = red_crop(oldimage, corners = corners, size = size, /nocrop)

end
