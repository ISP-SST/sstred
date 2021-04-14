; docformat = 'rst'

;+
; Find pixels with missing data and optionally set them to NaN,
; median(), or a given value.
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
; :Params:
; 
;   image : in, type=string
; 
;     The image.
; 
; 
; :Keywords:
; 
;   image_out : out, optional, type=array
;
;      The image, with the missing data set to the value indicated by
;      missing_type_used. 
; 
;   indx_data : out, optional, type=array
;
;      Indices of pixels identified as data.
;
;   indx_missing : out, optional, type=array
; 
;      Indices of pixels identified as missing data.
; 
;   inplace : in, optional, type=boolean
;
;      Change the padding inplace in the input image.
;
;   missing_type_used : out, optional, type=string
;
;     The type of value missing-data pixels were set to. One of 'nan',
;     'median' or 'value' (the latter if missing_value was used).
;
;   missing_type_wanted : in, optional, type=string, default='opposite'
;
;     What type of value to set missing-data pixels to. One of 'nan'
;     and 'median'. By default, set it to the opposite of what it
;     appears to be in the input image.
;
;   missing_value : in, optional, type=number
;
;     The value to set missing-data pixels to. If given,
;     missing_type_wanted is ignored.
;
; :History:
; 
;     2020-07-16 : MGL. First version.
; 
;     2021-04-08 : MGL. New keywords inplace and missing_value.
; 
;     2021-04-14 : MGL. Rewritten. New keywords Ndata and Nmissing.
; 
;-
pro red_missing, image $
                 , image_out = image_out $
                 , indx_data = indx_data $
                 , indx_missing = indx_missing $
                 , inplace = inplace $
                 , missing_value = missing_value $
                 , missing_type_wanted = missing_type_wanted $
                 , missing_type_used = missing_type_used $
                 , ndata = ndata $
                 , nmissing = nmissing $
                 , verbose = verbose

  inam = red_subprogram(/low, calling = inam1)

  undefine, Ndata
  undefine, Nmissing
  undefine, indx_data
  undefine, indx_missing
  undefine, image_out
  undefine, missing_type_used

  dims = size(image, /dim)
  Nx = dims[0]
  Ny = dims[1]
  
  ;; As long as we are dealing with un-cropped images that have been
  ;; derotated and shifted, it is likely that all missing data pixels
  ;; are connected. If cropped slightly, any of the corners or all of
  ;; them) may be disconnected from the others. If cropped a lot,
  ;; there could be missing data in only one or a few of the corners
  ;; but not necessarily all. 

  
  ;; First, find out what padding we already have, if any.

  ;; Is any of the corner pixels a NaN?
  currently_nans = ~( finite(image[ 0,  0]) and $
                      finite(image[ 0, -1]) and $
                      finite(image[-1,  0]) and $
                      finite(image[-1, -1]) $
                    )
  ;; Is any of the corner 2x2 sub-matrices constant?
  currently_constant = stddev(image[ 0: 1, 0: 1]) eq 0 or $
                         stddev(image[ 0: 1,-2:-1]) eq 0 or $
                         stddev(image[-2:-1, 0: 1]) eq 0 or $
                         stddev(image[-2:-1,-2:-1]) eq 0 


  ;; Some simple actions depending on detected padding
  
  case 1 of

    ~currently_nans and ~currently_constant : begin
      ;; If there is no padding, just set some keywords and return.
      Nmissing = 0
      Ndata = n_elements(image)
      if arg_present(image_out) then image_out = image 
      if arg_present(indx_data) then indx_data = lindgen(n_elements(image))
      return
    end

    currently_nans and currently_constant : begin
      ;; Could there be both kinds of padding? Warn about it because we
      ;; have not implemented support for this case yet.
      print, inam + ' : Both NaN and finite padding detected. Alert developers!'
      stop
    end

    currently_nans : begin
      if keyword_set(verbose) then print, inam + ' : Detected NaN padding'
    end

    currently_constant : begin
      if keyword_set(verbose) then print, inam + ' : Detected constant padding'
    end
    
  endcase

  
  ;; At this point we have ruled out cases where there is no padding
  ;; and cases where we don't know if the padding is NaNs or constant
  ;; values. From now on it should be one of those two.

  
  ;; Now, find out which pixels are padding and which are data.

  case 1 of
    
    currently_nans : begin

      ;; We define the padding as _all_ NaN pixels, regardless of
      ;; whether they are connected to a corner.

      indx_data = where(finite(image), Ndata, complement = indx_missing, Ncomplement = Nmissing)

    end

    currently_constant : begin

      ;; There is padding with a finite value. We don't know if
      ;; it's the median or something else. 
      
      ;; Find areas of constant values connected to the corners. 
  
      mask1 = image eq image[ 0, 0]                                   ; Same intensity as corner 1
      label = red_centerpic(mask1, xSize = Nx+8, ySize = Ny+8, z = 1) ; Add some rows and columns
      label = label_region(label)                                     ; Label regions
      label = label[4:-5, 4:-5]                                       ; Remove extra rows and columns
      mask1 = label eq label[0, 0]                                    ; Connected to corner 1

      mask2 = image eq image[ 0,-1]                                   ; Same intensity as corner 2
      label = red_centerpic(mask2, xSize = Nx+8, ySize = Ny+8, z = 1) ; Add some rows and columns
      label = label_region(label)                                     ; Label regions
      label = label[4:-5, 4:-5]                                       ; Remove extra rows and columns
      mask2 = label eq label[ 0,-1]                                   ; Connected to corner 2

      mask3 = image eq image[-1, 0]                                   ; Same intensity as corner 3
      label = red_centerpic(mask3, xSize = Nx+8, ySize = Ny+8, z = 1) ; Add some rows and columns
      label = label_region(label)                                     ; Label regions
      label = label[4:-5, 4:-5]                                       ; Remove extra rows and columns
      mask3 = label eq label[-1, 0]                                   ; Connected to corner 3

      mask4 = image eq image[-1,-1]                                   ; Same intensity as corner 4
      label = red_centerpic(mask4, xSize = Nx+8, ySize = Ny+8, z = 1) ; Add some rows and columns
      label = label_region(label)                                     ; Label regions
      label = label[4:-5, 4:-5]                                       ; Remove extra rows and columns
      mask4 = label eq label[-1,-1]                                   ; Connected to corner 4

      ;; Make a mask with all the padded areas. The individual corner
      ;; masks are padding only if they are larger than 1.
      mask = bytarr(Nx, Ny)   
      if round(total(mask1)) gt 1 then mask or= mask1
      if round(total(mask2)) gt 1 then mask or= mask2
      if round(total(mask3)) gt 1 then mask or= mask3
      if round(total(mask4)) gt 1 then mask or= mask4

      indx_missing = where(mask, Nmissing, complement = indx_data, Ncomplement = Ndata)

    end
    
  endcase
  

  
  ;; What kind of padding should we use in the output?
  
  case 1 of
    
    n_elements(missing_value) gt 0 : begin
      ;; A value was provided with the missing_value keyword
      missing_type_used = 'value'
    end
    
    n_elements(missing_type_wanted) gt 0 : begin
      ;; A type was provided with the missing_type_wanted keyword
      missing_type_used = missing_type_wanted
    end
    
    else : begin
      ;; No type or value specified, change to the opposite of what we
      ;; have
      case 1 of
        currently_nans     : missing_type_used = 'median'
        currently_constant : missing_type_used = 'nan'
      endcase
    end
    
  endcase

  ;; If we don't want to change the image or set image_out then
  ;; we can return now.
  if ~arg_present(image_out) and ~keyword_set(inplace) then return

  
  ;; Now actions depending on combinations of current padding and
  ;; padding to use.

  case 1 of
    
    currently_nans : begin

      case strlowcase(missing_type_used) of

        'nan' : begin
          ;; Padding is what we want, no need to change the input image
          if arg_present(image_out) then image_out = image 
          return
        end

        'median' : begin
          missing_value = median(image[indx_data])
        end

        'value' : begin
          ;; missing_value already set 
        end 

        else : stop             ; Should not happen!
        
      endcase                   ; missing_type_used

    end                         ; currently_nans
    
    currently_constant : begin

      case strlowcase(missing_type_used) of
        
        'nan' : begin
          missing_value = !Values.F_NaN
        end

        'median' : begin
          missing_value = median(image[indx_data])
        end

        'value' : begin
          ;; missing_value already set 
        end

        else : stop             ; Should not happen!

      endcase                   ; missing_type_used

    end                         ; currently_constant
    
  endcase                       ; currently


  ;; Finally, pad the input image and/or set image_out to the
  ;; correctly padded image.

  if keyword_set(inplace) then begin
    ;; Modify the input image
    image[indx_missing] = missing_value
    ;; Did we want image_out as well?
    if arg_present(image_out) then image_out = image
    ;; Done!
    return
  endif

  if arg_present(image_out) then begin
    ;; Make an image_out
    image_out = image
    image_out[indx_missing] = missing_value
    return
  endif

  ;; How did we get here?
  stop

end

dir = '/scratch/mats/2016.09.19/CRISP-aftersummer/cubes_nb/'
filename = dir+'nb_6302_2016-09-19T09:30:20_scans=2-8_stokes_corrected_im.fits'

red_fitscube_getframe, filename, image_orig, iframe = 0

if 1 then begin

  ;; Test inplace and missing_value. Works.

  image_orig = red_centerpic(image_orig, sz = 1000)
  image = image_orig
  
  red_show, image, w = 1
  
  red_missing, image $
               , ndata = ndata $
               , nmissing = nmissing $
               , image_out = image_med $
               , indx_data = indx_data $
               , indx_missing = indx_missing $
                                ;, missing_type_wanted = 'median' $
               , missing_type_used = missing_type_used
  red_show, image_med, w = 2
  
  red_missing, image $
               , ndata = ndata $
               , nmissing = nmissing $
               , /inplace $
               , indx_data = indx_data $
               , indx_missing = indx_missing $
               , missing_value = max(image) $
               , missing_type_used = missing_type_used

  red_show, image, w = 3

endif


if 0 then begin
  ;; No padding. This works!
  im = red_centerpic(image_orig,sz=500)
  red_missing, im $
             , ndata = ndata $
             , nmissing = nmissing $
             , image_out = image_med $
             , indx_data = indx_data $
             , indx_missing = indx_missing $
             , missing_value = 0. $
             , missing_type_used = missing_type_used
  print, Ndata, Nmissing
endif


if 0 then begin

  ;; Four corners with missing data, but disconnected from each other.
  ;; This works!

  image_orig = red_centerpic(image_orig, sz = 1000)
  image = image_orig
  
  red_show, image, w = 1
  
  red_missing, image $
               , ndata = ndata $
               , nmissing = nmissing $
               , image_out = image_med $
               , indx_data = indx_data $
               , indx_missing = indx_missing $
               , missing_type_wanted = 'median' $
               , missing_type_used = missing_type_used
  
  print, Ndata, Nmissing
  
  red_show, image_med, w = 2
  
  red_missing, image_med $
               , ndata = ndata $
               , nmissing = nmissing $
               , image_out = image_nan $
               , indx_data = indx_data $
               , indx_missing = indx_missing $
               , missing_type_wanted = 'nan' $
               , missing_type_used = missing_type_used
  
  print, Ndata, Nmissing
  
  red_show, image_nan, w = 3
endif

if 0 then begin
  ;; Three corners, connected. This works!
  image_orig = image_orig[0:1000, 0:1000]

  image = image_orig
  
  red_show, image, w = 1
  
  red_missing, image $
               , ndata = ndata $
               , nmissing = nmissing $
               , image_out = image_med $
               , indx_data = indx_data $
               , indx_missing = indx_missing $
               , missing_type_wanted = 'median' $
               , missing_type_used = missing_type_used
  
  print, Ndata, Nmissing
  
  red_show, image_med, w = 2
  
  red_missing, image_med $
               , ndata = ndata $
               , nmissing = nmissing $
               , image_out = image_nan $
               , indx_data = indx_data $
               , indx_missing = indx_missing $
               , missing_type_wanted = 'nan' $
               , missing_type_used = missing_type_used
  
  print, Ndata, Nmissing
  
  red_show, image_nan, w = 3
endif

if 0 then begin

  red_missing, image $
               , missing_type_wanted = 'median' $
               , missing_type_used = missing_type_used
endif


end
