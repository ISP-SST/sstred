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
;-
pro red_missing, image $
                 , image_out = image_out $
                 , indx_data = indx_data $
                 , indx_missing = indx_missing $
                 , inplace = inplace $
                 , missing_value = missing_value $
                 , missing_type_wanted = missing_type_wanted $
                 , missing_type_used = missing_type_used $
                 , verbose = verbose

  inam = red_subprogram(/low, calling = inam1)

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
  ;; but not necessarily all. The current logic should be good for the
  ;; first case and perhaps for the second. We need to make sure it
  ;; works for the second and third cases.


  ;; Heuristics for identifying the missing-data pixels. They would be
  ;; connected to the outermost rows and columns and they would all
  ;; have the same value, usually the median of the data pixels. We
  ;; can't simply assume all pixels with this value are missing data
  ;; as there *could* be interior pixels with this exact value.

  ;; What kind of padding do we have now?
  
  currently_constant = image[ 0,  0] eq image[ 0, -1] and $
                       image[ 0,  0] eq image[-1,  0] and $
                       image[ 0,  0] eq image[-1, -1]

  currently_nans = ~( finite(image[ 0,  0]) or $
                      finite(image[ 0, -1]) or $
                      finite(image[-1,  0]) or $
                      finite(image[-1, -1]) $
                    )


  case 1 of
    
    n_elements(missing_value) gt 0       : missing_type_used = 'value'
  
    n_elements(missing_type_wanted) gt 0 : missing_type_used = missing_type_wanted

    else : begin
      
      ;; Not specified, change to the opposite of what we have

      case 1 of
        
        currently_nans     : missing_type_used = 'median'
        
        currently_constant : missing_type_used = 'nan'
        
        else : begin
          if keyword_set(verbose) then print, inam + ' : Could not identify missing data pixels'
          if keyword_set(verbose) then print, 'Corner pixel intensities: ' $
                                              , image[ 0,  0] $
                                              , image[ 0, -1] $
                                              , image[-1,  0] $
                                              , image[-1, -1]
          missing_type_used = 'none'
        end
        
      endcase
    end
    
  endcase
    
  case strlowcase(missing_type_used) of

    'nan' : begin               ; Assume all corner pixels have the same value.
      
      if currently_nans then begin
        if keyword_set(verbose) then print, inam+' : Padding seems to be NaN already.'
        indx_data = where(finite(image), Ndata, complement = indx_missing, Ncomplement = Nmissing)
        image_out = image       ; No change
        return
      endif
      
      ;; Find pixels with the same value as the corners
      mask = image eq image[ 0,  0] 
      indx_missing = where(mask, Nmissing, complement = indx_data, Ncomplement = Ndata)

      if Nmissing eq 1 then Nmissing = 0 ; image[0,0] is only equal to itself
      
      ;; We want NaNs
      missing_value = !Values.F_NaN
      
    end

    'median' : begin

      if currently_constant then begin
        if keyword_set(verbose) then print, inam+" : Padding seems to be constant. Set to median."
        mask = image eq image[ 0,  0] 
      endif else begin
        if keyword_set(verbose) then print, inam+" : Padding seems to be NaN. Set to median."
        ;; Find non-finite pixels
        mask = ~finite(image)
      endelse
      
      ;; We want the median
      ;;indx_data = where(~mask, Ndata)
      indx_missing = where(mask, Nmissing, complement = indx_data, Ncomplement = Ndata)
      if Ndata eq 0 then stop
      missing_value = median(image(indx_data))
      
    end
 
    'value' : begin

      if currently_constant then begin
        if keyword_set(verbose) then print, inam+" : Padding seems to be constant. Set to the given value."
        mask = image eq image[ 0,  0] 
      endif else begin
        if keyword_set(verbose) then print, inam+" : Padding seems to be NaN. Set to the given value."
        ;; Find non-finite pixels
        mask = ~finite(image)
      endelse
      
      indx_missing = where(mask, Nmissing, complement = indx_data, Ncomplement = Ndata)
      if Ndata eq 0 then stop
       
    end

    'none' : begin

      Ndata = n_elements(image)
      indx_data = indgen(Ndata)
      Nmissing = 0
      
    end
    
    else :  stop
    
  endcase

  if arg_present(image_out) then image_out = image
    
  if Nmissing eq 0 then return

  ;; Use label_region to find the pixels that are connected to
  ;; the corners
  mask = red_centerpic(mask, xSize = Nx+8, ySize = Ny+8, z = 1) ; Add some rows and columns
  label = label_region(mask)                                    ; Label regions
  label = label[4:-5, 4:-5]                                     ; Remove extra rows and columns
  mask = label eq label[ 0, 0] $                                ; Connected to any of the corners
         or label eq label[ 0,-1] $
         or label eq label[-1, 0] $
         or label eq label[-1,-1] 
  ;; Change the value of those pixels to NaN
  indx_missing = where(mask, Nmissing)
  
  if keyword_set(InPlace) then if Nmissing gt 0 then $
     image[indx_missing] = missing_value

  if arg_present(image_out) then if Nmissing gt 0 then $
     image_out[indx_missing] = missing_value

end

dir = '/scratch/mats/2016.09.19/CRISP-aftersummer/cubes_nb/'
filename = dir+'nb_6302_2016-09-19T09:30:20_scans=2-8_stokes_corrected_im.fits'

red_fitscube_getframe, filename, image_orig, iframe = 0

im = red_centerpic(image_orig,sz=500)
red_missing, im $
             , image_out = image_med $
             , indx_data = indx_data $
             , indx_missing = indx_missing $
             , missing_value = 0. $
             , missing_type_used = missing_type_used

stop



;; Four corners with missing data, but disconnected from each other.
;; This works!
;image = red_centerpic(image, sz = 1000)


;; Three corners, connected. This does not work!
image_orig = image[0:1000, 0:1000]

image = image_orig

red_show, image, w = 1

red_missing, image $
             , image_out = image_med $
             , indx_data = indx_data $
             , indx_missing = indx_missing $
             , missing_type_wanted = 'median' $
             , missing_type_used = missing_type_used

red_show, image_med, w = 2

red_missing, image_med $
             , image_out = image_nan $
             , indx_data = indx_data $
             , indx_missing = indx_missing $
             , missing_type_wanted = 'nan' $
             , missing_type_used = missing_type_used

red_show, image_nan, w = 3


stop

red_missing, image $
             , missing_type_wanted = 'median' $
             , missing_type_used = missing_type_used


end
