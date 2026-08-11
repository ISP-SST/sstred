; docformat = 'rst'

;+
; Calculate hpln and hplt coordinates, optionally for the four corners
; of an image. 
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
;    time : in, type=float
;
;      The observation time of the image in seconds since midnight.
;
;    pointing : in, type="fltarr(2,Ntime)"
;
;      The pointing coordinates as a function of time, e.g., from PIG.
;
;    pointing_time : in, type="fltarr(Ntime)"
;
;      The time coordinates corresponding to the pointing array in
;      seconds since midnight.
;
;    hpln : out, type="dblarr(2,2)"
;
;       The HPLN coordinates of the FOV corners.
;
;    hplt : out, type="dblarr(2,2)"
;
;       The HPLT coordinates of the FOV corners.
;
; 
; 
; 
; 
; :Keywords:
; 
;    Nx : in, optional, type=integer
;
;      The image X size in pixels.
;
;    Ny : in, optional, type=integer
;
;      The image Y size in pixels.
;
;    image_scale : in, optional, type=float
;
;      The image scale in arcsec/pixel.
;   
; 
; 
; :History:
; 
;   2017-08-17 : MGL. First version.
; 
;   2017-10-23 : MGL. Change parameters Nx, Ny, and image_scale into
;                (optional) keywords.
; 
;-
pro red_wcs_hpl_coords, time, pointing, pointing_time, hpln, hplt $
                        , Nx = Nx, Ny = Ny, image_scale = image_scale

  ;; Coordinates for the center of the FOV:

  hpln = interpol(pointing[0, *], pointing_time, time)
  hplt = interpol(pointing[1, *], pointing_time, time)

  if n_elements(Nx) ne 0 and n_elements(Ny) ne 0 and n_elements(image_scale) ne 0 then begin
    
    ;; Now tabulate the corner coordinates, assuming the FOV is aligned
    ;; to solar coordinates. [The distance between the center of the FOV
    ;; and the centers of the corner pixels is pixelsize*(Nx-1)/2 and
    ;; pixelsize*(Ny-1), resp.]

    hpln = hpln + [[-1, 1],[-1, 1]] * double(image_scale) * (Nx-1)/2.d
    hplt = hplt + [[-1,-1],[ 1, 1]] * double(image_scale) * (Ny-1)/2.d

  endif
  
end
