; docformat = 'rst'

;+
; Calculate hpln and hplt coordinates for the four corners of an image.
;
; Todo: implement rotation of the FOV.
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
;    Nx : in, type=integer
;
;      The image X size in pixels.
;
;    Ny : in, type=integer
;
;      The image Y size in pixels.
;
;    image_scale : in, type=
;
;
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
;   
;   
;   
; 
; 
; :History:
; 
;   2017-08-17 : MGL. First version.
; 
; 
; 
; 
;-
pro red_wcs_hpl_coords, time, pointing, pointing_time, Nx, Ny, image_scale, hpln, hplt

  ;; Coordinates for the center of the FOV:

  hpln = interpol(pointing[0, *], pointing_time, time)
  hplt = interpol(pointing[1, *], pointing_time, time)
 
  ;; Now tabulate the corner coordinates, assuming the FOV is aligned
  ;; to solar coordinates. [The distance between the center of the FOV
  ;; and the centers of the corner pixels is pixelsize*(Nx-1)/2 and
  ;; pixelsize*(Ny-1), resp.]

  hpln = hpln + [[-1, 1],[-1, 1]] * double(image_scale) * (Nx-1)/2.d
  hplt = hplt + [[-1,-1],[ 1, 1]] * double(image_scale) * (Ny-1)/2.d
  
end
