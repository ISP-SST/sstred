; docformat = 'rst'

;+
; Add BoundingBox keywords.
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
;    filename : in, type=string
; 
;      The path of the fitscube file. 
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
; 
; 
;-
pro red_fitscube_addboundingbox, filename

  ;; Read WCS cordinates from the file
  red_fitscube_getwcs, filename, coordinates = coordinates
  tags = tag_names(coordinates)

  ;; Get the WCS coordinates in order
  hdr = headfits(filename)
  ctype = fxpar(hdr, 'CTYPE*')

  ;; CBBMINnn/CBBMAXnn - Coordinate Bounding Box Min/Max for koordinat nn

  anchor = 'FILENAME'
  
  for icoord = 0, n_elements(ctype)-1 do begin

    coord_name = strtrim((strsplit(ctype[icoord], '-', /extract))[0], 2)
    
    case coord_name of

      'STOKES' : begin
        ;; This is not tabulated but we can get min and max from the
        ;; appropriate NAXIS keyword
        mn = 1                                        ; min always Stokes I = 1
        mx = fix(fxpar(hdr, 'NAXIS'+strtrim(icoord+1, 2))) ; max 1 or 4 (Stokes I or V)
      end

      'UTC' : begin
        indx = where(tags eq 'TIME', Nwhere)
        if Nwhere eq 0 then stop
        mn = min(coordinates.(indx[0]), max = mx)
      end

      else : begin
        indx = where(tags eq coord_name, Nwhere)
        if Nwhere eq 0 then stop
        mn = min(coordinates.(indx[0]), max = mx)
      end
      
    endcase

    mncomment = 'Coordinate BB min for '+coord_name+' coordinate'
    mxcomment = 'Coordinate BB max for '+coord_name+' coordinate'

    if 1 then begin
      red_fitsaddkeyword, hdr, 'CBBMIN'+strtrim(icoord+1, 2), mn, mncomment, anchor = anchor
      red_fitsaddkeyword, hdr, 'CBBMAX'+strtrim(icoord+1, 2), mx, mxcomment, anchor = anchor
    endif else begin
      print
      print, 'CBBMIN'+strtrim(icoord+1, 2), mn, '  '+mncomment
      print, 'CBBMAX'+strtrim(icoord+1, 2), mx, '  '+mncomment
    endelse
  endfor                        ; icoord

  stop
  red_fitscube_newheader, filename, hdr
  
  
end

filename = '/scratch/mats/NEW/2024-04-21/CRISP-merged/cubes_nb/nb_6173_2024-04-21T10:46:39_10:46:39=24,25_stokes_corrected_im.fits'

red_fitscube_addboundingbox, filename

end
