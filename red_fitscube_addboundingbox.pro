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
;    anchor : in, out, optional, type=string
; 
;      See red_fitsaddkeyword.pro
; 
;    checksum : in, optional, type=boolean
;   
;      Check the CHECKSUM keyword before adding new keywords and
;      update it afterwards.
; 
;    hdr : in, out, optional, type=strarr
; 
;      A header to be updated. If present, don't read the header from
;      the file and don't write it either. Also implies checksum=0.
; 
; :History:
; 
;    2026-05-12 : MGL. First version.
; 
;-
pro red_fitscube_addboundingbox, filename, anchor = anchor, checksum = checksum, hdr = hdr

  hdr_present = arg_present(hdr) 

  if ~hdr_present then begin
    ;; We didn't provide a header, so read it from the file
    hdr = headfits(filename)
  endif else begin
    ;; We provided a header so we will not write it to the file. Leave
    ;; doing a checksum test to the calling program.
    checksum = 0
  endelse
  
  if keyword_set(checksum) then begin
    checksum_status = fits_test_checksum(hdr, /trust_datasum)
    if checksum_status eq -1 then begin
      red_message, ['CHECKSUM keyword does not have the' $
                    , 'correct value, indicating possible data' $
                    , 'corruption.']
      hgrep,hdr,'CHECK'
      red_message, 'Will not add BoundingBox keywords.'
      return
    endif
  endif
  
  ;; Read WCS cordinates from the file
  red_fitscube_getwcs, filename, coordinates = coordinates
  tags = tag_names(coordinates)

  ;; Get the WCS coordinates in order
  ctype = fxpar(hdr, 'CTYPE*')

  ;; CBBMINnn/CBBMAXnn - Coordinate Bounding Box Min/Max for koordinat nn

  if n_elements(anchor) eq 0 then anchor = 'FILENAME'
  
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

      'WAVE' : begin
        if strmatch(file_basename(filename), 'wb_*') then begin
          ;; USE WAVEMIN/MAX för WB data
          mn = fxpar(hdr, 'WAVEMIN')
          mx = fxpar(hdr, 'WAVEMAX')
        endif else begin
          ;; Use min and max of WCS coordinate
          indx = where(tags eq 'WAVE', Nwhere)
          if Nwhere eq 0 then stop
          mn = min(coordinates.(indx[0]), max = mx)
        endelse
      end
      
      else : begin
        indx = where(tags eq coord_name, Nwhere)
        if Nwhere eq 0 then stop
        mn = min(coordinates.(indx[0]), max = mx)
      end
      
    endcase

      
    mncomment = 'Coordinate BB min for '+coord_name+' coordinate'
    mxcomment = 'Coordinate BB max for '+coord_name+' coordinate'

    red_fitsaddkeyword, hdr, 'CBBMIN'+strtrim(icoord+1, 2), mn, mncomment, anchor = anchor
    red_fitsaddkeyword, hdr, 'CBBMAX'+strtrim(icoord+1, 2), mx, mxcomment, anchor = anchor

    if icoord eq 0 then begin
      print
      print, 'Adding BoundingBox keywords:'
      print
    endif
    print, 'CBBMIN'+strtrim(icoord+1, 2), mn, '  ' + mncomment
    print, 'CBBMAX'+strtrim(icoord+1, 2), mx, '  ' + mncomment
    
  endfor                        ; icoord
  print
  
  if ~hdr_present then begin  ;; Write the new header to the file
    red_fitscube_newheader, filename, hdr
  endif
  
  if keyword_set(checksum) then begin
    ;; If there is a main header checksum, then re-calculate it
    ;; The BB keywords are meant to be added before checksums are added
    ;; to the file.
    if checksum_status eq 1 then begin
      ;; CHECKSUM keyword is present with correct values. Update the
      ;; main HDU checksum since we added keywords.
      hgrep,hdr,'CHECK'
      red_fitscube_checksums, filename $
                              , hdus = 0 $
                              , /trust_datasum 
      hdr = headfits(filename)
      hgrep,hdr,'CHECK'
    endif
  endif
  
end

filename = '/scratch/mats/NEW/2024-04-21/CRISP-merged/cubes_nb/nb_6173_2024-04-21T10:46:39_10:46:39=24,25_stokes_corrected_im.fits'

red_fitscube_addboundingbox, filename

end
