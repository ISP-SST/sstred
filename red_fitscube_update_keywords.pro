; docformat = 'rst'

;+
; Update a fitscube for GaVO DaCHS if needed.
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
;   filename : in, type=string
; 
;     The file to update.
; 
; 
; :Keywords:
; 
;   do_update : in, optional, type=boolean
;   
;     Do the update even if the GaVO DaCHS keywords already seem to be
;     there.
; 
;   out_dir : in, optional, type=string 
; 
;     Put the updated file here instead of in the original file's
;     directory (default if /remove_original) or the local directory
;     (default otherwise).
; 
;   overwrite : in, optional, type=boolean
;
;     Allow overwriting of existing files.
; 
;   remove_original : in, optional, type=boolean 
; 
;     Remove the original file.
; 
; 
; 
; :History:
; 
;   2026-07-01 : MGL. First version.
; 
;-
pro red_fitscube_update_keywords, filename $
                                  , do_update = do_update $
                                  , out_dir = out_dir $
                                  , remove_original = remove_original

  hdr = headfits(filename)

  is_wb = strmid(file_basename(filename), 0, 2) ne 'wb'
  
  ;; Check if keywords need to be updated, exit if not.
  keys = ['CROTA' $
          , 'CBBMIN'+strtrim(indgen(5)+1, 2) $
          , 'CBBMAX'+strtrim(indgen(5)+1, 2) $
          , 'OBS_MODE' $
          , 'POINT_ID' $
         ]
  if is_wb then begin
    keys = [keys $
            , 'RESOLVPW' $
            , 'CSPMIN3' $
            , 'CSPMAX3' $
           ]
  endif
  Nkeys = n_elements(keys)
  for ikey = 0, Nkeys-1 do begin
    val = red_fitsgetkeyword(hdr, keys[ikey], count = cnt)
    if cnt eq 0 then begin
      do_update = !true
      print, 'Missing: ', keys[ikey]
    endif
  endfor                        ; ikey

  if ~keyword_set(do_update) then begin
    print, 'Keywords are fine: ' + file_basename(filename)
    return
  endif

  
  ;; Check that CHECKSUM is ok.
  checksum_status = fits_test_checksum(hdr, /trust_datasum)
  if checksum_status eq 1 then begin
    print
    red_message, 'Source file checksum present and correct'
    hgrep,hdr,'CHECK'
    print
  endif else begin
    red_message, 'Source file checksum incorrect or not present. Skip this file: ' $
                 + filename
    return
  endelse

  instrument = strtrim(fxpar(hdr, 'INSTRUME'), 2)

  
  old_dir = file_dirname(filename)

  if n_elements(out_dir) ne 0 then begin
    new_dir = out_dir + '/'
  endif else begin
    if keyword_set(remove_original) then begin
      new_dir = old_dir + '/'
    endif else begin
      new_dir = './'
    endelse
  endelse
  tmpname = new_dir + 'tmp' + cgTimeStamp() + '.fits'

  ;; Make a copy of the file to do the work on
  file_copy, filename, tmpname, overwrite = overwrite

  ;; Infer CROTA from the file and give as parameter to
  ;; red_fitscube_addgavokeywords
  undefine, crota               ; Default: do not write the keyword
  
  procs = fxpar(hdr,'PRPROC*')
  indx = where(strmatch(procs, '*::make_scan_cube'), Nwhere)
  if Nwhere gt 0 then begin

    ;; Scan cubes are always rotated to Solar N up
    crota = 0.0

  endif else begin

    ;; For WB and NB cubes we need to find out if make_wb_cube was
    ;; called with the /subtract_meanang option. With this option,
    ;; the average rotation angle was subtracted in order to
    ;; minimize no-data real estate in the final cube.)

    ;; Are we dealing with a WB cube?
    indx = where(procs eq 'red::make_wb_cube', Nwhere)

    if Nwhere eq 0 then begin

      ;; If this is a NB file, we need to find the header of the WB
      ;; file used to make the NB file.
      
      ;;  help, sourcefile
      indx = where(strmatch(procs, '*::make_nb_cube'), Nwhere)
      prpara = json_parse(fxpar(hdr,'PRPARA'+strtrim(indx[0]+1,2)))
      wbfile = file_search(file_dirname(filename) + '/' + file_basename(prpara['WCFILE'], '_im.fits')+'*im.fits')
      if file_test(wbfile) then begin
        wbhdr = headfits(wbfile)
        procs = fxpar(wbhdr,'PRPROC*')
        indx = where(procs eq 'red::make_wb_cube', Nwhere)
      endif
      
    endif else begin

      ;; If this is a WB file, we use its own header
      
      wbhdr = hdr
      
    endelse
    
    if Nwhere gt 0 then begin
      
      ;; Now look in the PRPARA keyword for the /subtract_meanang
      ;; keyword.

      prpara = json_parse(fxpar(wbhdr,'PRPARA'+strtrim(indx[0]+1,2)))
      keys = (prpara.keys()).toarray()
      
      if total(keys eq 'SUBTRACT_MEANANG') eq 0 then begin
        
        ;; make_wb_cube was called without subtract_meanang keyword,
        ;; so the cube should be rotated to Solar N up.
        crota = 0.0
        
      endif else begin
        
        ;; Was subtract_meanang set to true or false?
        
        if prpara['SUBTRACT_MEANANG'] then begin
          
          ;; TRUE: Here we need to figure out what the subtracted
          ;; mean angle was. Calculate the angles the same way as in
          ;; make_wb_cube.
          
          rotation = prpara['ROTATION']
          red_fitscube_getwcs, destfile, coordinates=coordinates
          time = coordinates.time
          time = mean(mean(mean(time,dim=1),dim=1),dim=1)
          date = (strsplit(fxpar(hdr, 'DATE-BEG'), 'T', /extract))[0]
          ang = red_lp_angles(time, date[0], /from_log, offset_angle = rotation)
          crota = median(ang) * 180./!pi
          
        endif else begin
          
          ;; FALSE: the cube should be rotated to Solar N up.
          crota = 0.0
          
        endelse
        
      endelse

    endif
    
  endelse

  
  ;; Update the header with fitscube_addgavo
  red_fitscube_addgavokeywords, tmpname $
                                , crota = crota $
                                , hdr = hdr
  
  ;; Add stepinfo
  red_headerinfo_addstep, hdr
  datestamp = red_timestamp(/utc,/iso)
  fxaddpar, hdr, 'DATE', datestamp

  
  ;; Change filename (and FILENAME), replacing the creation date part.
  splt = strsplit(file_basename(filename), 'export', /extract)
  pos = strpos(file_basename(filename), 'export')
  outfile = new_dir + '/' + strmid(file_basename(filename),0,pos) + 'export' + fxpar(hdr, 'DATE') + '_im.fits'
  file_move, tmpname, outfile, overwrite = overwrite
  fxaddpar, hdr, 'FILENAME', file_basename(outfile)
 
  ;; Update CHECKSUM and write the header
  fits_add_checksum, hdr
  red_fitscube_newheader, outfile, hdr
  


  ;; Copy/rename related files, like movies and thumbnails.
  other_files = file_search(old_dir+'/'+file_basename(filename, '_im.fits')+'*', count = Nfiles)
  for ifile = 0, Nfiles-1 do begin
    if filename eq other_files[ifile] then continue ; This is not *other* file
    splt = strsplit(other_files[ifile], '.', /extract)
    extension = splt[-1]
    new_other = new_dir + '/' + file_basename(outfile, '.fits') + '.' + extension
    if keyword_set(remove_original) then begin
      ;; Rename
      file_move, other_files[ifile], new_other, overwrite = overwrite
    endif else begin
      ;; Copy
      file_copy, other_files[ifile], new_other, overwrite = overwrite
    endelse
  endfor
  
  
  if keyword_set(remove_original) then begin
    
    ;; Remove the old FITS file
    file_delete, filename
    
    ;; Add to lists of removed and new files so we can update info to
    ;; archive databases.
  endif
  

  
end


oldfile = '/storage/science_data/2023-10-20/CRISP/nb_6173_2023-10-20T10:51:48_10:51:48=0-5_stokes_corrected_export2023-12-25T11:22:29_im.fits'

 red_fitscube_update_keywords, oldfile $
                               , do_update = do_update $
                               , out_dir = out_dir 

end
