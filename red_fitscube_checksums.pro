; docformat = 'rst'

;+
; Add (or check) DATASUM and CHECKSUM keywords for all HDUs in a
; fitscube file.
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
; :Params:
; 
;    filename : in, type=string
; 
;       The name of the file.
; 
; 
; :Keywords:
; 
;    hdus : in, optional, type=array
;
;      Positive integer array specifying which HDUs to check. Default
;      is to check all HDUs, including Main (0).
; 
; 
; :History:
; 
;    2020-07-13 : MGL. First version.
; 
;-
pro red_fitscube_checksums, filename $
                            , hdus = hdus

  ;; Name of this subprogram
  inam = red_subprogram(/low, calling = inam1)

  main_hdr = headfits(filename)
  fits_open, filename, fcb, /update 

  Nhdus = N_elements(hdus)
  if Nhdus EQ 0 then begin 
    Nhdus = fcb.nextend + 1     ; +1 to include the main HDU = 0
    hdus = indgen(Nhdus)
  endif else begin 
    if max(hdus) GT fcb.nextend then $
       print, inam +  ' : HDU ' + strtrim(max(hdus),2) + ' does not exist'
    return
  endelse

  for ihdu = 0, Nhdus-1 do begin

    if hdus[ihdu] eq 0 then hduname = 'Main' else hduname = fcb.extname[hdus[ihdu]]
    
    red_progressbar, ihdu, Nhdus, 'Processing checksums of HDU #'+strtrim(hdus[ihdu], 2)+' : '+hduname
    
    ;; Test the checksum/datasum status. Do the keywords exist and are
    ;; they correct?
    
    if hdus[ihdu] eq 0 then begin

      ;; Main HDU - don't read the data, can be really large
      ;; We need an incremental check of DATASUM for the main HDU. The
      ;; following line only checks the header for corruption.

      dsum = strtrim(red_fitscube_datasum(filename), 2)
      datasum = strtrim(fxpar(main_hdr, 'DATASUM', count = Ndatasum),2)
      if Ndatasum eq 0 then begin
        ;; No DATASUM keyword
        checksum_status = 0     
        ;; Setting checsum_status=0 here *could* hide a faulty
        ;; CHECKSUM keyword and hence file corruption. But we
        ;; shouldn't add one unless we first add DATASUM.
      endif else begin
        if dsum ne datasum then begin
          ;; Incorrect DATASUM keyword
          checksum_status = -1
        endif else begin
          ;; DATASUM keyword present and correct
          checksum_status = fits_test_checksum(main_hdr, errmsg = errmsg, /trust_datasum)
        endelse
      endelse
      ext_header = main_hdr

    endif else begin

      ;; All other HDUs

      fits_read, fcb, ext_data, ext_header, /no_pdu, /NoSCALE, exten = hdus[ihdu]
      checksum_status = fits_test_checksum(ext_header, ext_data, errmsg = errmsg)

    endelse

    ;; Now take action depending on the checksum/datasum status.

    case checksum_status of

      1 : begin

        ;; CHECKSUM (and DATASUM) keywords are present with correct
        ;; values. Nothing to do.
        
      end

      0 : begin
        
        ;; Missing DATASUM or CHECKSUM keyword. Assume no corruption
        ;; and that they just weren't added yet. Calculate and add
        ;; DATASUM and CHECKSUM keywords.
        
        if hdus[ihdu] eq 0 then begin
          ;; Main HDU 
          datasum = red_fitscube_datasum(filename) ; Calculate datasum frame by frame
          fxaddpar, ext_header, 'DATASUM', datasum
          fits_add_checksum, ext_header ; Add 
          modfits, fcb, 0, ext_header
        endif else begin
          ;; All other HDUs
          fits_add_checksum, ext_header, ext_data
          modfits, fcb, 0, ext_header, exten_no = hdus[ihdu]
        endelse

      end
      
      -1 : begin

        ;; CHECKSUM or DATASUM keyword does not have the correct value
        ;; indicating possible data corruption.

        print, inam + ' : CHECKSUM or DATASUM keyword does not have the correct value.'
        print, inam + ' : This file may be corrupt: ' + filename
        print, 'Error message from fits_test_checksum:',  errmsg

        ;; Override:
        ;; fits_add_checksum, ext_header
        ;; modfits, fcb, 0, ext_header, exten_no = hdus[ihdu]
        ;; fits_close, fcb
        ;; retall
        
        stop
      end
      
    endcase
    
  endfor                        ; ihdu
  
  fits_close, fcb
  
end


filename = 'cubes_raw/wb_3950_2016-09-19T09:28:36_scans=3-8_raw_im.fits'

filename = 'cubes_wb/wb_3950_2016-09-19T09:28:36_scans=0-9_meanang_corrected_im.fits'

red_fitscube_checksums, filename
;red_fitscube_add_checksums, filename, hdus = 0

end

h = headfits(filename)
print, fits_test_checksum(h, errmsg = errmsg, /trust)
fits_add_checksum, h
print, fits_test_checksum(h, errmsg = errmsg, /trust)

end
