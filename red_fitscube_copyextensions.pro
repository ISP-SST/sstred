; docformat = 'rst'

;+
; Copy extensions from one fitscube file to another.
;
; No check for unique extension names is done. If an extension name is
; asked for, all extensions with that name are copied (or ignored). It
; doesn't matter if an extension name is asked for multiple times.
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
;    infile : in, type=string
; 
;       The name of the original file, from which to copy the
;       extensions.
; 
;    outfile : in, type=string
; 
;       The name of the target file, to which the extensions should be
;       copied.
; 
; 
; :Keywords:
; 
;     anchor : in, out, optional, type=string
;
;        Passed on to red_fitscube_addvarkeyword for extensions that
;        are variable keywords.
; 
;     ext_list : in, out, optional, type=array, default='All extensions' 
;
;        Names or indices of extensions to be copied (default) or
;        ignored (if /ignore keyword is set). 
; 
;     ext_regex : in, optional, type=strarr
; 
;        Add to ext_list all existing extensions that match any of
;        these regular expressions.
; 
;     ext_statistics : in, optional, type=boolean
; 
;        Add to ext_reqex regular expressions that match extensions
;        corresponding to statistics variable keywords.
; 
;     ignore : in, optional, type=boolean 
;
;        Interpret ext_list as a list of extensions to ignore.
; 
; :History:
; 
;    2020-10-16 : MGL. First version.
; 
;    2022-03-30 : MGL. New keywords ext_regex and ext_statistics.
; 
;-
pro red_fitscube_copyextensions, infile, outfile $
                                 , anchor = anchor $
                                 , ext_list = ext_list $
                                 , ext_regex = ext_regex $
                                 , ext_statistics = ext_statistics $
                                 , ignore = ignore 


  ;; Name of this subprogram
  inam = red_subprogram(/low, calling = inam1)

  main_hdr = headfits(infile)
  fits_open, infile, fcb_in

  
  if keyword_set(ext_statistics) then red_append, ext_regex, ['VAR-EXT-DATA*', 'VAR-EXT-NPIXELS']
  
  for iregex = 0, n_elements(ext_regex)-1 do begin
    indx = where(strmatch(fcb_in.extname, ext_regex[iregex]), Nmatch)
    if Nmatch gt 0 then begin
      red_append, ext_list, fcb_in.extname[indx]
    endif
  endfor                        ; iregex

  ;; Don't repeat items in ext_list
  if n_elements(ext_list) gt 0 then ext_list = ext_list[uniq(ext_list,sort(ext_list))]
  
  Nlist = n_elements(ext_list)  
  
  if Nlist eq 0 then begin

    ;; Default: copy all extensions

    Next = fcb_in.Nextend 
    ext_numbers = indgen(Next) + 1

  endif else begin

    ;; The ext_list keyword was given. First check it and convert to
    ;; numbers if needed.
    
    case !true of

      isa(ext_list, /integer) : begin
        ;; Integers are interpreted as extension numbers, this is what
        ;; we want.

        ;; In range?
        if max(ext_list) gt fcb_in.Nextend or min(ext_list) lt 1 then begin
          print, inam + ' : At least one extension number in ext_list is out of range [1,' $
                 + strtrim(fcb_in.Nextend) + ']:'
          print, ext_list
          stop
        endif

        ;; Unique?
        
        ext_numbers = ext_list

      end

      isa(ext_list, /string) : begin

        ;; Match ext_list with the existing extension names 
        
        match2, ext_list, fcb_in.extname, suba, subb

        ;; Any non-existing extensions asked for?
        if min(suba) le 0 then begin
          print, inam + ' : At least one extension name in ext_list is not in list of existing extensions:'
          print, fcb_in.extname
          print
          print, ext_list
          stop
        endif
        
        ;; We pick the indices i, where fcb_in.extname[i] matches any
        ;; of the names in ext_list.
        ext_numbers = where(subb ne -1)        
        
      end

      else : stop               ; Anything but integers or strings that makes sense?
      
    endcase 

    if keyword_set(ignore) then begin
      ;; Ignore extensions in ext_list
      match2, ext_numbers, indgen(fcb_in.Nextend) + 1, suba, subb
      ext_numbers = where(subb eq -1, Nwhere) + 1
      if Nwhere eq 0 then stop
    endif
    
    Next = n_elements(ext_numbers)

  endelse

  
;  if n_elements(ext_list) gt 0 then print, ext_list
;  print
;  print, ext_numbers
;  print
;  print, fcb_in.extname[ext_numbers]
;  
;  stop
  
  ;; Loop through the extension numbers, copying the extensions one by
  ;; one. 
  for iext = 0, Next-1 do begin

    red_progressbar, iext, Next, 'Copying extensions : #' $
                     + strtrim(ext_numbers[iext], 2) + ', ' $
                     + fcb_in.extname[ext_numbers[iext]]
    
    case !true of

      strmid(fcb_in.extname[ext_numbers[iext]], 0, 8) eq 'VAR-EXT-' : begin
        ;; For a variable keyword, there is a BINTABLE extension to
        ;; copy but we need to also update the VAR_KEYS main header
        ;; keyword. This is taken care of by
        ;; red_fitscube_addvarkeyword.
        var_key = strmid(fcb_in.extname[ext_numbers[iext]], 8)
        red_fitscube_addvarkeyword, outfile, var_key $
                                    , anchor = anchor $
                                    , old_filename = infile
      end

      fcb_in.xtension[ext_numbers[iext]] eq 'BINTABLE' : begin
        ;; Other BINTABLE extension
        red_fits_copybinext, infile, outfile, ext_numbers[iext]
        ;;  red_fits_copybinext, infile, outfile, fcb_in.extname[ext_numbers[iext]]
      end

      fcb_in.xtension[ext_numbers[iext]] eq 'IMAGE' : begin
        ;; IMAGE extension, e.g., WCSDVARR.
        ext_data = mrdfits(infile, ext_numbers[iext], ext_hdr $
                           , status = status, /silent)
        if status ne 0 then stop
        writefits, outfile, ext_data, ext_hdr, /append
      end

      fcb_in.xtension[ext_numbers[iext]] eq 'TABLE' : begin
        ext_data = mrdfits(infile, ext_numbers[iext], ext_hdr $
                           , status = status, /silent)
        if status ne 0 then stop
        stop             ; When I did this (without ext_hdr!) with raw data files, redux threw an error when trying to read the extension.
        mwrfits, ext_data, ext_hdr, outfile, /ascii, status=outstatus, /silent
      end
      
      else : begin
        ;; Other extensions. Add more CASEs if we get into trouble
        ;; with any particular type.

        ;; Stop because we should be aware if we are using anything
        ;; but BINTABLE, TABLE, or IMAGE extensions.
        stop
        
        ;; Read the extension
        fits_read, fcb_in, ext_data, ext_header $
                 , /no_pdu, /NoSCALE $
                 , exten_no = ext_numbers[iext] $
                 , extver = extver $
                 , xtension = xtension
        help, xtension
        print, fxpar(ext_header, 'XTENSION')

        if 0 then begin
          ;; Checksums?
          checksum_status = fits_test_checksum(ext_header, ext_data $
                                               , errmsg = errmsg)

          if checksum_status eq -1 then begin
            
            ;; CHECKSUM or DATASUM keyword do not have the correct value
            ;; indicating possible data corruption.

            print, inam + ' : CHECKSUM or DATASUM keyword does not have the correct value.'
            print, inam + ' : This extension may be corrupt.'
            print, 'Error message from fits_test_checksum:',  errmsg

            stop

            ;; Remove from list of copied extensions
                                ;           continue            ; Proceed with the next extension
            
          endif
        endif
        
        ;; Does this extension exist in the outfile already? Update the
        ;; existing extension or write another with the same name?
      
        ;; Write the extension
        fits_write, outfile, reform(ext_data, fxpar(ext_header,'NAXIS*')), ext_header

;        fits_open, outfile, fcb_out
;        help,fcb_out
;        fits_read, fcb_in, ext_data_out, ext_header_out $
;                   , /no_pdu, /NoSCALE $
;                   , extname = fcb_in.extname[ext_numbers[iext]] $
;                   , extver = extver_out $
;                   , xtension = xtension_out
;        fits_close, fcb_out
;        
;        print, ext_header + ext_header_out
;        help, ext_data,ext_data_out, extver, extver_out, xtension, xtension_out
;        
;        stop
        
      end

    endcase
    
  endfor                        ; iext
  
  fits_close, fcb_in
  
end


wdir = '/scratch/mats/2016.09.19/CRISP-aftersummer/'
cd, wdir

filename = 'cubes_wb/wb_6302_2016-09-19T09:30:20_scans=2-8_corrected_im.fits'
filename = 'cubes_nb/nb_6302_2016-09-19T09:30:20_scans=2-8_stokes_corrected_im.fits'

file_copy, filename, filename+'.bak', /overwrite


red_fitscube_copyextensions, filename, filename+'.bak'  $
                             , ext_list = ['WCS-TAB'] $
                             , ignore = 1
end


hdr = headfits(filename)
fxaddpar, hdr, 'TEST', 22

red_fitscube_newheader, filename+'.bak', hdr, Nframes_max = 1

hdr_orig = headfits(filename)
hdr_copy = headfits(filename+'.bak')

fits_open, filename,  fcb_orig
fits_open, filename+'.bak', fcb_copy, /update 

hprint, fcb_orig.extname
print
hprint, fcb_copy.extname

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
