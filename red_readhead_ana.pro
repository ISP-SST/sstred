; docformat = 'rst'

;+
; Return a fits-type header from an ANA fz format file.
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, ISP
; 
; 
; :Returns:
; 
; 
; :Params:
; 
;    fname : in, type=string
;
;       The name of the data file.
;
; :Keywords:
;
;    date_beg : out, type=strarr
;
;       The timestamps for exposure start.
; 
;    framenumbers : out, type=intarr
;
;       The framenumbers extracted from the file metadata.
;    
; 
; :History:
; 
;    2017-03-10 : MGL. Moved reading of headers from ANA fz format
;                 files from red_readhead.pro.
; 
;    2018-06-04 : MGL. New keywords date_beg and framenumbers.
; 
;    2019-04-05 : MGL. Import code from red_anahdr2fits.
; 
;-
function red_readhead_ana, fname, $
                           date_beg = date_beg, $
                           framenumbers = framenumbers

  compile_opt idl2

  anaheader = fzhead(fname)

  if n_elements(anaheader) eq 0 then begin
    print, 'No ana header in '+fname
    return, ''
  endif 

  
  if strmatch( anaheader, "SIMPLE*" ) gt 0 then begin
    
    ;; It's already a fits-type header, split into strarr with length 80
    len = strlen(anaheader)
    
    header = strmid(anaheader,indgen(len/80 + (len MOD 80 ne 0))*80,80)

    framenumbers = fxpar(header, 'FRAMENUM', count = count)
    if count eq 0 then undefine, framenumbers

    date_beg = fxpar(header, 'DATE-BEG', count = count)
    if count eq 0 then undefine, date_beg

    header = red_meta2head(header, meta={filename:fname})
    
  endif else begin

    filename  = file_basename(fname)
    directory = file_dirname(fname)

    timeregex = '[0-2][0-9]:[0-5][0-9]:[0-6][0-9]'
    dateregex = '20[0-2][0-9][.-][01][0-9][.-][0-3][0-9]'
    
    ;; Do what we can with a "raw" fz file header.

    ;; Get size info from the file
    openr, llun, fname, /get_lun
    ;; Read the first few bytes to get data type and number of
    ;; dimensions:
    ah = bytarr(192) 
    readu, llun, ah
    case ah[7] of
      0: datatype = 1           ; ANA_BYTE
      1: datatype = 2           ; ANA_WORD
      2: datatype = 3           ; ANA_LONG
      3: datatype = 4           ; ANA_FLOAT
      4: datatype = 5           ; ANA_DOUBLE
      else: datatype = 2        ; default?
    endcase

    ;; Read bytes 192-255 as longs to get the dimensions
    bh = lonarr(16)
    readu, llun, bh
    naxis = bh[0:ah[8]-1]

    ;; Start making a header
    red_mkhdr, header, datatype, naxis
    anchor = 'DATE'

    red_fitsaddkeyword, anchor = anchor, header, 'FILENAME', filename, ''
    
    ;; Close the file
    free_lun, llun

    ;; A clue as to what kind of ana file this is:
    tspos = strpos(anaheader, 'Ts=')

    if tspos ne -1 then begin

      ;; This is probably a raw data file written with Michiel's
      ;; acquisition code

      date_obs = stregex(directory, dateregex, /extract) + 'T' $
                 + stregex(directory, timeregex, /extract)
      if date_obs ne 'T' then begin
        red_fitsaddkeyword, anchor = anchor, header $
                            , 'DATE-OBS', date_obs, 'Data set timestamp'
      endif else begin
        ;; If there is a properly formatted date+'T'+timestamp in the
        ;; file name, it should be the date-obs. I.e., corresponding
        ;; to the timestamp of the data collection directory.
        date_obs = stregex(filename, dateregex+'T'+timeregex, /extract) 
        if date_obs ne '' then red_fitsaddkeyword, anchor = anchor, header $
           , 'DATE-OBS', date_obs $
           , 'Data set timestamp'
      endelse
      

      
      ;; We have already found the starting time
      date_beg = strmid(anaheader, tspos+3, 26)
      strput,date_beg,'-',4
      strput,date_beg,'-',7
      strput,date_beg,'T',10
      red_fitsaddkeyword, anchor = anchor, header $
                          ,'DATE', date_beg, 'Creation of fz file'
      red_fitsaddkeyword, anchor = anchor, header $
                          ,'DATE-BEG', date_beg, 'Begin exposure'

      ;; There should be an ending time as well but don't just assume
      ;; there is.
      tepos = strpos(anaheader, 'Te=')
      if tepos ne -1 then begin
        date_end = strmid(anaheader, tepos+3, 26)
        strput,date_end,'-',4
        strput,date_end,'-',7
        strput,date_end,'T',10
        red_fitsaddkeyword, anchor = anchor, header $
                            ,'DATE-END', date_end, 'End exposure' 
      endif
      
      if tspos ne -1 and tepos ne -1 then begin
        exptime = red_time2double(strmid(date_end, 11)) $
                  - red_time2double(strmid(date_beg, 11))
        red_fitsaddkeyword, anchor = anchor, header $
                            , 'XPOSURE', exptime, '[s] Exposure length'
      end

      
      ;; Determine the instrument
      tmp = stregex(anaheader, 'Camera ([IVXL]*) \[CRISP-([WRTD])\]', /extract, /subexpr)
      case 1 of

        tmp[0] ne '' : begin
          ;; CRISP data
          red_fitsaddkeyword, anchor = anchor, header $
                              , 'INSTRUME', 'CRISP', ' Name of instrument'
          ;; Find DETECTOR and CAMERA
          if n_elements(tmp) eq 3 then begin
            red_fitsaddkeyword, anchor = anchor, header $
                                , 'DETECTOR', 'cam'+tmp[1], 'Detector identifier'
            camera = 'Crisp-'+tmp[2] 
            red_fitsaddkeyword, anchor = anchor, header $
                                , 'CAMERA', camera
          endif
        end
        
        strmatch(anaheader,'*"Spectrograph*"*') : begin
          ;; TRIPPEL data 
          red_fitsaddkeyword, anchor = anchor, header $
                              , 'INSTRUME', 'TRIPPEL', ' Name of instrument'
          ;; Find DETECTOR 
          campos = strpos(anaheader, '"Camera')
          if campos ne -1 then begin
            detector = 'cam' + (strsplit(strmid(anaheader, campos+8), ' ', /extract))[0]
            red_fitsaddkeyword, anchor = anchor, header $
                                , 'DETECTOR', detector, 'Camera identifier'
          end
        end

        else : begin
          ;; Instrument unknown. What do slitjaw and other supporting
          ;; TRIPPEL data look like? Need to get the detector from the
          ;; file name? Are there defined camera names?          
        end
      endcase

      ;; The frame number is usually the last part of the file name.
      ;; Seven digits long.
      lastpart = (strsplit(fname, '[._]', /extract))[-1]
      if strlen(lastpart) eq 7 then begin
        framenumbers = long(lastpart)
        red_fitsaddkeyword, anchor = anchor, header $
                            , 'FRAMENUM', framenumbers, 'Frame number'
      end

      ;; Scan numbers are 5 digits long.
      scannumber = ((stregex(filename, '(_|\.|^)([0-9]{5})(_|\.|$)', /extr, /subexp))[2,*])[0]
      if scannumber ne '' then red_fitsaddkeyword, anchor = anchor, header $
         , 'SCANNUM', long(scannumber) $
         , 'Scan number'

      lcstate = ((stregex(filename, '(\.)lc(.?)(\.)', /extr, /subexp))[2,*])[0]
      if lcstate ne '' then red_fitsaddkeyword, anchor = anchor, header $
                                                , 'LCSTATE', long(lcstate) $
                                                , 'Liquid crystal state'

      lpstate = ((stregex(filename, '(\.)LP(.?.?.?)(\.)', /extr, /subexp))[2,*])[0]
      if lpstate ne '' then red_fitsaddkeyword, anchor = anchor, header $
                                                , 'LPSTATE', long(lpstate) $
                                                , 'Linear polarizer state'
 
      qwstate = ((stregex(filename, '(\.)qw(.?.?.?)(\.)', /extr, /subexp))[2,*])[0]
      if qwstate ne '' then red_fitsaddkeyword, anchor = anchor, header $
                                                , 'QWSTATE', long(qwstate) $
                                                , 'Quarterwave plate state'
      
      
      tuninfo = stregex(filename, '([0-9][0-9][0-9][0-9])[._]([0-9][0-9][0-9][0-9])_([+-][0-9]*)' $
                        , /extract, /subexpr)
      
      if tuninfo[0] ne '' then begin
        prefilter = strtrim(tuninfo[1], 2)
        red_fitsaddkeyword, anchor = anchor, header $
                            , 'FILTER1', prefilter $
                            , 'CRISP prefilter'
        if camera eq 'Crisp-T' or camera eq 'Crisp-R' then begin
          ;; NB
          red_fitsaddkeyword, anchor = anchor, header $
                              , 'STATE', tuninfo[0] $
                              , 'NB tuning state'        
        endif else begin
          ;; WB
          red_fitsaddkeyword, anchor = anchor, header $
                              , 'STATE', prefilter+'_'+prefilter+'_+0000' $
                              , 'WB tuning state'
        endelse
        ;; Info inferred from the tuning info
        case prefilter of
          '5173' : begin        ; CRISP
            wavelnth = 517.3e-9
            fwhm = 0.30e-9
            waveband = 'Mg b 5173'
          end
          '5250' : begin        ; CRISP
            wavelnth = 525.0e-9
            fwhm = 0.30e-9
            waveband = 'Fe I 5250'
          end
          '5382' : begin        ; CRISP
            wavelnth = 538.2e-9
            fwhm = 0.33e-9
            waveband = 'C I 5382'
          end
          '5578' : begin        ; CRISP
            wavelnth = 557.8e-9
            fwhm = 0.30e-9
            waveband = 'Fe I 5578'
          end
          '5897' : begin        ; CRISP
            wavelnth = 589.7e-9
            fwhm = 0.38e-9
            waveband = 'Na D 5897'
          end
          '6173' : begin        ; CRISP - Alluxa filter
            wavelnth = 617.38e-9
            fwhm = 0.45e-9
            waveband = 'Fe I ' + strtrim(prefilter,2)
          end
          '6174' : begin        ; CRISP
            wavelnth = 617.4e-9
            fwhm = 0.43e-9
            waveband = 'Fe I 6174'
          end
          '6302' : begin        ; CRISP
            wavelnth = 630.26e-9
            fwhm = 0.45e-9
            waveband = 'Fe I 6301+6302' 
          end
          '6563' : begin        ; CRISP
            wavelnth = 656.38e-9
            fwhm = 0.49e-9
            waveband = 'H-alpha 6563' 
          end
          '7772' : begin        ; CRISP
            wavelnth = 777.25e-9
            fwhm = 0.7e-9
            waveband = 'O I 7772'
          end
          '8542' : begin        ; CRISP
            wavelnth = 854.1e-9
            fwhm = 0.93e-9
            waveband = 'Ca II 8542'
          end
          else :
        endcase
        if n_elements(wavelnth) ne 0 then begin
          red_fitsaddkeyword, anchor = anchor, header $
                              , 'WAVEUNIT', -9, 'Wavelength unit 10^WAVEUNIT m = nm'
          red_fitsaddkeyword, anchor = anchor, header $
                              , 'WAVELNTH', wavelnth*1e9, '[nm] Prefilter peak wavelength'
          red_fitsaddkeyword, anchor = anchor, header $
                              , 'WAVEMIN', (wavelnth-fwhm/2)*1e9, '[nm] Prefilter min wavelength (0.5 peak)'
          red_fitsaddkeyword, anchor = anchor, header $
                              , 'WAVEMAX', (wavelnth+fwhm/2)*1e9, '[nm] Prefilter max wavelength (0.5 peak)'
          red_fitsaddkeyword, anchor = anchor, header $
                              , 'WAVEBAND', waveband
        endif
      endif else begin

        ;; File name does not have a filter_line_tuning part. Some
        ;; polcal data is like that. Look for just the prefilter tag
        ;; in this case.
        prefilter = ((stregex(filename, '(_|\.|^)([0-9]{4})(_|\.|$)', /extr, /subexp))[2,*])[0]
        if prefilter ne '' then $
           red_fitsaddkeyword, anchor = anchor, header $
                               , 'FILTER1', prefilter $
                               , 'CRISP prefilter'
        
      endelse

    
    
      ;; Should extract more info from anaheader: states of prefilter,
      ;; liquid crystals, LRE, and HRE. But first find out what
      ;; keywords to use for them in the FITS header.
      
    endif else begin
      
      ;; Assume this is fz format output from momfbd
      
      ;; Header date 
      dpos = strpos(anaheader, 'DATE=')
      if dpos ne -1 then begin
            date = (red_strreplace(strmid(anaheader, dpos+5, 19), ' ', 'T'))[0]
            red_fitsaddkeyword, header, 'DATE', date, ''
      endif

      case 1 of
        strmatch(anaheader,'*CRISP-W*',/FOLD) : red_fitsaddkeyword, anchor = anchor, header, 'CAMERA', 'Crisp-W' 
        strmatch(anaheader,'*CRISP-R*',/FOLD) : red_fitsaddkeyword, anchor = anchor, header, 'CAMERA', 'Crisp-R' 
        strmatch(anaheader,'*CRISP-T*',/FOLD) : red_fitsaddkeyword, anchor = anchor, header, 'CAMERA', 'Crisp-T' 
        strmatch(anaheader,'*CRISP-D*',/FOLD) : red_fitsaddkeyword, anchor = anchor, header, 'CAMERA', 'Crisp-D'
        else:
      endcase
      
      ;; Observations date
      dpos = strpos(anaheader, 'DATE_OBS')
      if dpos ne -1 then begin
        date_obs = strmid(anaheader, dpos+9, 10)
        ;; Would like to add a time but TIME_OBS is usually empty
        tpos = strpos(anaheader, 'TIME_OBS')
        if tpos ne -1 then begin
          time_obs = strmid(anaheader, tpos+9,dpos-(tpos+9))
          if strlen(time_obs) gt 1 then date_obs += 'T' + time_obs
        endif 
        red_fitsaddkeyword, anchor = anchor, header, 'DATE-AVG', date_obs, '', after = 'DATE'
      endif

      header = red_meta2head(header, meta={filename:fname})
      
    endelse

  endelse

  
  red_fitsaddkeyword, header, 'SOLARNET', 0.5,  format = 'f3.1' $
                      , 'Fully SOLARNET-compliant=1.0, partially=0.5', before = 'DATE'


  ;; Get what other info is possible from the anaheader:
  ;; header = red_anahdr2fits(anaheader, datatype = datatype, naxisx = naxisx)


  return, header

end




