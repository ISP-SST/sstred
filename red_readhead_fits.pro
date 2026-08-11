; docformat = 'rst'

;+
; Return the header from a FITS format file, taking some
; pipeline-specific issues into account. 
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
;    select_frame : in, type=int
;
;       Only get meta-data for a specific frame. (TODO: implement)
;
;    status : out, optional, type=integer
;
;       0 is succeeded, -1 otherwise.
; 
; :History:
; 
;    2017-03-10 : MGL. Moved reading of headers from ANA fz format
;                 files from red_readhead.pro.
; 
;    2017-03-13 : MGL. Deal with camera software OBS_PHDU bug.
; 
;    2017-07-06 : THI. Use rdx_readhead to get header.
; 
;    2017-09-01 : THI. Get date_beg and framenumbers from file.
; 
;    2017-12-01 : MGL. New keyword status, use status from rdx_readhead.
;
; 
;-
function red_readhead_fits, fname, $
                            date_beg = date_beg, $
                            extension = extension, $
                            framenumbers = framenumbers, $
                            select_frame = select_frame, $
                            silent = silent, $
                            status = status

  compile_opt idl2

  if ~file_test(fname) then begin
    status = -1
    return, 0
  endif
  
  ;; Are we here for an extension or the primary HDU?
  if n_elements(extension) eq 0 then begin
    ;; primary

    header = rdx_readhead(fname,date_beg=date_beg,framenumbers=framenumbers, status = status)

    fxaddpar, header, 'TIMESYS', 'UTC', after = 'SOLARNET'    

    if status ne 0 then return, 0
    
    pref = fxpar( header, red_keytab('pref'), count=Npref )

    if Npref eq 1 then begin

      ;; The prefilter is given in the header as FILTER1, all is well
      prefilter = pref

    endif else begin

      ;; Try to find out the prefilter from other information

      camera = fxpar(header, red_keytab('camera'), count=Ncamera)
      state = fxpar(header, 'STATE', count=Nstate )

      if Nstate gt 0 && state eq '' then begin
        ;; We have no use for an empty STATE
        sxdelpar, header, 'STATE'
        Nstate = 0
      endif
      
      ;; Need to re-implement getting it from the file name in the
      ;; case Nstate eq 0? Ignore for now!
      
      if Nstate gt 0 then begin
        
        
        ;; So now we know there is a state keyword. But does it have
        ;; the tuning info in readable form?
        state_split = strsplit( state, '_',  /extr )
        
        ;; Automatic-mosaic data have state strings that start with
        ;; "mosNN", where "NN" is a two-digit number. Remove it!
        if strmatch(state_split[0], 'mos[0-9][0-9]') then begin
          state_split = state_split[1:*]
          fxaddpar, header, 'STATE', strjoin(state_split, '_'), 'Shortened by red_readhead_fits'  
        endif
        
        Nsplit = n_elements(state_split)        

        state1 = state_split[0]
        
        if strmatch(state1,'[0-9][0-9][0-9][0-9]') then begin

          ;; The initial segment is four digits, has to be the
          ;; prefilter tag!
          prefilter = state1

        endif else begin
          
          ;; This is probably chromis data with the state given in the
          ;; something like 'wheel00002_hrz32600'. Or it could be
          ;; CRISP darks.

          case camera of

            'Chromis-N' : begin
              case state1 of
                'wheel00001' : prefilter = '3925' ; Ca II K blue wing
                'wheel00002' : prefilter = '3934' ; Ca II K core
                'wheel00003' : prefilter = '3969' ; Ca II H core
                'wheel00004' : prefilter = '3978' ; Ca II H red wing
                'wheel00005' : prefilter = '3999' ; Ca II H continuum
                'wheel00006' : prefilter = '4862' ; H-beta core
                else :
              endcase
            end

            'Chromis-W' : begin
              case state1 of
                'wheel00006' : prefilter = '4846' ; H-beta continuum
                else: prefilter = '3950'          ; Ca II HK wideband
              endcase
            end
            
            'Chromis-D' : begin
              case state1 of
                'wheel00006' : prefilter = '4846' ; H-beta continuum
                else: prefilter = '3950'          ; Ca II HK wideband
              endcase
            end

            else : begin
              ;; Could be CRISP darks with an lc tag in the STATE
              ;; keyword. Remove if that's the case.
              if strmatch(state, 'lc?') then sxdelpar, header, 'STATE'
            end
          endcase
        endelse
        
        ;; LCSTATE, LPSTATE, QWSTATE
        case Nsplit of

          6 :  begin
            ;; Polcal data, state_split[3:*] = [LP,QW,LC]
            lp = long(strmid(state_split[3], 2))
            qw = long(strmid(state_split[4], 2))
            lc = long(strmid(state_split[5], 2))
            fxaddpar, header, 'LPSTATE', lp, 'Extracted from state keyword'            
            fxaddpar, header, 'QWSTATE', qw, 'Extracted from state keyword'            
            fxaddpar, header, 'LCSTATE', lc, 'Extracted from state keyword'           
            fxaddpar, header, 'STATE', strjoin(state_split[0:2], '_'), 'Shortened by red_readhead_fits'            
          end

          5 :  begin
            ;; Polcal data, state_split[3:*] = [LP,QW], LCSTATE
            ;; already in header
            lp = long(strmid(state_split[3], 2))
            qw = long(strmid(state_split[4], 2))
            fxaddpar, header, 'LPSTATE', lp, 'Extracted from state keyword'            
            fxaddpar, header, 'QWSTATE', qw, 'Extracted from state keyword'            
            fxaddpar, header, 'STATE', strjoin(state_split[0:2], '_'), 'Shortened by red_readhead_fits'            
          end
          
          4 : begin
            ;; Other polarimetric data, state_split[3] = [LC]
            lc = long(strmid(state_split[3], 2))
            fxaddpar, header, 'LCSTATE', lc, 'Extracted from state keyword'  
            fxaddpar, header, 'STATE', strjoin(state_split[0:2], '_'), 'Shortened by red_readhead_fits'  
          end

          else :
          
        endcase
        
      endif 

    endelse 
  
    ;; Now write the prefilter to the header, if we found one.
    ;; (It's not defined for darks!)
    if n_elements(prefilter) gt 0 then begin
      fxaddpar, header, red_keytab('pref'), prefilter, 'Extracted from state keyword'
    endif

    ;; Now rewrite the STATE keyword, the pipeline expects it not to
    ;; have the LP,QW,LC tags
;    
;    ;; Hack to get the prefilter from the file name in data
;    ;; from 2016.08.30.
;    pref = fxpar( header, red_keytab('pref'), count=count )
;    if count eq 0 then begin
;      state = fxpar(header, 'STATE', count=count )
;      if count eq 0 then begin
;        ;; Try to read from file name
;        fname_split = strsplit(file_basename(fname,'.fits'),'_',/extr)
;        if fname_split[-1] eq 'cavityfree.flat' then begin
;          prefilter = (fname_split)[1]
;          fxaddpar, header, red_keytab('pref'), prefilter, 'Extracted from file name'
;        endif else if n_elements(fname_split) gt 4 then begin
;          ;; Shorter and it might be a dark
;          prefilter = (fname_split)[-1]
;          ;; Translate to previously used filter names
;          case prefilter of
;            'hbeta-core' : prefilter = '4862'
;            'hbeta-cont' : prefilter = '4846'
;            'cah-core'   : prefilter = '3969'
;            else:
;          endcase
;          fxaddpar, header, red_keytab('pref'), prefilter, 'Extracted from file name'
;        endif
;      endif else begin          ; STATE keyword exists but not FILTER1 (e.g. 2016.08.31 data)
;
;        state_split = strsplit( state, '_',  /extr )
;
;        if n_elements(state_split) gt 1 then begin
;          
;          if strmatch(state1,'[0-9][0-9][0-9][0-9]') then prefilter = state1 else begin          
;
;            state1 = state_split[0] ;  for 2016.08.31:  state = 'wheel00002_hrz32600'
;            camera = fxpar(header, red_keytab('camera'), count=count)
;
;            if count gt 0 then begin
;              if camera eq 'Chromis-N' then begin
;                
;                ;; Chromis-N
;                case state1 of
;                  'wheel00001' : prefilter = '3925' ; Ca II K blue wing
;                  'wheel00002' : prefilter = '3934' ; Ca II K core
;                  'wheel00003' : prefilter = '3969' ; Ca II H core
;                  'wheel00004' : prefilter = '3978' ; Ca II H red wing
;                  'wheel00005' : prefilter = '3999' ; Ca II H continuum
;                  'wheel00006' : prefilter = '4862' ; H-beta core
;                  else :
;                endcase
;              endif else begin
;                ;; Chromis-W and Chromis-D
;                case state1 of
;                  'wheel00006' : prefilter = '4846' ; H-beta continuum
;                  else: prefilter = '3950'          ; Ca II HK wideband
;                endcase
;              endelse
;            endif
;          endelse
;          if n_elements(prefilter) gt 0 then begin
;            ;; Not defined for darks
;            fxaddpar, header, red_keytab('pref'), prefilter, 'Extracted from state keyword'
;          endif
;        endif else begin
;          
;          state1 = state_split[0] ;  for 2016.08.31:  state = 'wheel00002_hrz32600'
;          camera = fxpar(header, red_keytab('camera'), count=count)
;
;          case state of
;            'hbeta-core' : prefilter = '4862'
;            'hbeta-cont' : prefilter = '4846'
;            'cah-core'   : prefilter = '3969'
;            'wheel00005' : prefilter = '3950' ; WB: Ca II HK continuum
;            'wheel00006' : prefilter = '4846' ; WB: H-beta continuum
;            else:
;          endcase
;          if n_elements(prefilter) gt 0 then begin
;            ;; Not defined for darks
;            fxaddpar, header, red_keytab('pref'), prefilter, 'Extracted from state keyword'
;          endif
;        endelse
;      endelse
;    endif


    ;; Correct OBS_PHDU bug in old version of camera software
    OBS_PHDU = sxpar( header, 'OBS_PHDU', comment = pcomment, count=count )
    if count gt 0 then begin
      if n_elements(pcomment) eq 0 || strtrim(pcomment,2) eq '' then $
         pcomment = ' Observational Header and Data Unit'
      fxaddpar, header, 'OBS_HDU', OBS_PHDU, pcomment, after = 'OBS_PHDU'
      sxdelpar, header, 'OBS_PHDU'
    endif

    ;; Correct also for OBS_SHDU, which used to be SOLARNET standard.
    OBS_SHDU = sxpar( header, 'OBS_SHDU', comment = scomment, count=count )
    if count gt 0 then begin
      if n_elements(scomment) eq 0 || strtrim(scomment,2) eq '' then $
         scomment = ' Observational Header and Data Unit'
      fxaddpar, header, 'OBS_HDU', OBS_SHDU, scomment, after = 'OBS_SHDU'
      sxdelpar, header, 'OBS_SHDU'
    endif


    if n_elements(select_frame) ne 0 then begin
      ;; We may want to change or remove some header keywords
      ;; here, like FRAME1, CADENCE, and DATE-END.
    endif
  endif else begin
    
    ;; EXTENSION header
    elun = fxposit(fname,extension,/readonly,/no_fpack)
    fxhread,elun,header
    free_lun,elun

    status = 0                  ; Set to -1 if we detect a problem!
    
  endelse

  return, header

end



;; Test compressed file
fname = '/nadir-scratch/tomas/wfwfs/ffov_12_compressed/sst_camXXXVI_00000_0000000.fits'
header = red_readhead_fits(fname,date_beg=date_beg,framenumbers=framenumbers, status = status)
h = headfits(fname)
hh = red_readhead(fname,date_beg=date_beg,framenumbers=framenumbers, status = status)
stop

;; Test broken file

dir='/storage/sand05n/Incoming/2017.04.20/CHROMIS-flats/18:34:10/Chromis-N/'
dir='/data/2017/2017.04/2017.04.20/CHROMIS-flats/18:34:10/Chromis-N/'
fname = 'sst_camXXX_00004_0036400_wheel00006_hrz34410.fits' ; Broken file
fname = 'sst_camXXX_00004_0034000_wheel00006_hrz33138.fits' ; OK file

header = red_readhead_fits(dir+fname,date_beg=date_beg,framenumbers=framenumbers, status = status)
print, status

end
