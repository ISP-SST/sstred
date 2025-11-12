; docformat = 'rst'

;+
; Find the used geometrical pinhole alignment model from a file.
; 
; The file can be a momfbd cfg file, a momfbd output file, or an
; associated fitsheader file.
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
;   Either "projective" or "polywarp".
; 
; :Params:
; 
;   filename : in, type=string
; 
;     The file to be examined for the used alignment model.
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
;   2025-02-19 : MGL. First version.
; 
;-
function red_align_model, filename

  if ~file_test(filename) then begin
    red_message, 'File does not exist: ' + filename
    return, ''                  ; No such file
  endif
  
  extension = (strsplit(filename, '.', /extract))[-1]

  if extension eq 'cfg' then begin
    ;; Look for the ALIGN_MAP type keywords, if any.
    cfg = redux_readcfg(filename)
    
    ;; The old projective transform model?
    tmp = redux_cfggetkeyword(cfg, 'OBJECT0.CHANNEL0.ALIGN_MAP', count = cnt)
    if cnt gt 0 then return, 'projective'

    ;; The new polynomial warp model?
    tmp = redux_cfggetkeyword(cfg, 'OBJECT0.CHANNEL0.ALIGN_MAP_X', count = cnt)
    if cnt gt 0 then return, 'polywarp'

    ;; No such alignment in the file.
    return, ''            

  endif

  case extension of

    'fitsheader' : hdr = headfits(filename)
 
    'momfbd' : hdr = red_readhead(filename)

    'fits' : hdr = red_readhead(filename)

    else : return, ''

  endcase
  
  ;; Now look in the step info

  ;; Check the MOMFBD step
  prpara = red_headerinfo_getstep(hdr, prstep = 'MOMFBD  ', prkey = 'PRPARA', count = cnt, /silent)
  if cnt eq 0 then begin
    ;; Some old files may have this instead
    prpara = red_headerinfo_getstep(hdr, prstep = 'MOMFBD image restoration', prkey = 'PRPARA', count = cnt, /silent)
    if cnt eq 0 then begin      
      ;; No MOMFBD processing step. This could be output from the
      ;; bypass_momfbd method. In this case the cfg file is in the
      ;; prstep for that.
      prpara = red_headerinfo_getstep(hdr, prstep = 'CONCATENATION,SPATIAL-ALIGNMENT,DESTRETCHING' $
                                      , prkey = 'PRPARA', count = cnt, /silent)
      if cnt gt 0 then begin
        keys = prpara.keys()  
        tmp = prpara[keys[0]]         
        params = json_parse(tmp.value)
        keys = params.keys()
        if max(keys eq 'CFGFILE') gt 0 then begin
          return, red_align_model(params['CFGFILE'])
        endif
      endif
      return, ''
    endif
  endif
  
  keys = prpara.keys()          ; Just a single key
  tmp = prpara[keys[0]]         
  params = json_parse(tmp.value)
  keys = params.keys()

  if max(keys eq 'ALIGN_MAP') gt 0 then return, 'projective'   ; Projective transform
  if max(keys eq 'ALIGN_MAP_X') gt 0 then return, 'polywarp'   ; Polywarp tranform
  return, ''                                                   ; No alignment
  
end


filenames = ['/scratch/mats/2024-06-05/CRISP-polywarp/momfbd_nopd/13:02:16/6173/cfg/momfbd_reduc_6173_00000.cfg' $
             , '/scratch/mats/2024-06-05/CRISP-polywarp/momfbd_nopd/13:02:16/6173/cfg/results/camXXXI_2024-06-05T13:02:16_00000_12.00ms_G00.00_6173_6173_+70_lc3.fitsheader' $
             , '/scratch/mats/2024-06-05/CRISP-polywarp/momfbd_nopd/13:02:16/6173/cfg/results/camXXXI_2024-06-05T13:02:16_00000_12.00ms_G00.00_6173_6173_+70_lc3.momfbd' $
             , '/scratch/mats/2016.09.19/CRISP-aftersummer/momfbd_nopd/09:28:36/6302/cfg/momfbd_reduc_6302_00000.cfg' $
             , '/scratch/mats/2016.09.19/CRISP-aftersummer/momfbd_nopd/09:28:36/6302/cfg/results/camXIX_2016-09-19T09:28:36_00000_6302_6302_-990_lc3.momfbd' $
             , 'dummy_file_name' $
            ]

for i = 0, n_elements(filenames)-1 do begin

  print, '-----'
  print, filenames[i]
  print, red_align_model(filenames[i])
  
endfor

end
