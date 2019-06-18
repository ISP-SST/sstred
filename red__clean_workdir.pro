; docformat = 'rst'

;+
; Clean a work directory to save disk space while preserving files
; that are a record of the processing. Use with care!
;
; Files that are deleted include data in the following subdirecories:
; darks, flats, gaintables, pinhs, cmap_intdif, momfbd*/.../results,  
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
; 
; 
; 
; 
; 
; :Keywords:
; 
;   prefilters : in, optional, type=strarr
;
;     Limit the cleaning to the specified prefilters.
;
;   timestamps : in, optional, type=strarr
;   
;     Limit the cleaning to the specified timestamps. 
; 
; 
; :History:
; 
;   2019-06-17 : MGL. First version.
; 
;-
pro red::clean_workdir, timestamps = timestamps, prefilters = prefilters

  inam = red_subprogram(/low, calling = inam1)

  case 1 of

    n_elements(prefilters) gt 0 and n_elements(timestamps) gt 0 : begin
      
      print, inam + ' : Specifying both prefilters and timestamps not yet implemented.'
      return
      
    end

    ;; --------------------
    
    n_elements(prefilters) gt 0 : begin
      for ipref = 0, n_elements(prefilters)-1 do begin

        pref = prefilters[ipref]
        searchstring = [ 'flats/cam*_'+pref+'_????_[+-]*' $
                         , 'gaintables/cam*_'+pref+'_????_[+-]*' $
                         , 'gaintables/??:??:??/cam*_?????_'+pref+'_????_[+-]*' $
                         , 'cmap_intdif/*/cam*.'+pref+'.intdif.*' $
                         , '*mfbd*/??:??:??/'+pref $
                       ] 

        file_list = file_search(searchstring, count = Nfiles)

        if Nfiles gt 0 then begin
          print
          print, inam + ' : Will delete files associated with prefilter ' + pref + ':'
          print, file_list, format = '(a0)'
          print
          ;;file_delete, file_list, /allow_nonexistent, /quiet, /recursive
        endif else begin
          print, inam + ' : There are no files associated with prefilter ' + pref 
        endelse

      endfor                    ; ipref
    end
    
    ;; --------------------
    
    n_elements(timestamps) gt 0 : begin
      for itime = 0, n_elements(timestamps)-1 do begin

        timestamp = timestamps[itime]
        searchstring = ['gaintables/', 'cmap_intdif/', '*mfbd*/'] + timestamp

        dir_list = file_search(searchstring, count = Nfiles)

        if Nfiles gt 0 then begin
          print
          print, inam + ' : Will delete files associated with timestamp ' + timestamp + ':'
          print, dir_list, format = '(a0)'
          print
          ;;file_delete, dir_list, /allow_nonexistent, /quiet, /recursive
        endif else begin
          print, inam + ' : There are no files associated with timestamp ' + timestamp 
        endelse

      endfor                    ; itime
    end

    ;; --------------------
    
    else : begin
      dir_list = self.out_dir + ['darks' $
                                 , 'flats' $
                                 , 'data' $
                                 , 'gaintables' $
                                 , 'pinhs' $
                                 , 'polcal_sums' $
                                 , 'polcal_cubes' $
                                ]

      momfbd_output = file_search('*mfbd*/??:??:??/????/cfg/results', count = Nmomfbd)
      if Nmomfbd gt 0 then red_append, dir_list, momfbd_output

      print
      print, inam + ' : Will delete the following directories' 
      print, dir_list, format = '(a0)'
      print
      ;;file_delete, dir_list, /allow_nonexistent, /quiet, /recursive
    end
    
  endcase
  
end

a = crispred(/dev)

a -> clean_workdir, prefilters = '8542'

stop

a -> clean_workdir, timestamps = '09:30:20'

stop

a -> clean_workdir

end
