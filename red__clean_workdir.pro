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
      
      print, inam + ' : Specifying both prefilters and timestamps not implemented (yet).'
      return
      
    end

    ;; --------------------
    
    n_elements(prefilters) gt 0 : begin
      for ipref = 0, n_elements(prefilters)-1 do begin

        pref = prefilters[ipref]
        searchstring = [ 'flats/cam*_'+pref+'_????_[+-]*' $
                         , 'pinhs/cam*_'+pref+'_????_[+-]*.pinh.fits' $
                         , 'polcal_cubes/cam*_'+pref+'_polcalcube.fits' $
                         , 'gaintables/cam*_'+pref+'_????_[+-]*' $
                         , 'gaintables/??:??:??/cam*_?????_'+pref+'_????_[+-]*' $
                         , 'cmap_intdif/*/cam*.'+pref+'.intdif.*' $
                         , '*mfbd*/??:??:??/'+pref+'/cfg/results' $
                       ] 

        file_list = file_search(searchstring, count = Nfiles)

        if Nfiles gt 0 then begin

          openw, /get_lun, lun, 'clean_workdir.sh'
          printf, lun, '#!/bin/bash'
          for ii = 0, Nfiles-1 do printf, lun, 'rm '+file_list[ii]
          free_lun, lun
          spawn, 'chmod o+x clean_workdir.sh'
          
        endif else begin
          print, inam + ' : There are no files associated with prefilter ' + pref
          return
        endelse

      endfor                    ; ipref
    end
    
    ;; --------------------
    
    n_elements(timestamps) gt 0 : begin
      for itime = 0, n_elements(timestamps)-1 do begin

        timestamp = timestamps[itime]
        searchstring = ['gaintables/'+ timestamp $
                        , 'cmap_intdif/'+ timestamp $
                        , '*mfbd*/'+timestamp+'/????/cfg/results' $
                       ] 

        dir_list = file_search(searchstring, count = Nfiles)

        if Nfiles gt 0 then begin
          
          openw, /get_lun, lun, 'clean_workdir.sh'
          printf, lun, '#!/bin/bash'
          for ii = 0, Nfiles-1 do printf, lun, 'rm '+dir_list[ii]
          free_lun, lun
          spawn, 'chmod o+x clean_workdir.sh'
          
        endif else begin
          print, inam + ' : There are no files associated with timestamp ' + timestamp 
          return
        endelse
        
      endfor                    ; itime
    end
    
    ;; --------------------
    
    else : begin
      dir_list = ['darks' $
                  , 'flats' $
                  , 'data' $
                  , 'gaintables' $
                  , 'pinhs' $
                  , 'polcal_sums' $
                  , 'polcal_cubes' $
                 ]

      momfbd_output = file_search('*mfbd*/??:??:??/????/cfg/results', count = Nmomfbd)
      if Nmomfbd gt 0 then red_append, dir_list, momfbd_output

      Nfiles = n_elements(dir_list)
      
      if Nfiles gt 0 then begin

        openw, /get_lun, lun, 'clean_workdir.sh'
        printf, lun, '#!/bin/bash'
        for ii = 0, Nfiles-1 do printf, lun, 'rm '+dir_list[ii]
        free_lun, lun
        spawn, 'chmod o+x clean_workdir.sh'
        
      endif else return
    endcase
    
  endcase

  print
  print, inam + ' : Please inspect the new file "clean_workdir.sh"'
  print, inam + ' : Edit as needed and then execute it at your own risk.'
  print
  print, inam + ' : Take extra care not to delete your only copies of summed calibration data!'
  print
  
end

a = crispred(/dev)

a -> clean_workdir, prefilters = '8542'

stop

a -> clean_workdir, timestamps = '09:30:20'

stop

a -> clean_workdir

end
