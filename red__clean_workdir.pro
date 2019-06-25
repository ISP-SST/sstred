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
      searchstring = [ 'flats/cam*_'+prefilters+'_????_[+-]*' $
                       , 'pinhs/cam*_'+prefilters+'_????_[+-]*.pinh.fits' $
                       , 'polcal_cubes/cam*_'+prefilters+'_polcalcube.fits' $
                       , 'polcal_sums/Crisp-?/cam*_lp*_qw*_'+prefilters+'_lc?.pols.*' $
                       , 'gaintables/cam*_'+prefilters+'_????_[+-]*' $
                       , 'gaintables/??:??:??/cam*_?????_'+prefilters+'_????_[+-]*' $
                       , 'cmap_intdif/*/cam*.'+prefilters+'.intdif.*' $
                     ] 
      file_list = file_search(searchstring, count = Nfiles)

      searchstring = ['*mfbd*/??:??:??/'+prefilters+'/cfg/results'] 
      dir_list = file_search(searchstring, count = Ndirs)

      
      if Nfiles+Ndirs gt 0 then begin
        openw, /get_lun, lun, 'clean_workdir.sh'
        printf, lun, '#!/bin/bash'
        for ii = 0, Nfiles-1 do printf, lun, 'rm '+file_list[ii]
        for ii = 0, Ndirs-1 do printf, lun, 'rm -r '+dir_list[ii]
        free_lun, lun
      endif else begin
        print, inam + ' : There are no files or directories associated with prefilter(s) ' + prefilters
        return
      endelse

    end
    
    ;; --------------------
    
    n_elements(timestamps) gt 0 : begin
      searchstring = ['gaintables/'+ timestamps $
                      , 'cmap_intdif/'+ timestamps $
                      , '*mfbd*/'+timestamps+'/????/cfg/results' $
                     ] 

      dir_list = file_search(searchstring, count = Nfiles)

      if Nfiles gt 0 then begin
        openw, /get_lun, lun, 'clean_workdir.sh'
        printf, lun, '#!/bin/bash'
        for ii = 0, Nfiles-1 do printf, lun, 'rm -r '+dir_list[ii]
        free_lun, lun
      endif else begin
        print, inam + ' : There are no directories associated with timestamp(s) ' + timestamps
        return
      endelse
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
        for ii = 0, Nfiles-1 do printf, lun, 'rm -r '+dir_list[ii]
        free_lun, lun
        
      endif else return
    endcase
    
  endcase

  ;; Make it executable
  spawn, 'chmod a+x clean_workdir.sh'

  print
  print, inam + ' : Please inspect the new file "clean_workdir.sh"'
  print, inam + ' : Edit as needed and then execute it at your own risk.'
  print
  print, inam + ' : Take extra care not to delete your only copies of summed calibration data!'
  print, inam + ' : If you have soft links to directories in another workdir, and not just to their contents, then you do not want to rm -r the links!'
  print
  
end

; Todo: Examine all found files with file_info() and use this info to
; warn about files without the needed permissions and to do the right
; thing with symbolic links. And maybe other things.

a = crispred(/dev)

a -> clean_workdir, prefilters = '8542'
a -> clean_workdir, prefilters = ['6302', '8542']

stop

a -> clean_workdir, timestamps = '09:30:20'

stop

a -> clean_workdir

end
