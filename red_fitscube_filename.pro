; docformat = 'rst'

;+
; Construct a file name for a fitscube.
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
;    A fitscube filename constructed by use of the input data. 
; 
; :Params:
; 
;    filetype : in, type=string
; 
;       One of "wb" and "nb".
; 
;    prefilter : in, type=string
; 
;       The four-digit string identifying the prefilter.
; 
;    timestamps : in, type=strarr
; 
;       The timestamp directories contributing data to the fitscube. 
; 
;    scannos : in, type=strarr
; 
;       The scannumbers for each timestamp directory represented as a
;       comma- and dash-delimited string. Should have the same number
;       of elements as the timestamps parameter.
;
;    point_id : in, type=string
; 
;       A string identifying the observation. Typically the DATE-OBS
;       or STARTOBS keyword of data from the first timestamp directory
;       with the same pointing.
; 
; 
; :Keywords:
; 
;   datatags : in, optional, type=strarr
; 
;      Tags to be added to the filename. Might be used to tell
;      different versions apart.
;   
;   oldname : in, optional, type=boolean
;
;      For data from a single datestamp directory, construct the
;      filename as before multiple directories were implemented.
; 
; 
; :History:
; 
;   2021-11-28 : MGL. First version.
; 
;-
function red_fitscube_filename, filetype, prefilter, timestamps, scannos, point_id $ 
                                , datatags = datatags $
                                , oldname = oldname

  Nstamps = n_elements(timestamps)
  if Nstamps eq 0 then stop
  if n_elements(scannos) ne Nstamps then stop
  
  if keyword_set(oldname) and Nstamps eq 1 then begin
 
    midparts = timestamps[0] + '_scans=' + scannos[0]

  endif else begin

    for istamp = 0, Nstamps-1 do red_append, midparts, timestamps[istamp] + '=' + scannos[istamp]

  endelse

  if n_elements(datatags) eq 0 then begin
    filename = strjoin([filetype, prefilter, point_id, midparts], '_') + '.fits'    
  endif else begin
    filename = strjoin([filetype, prefilter, point_id, midparts, datatags], '_') + '.fits'    
  endelse
  
  return, filename
  

end

