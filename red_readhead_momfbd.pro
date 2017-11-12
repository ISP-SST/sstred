; docformat = 'rst'

;+
; Return a fits-type header from a MOMFBD format file.
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
;    version : out, optional, type=string
;
;       The version of the momfbd (or redux) program that wrote the
;       file.
; 
; :History:
; 
;    2017-03-10 : MGL. Moved reading of headers from MOMFBD format
;                 files from red_readhead.pro.
; 
;    2017-03-17 : MGL. New keyword "version".
;
;    2017-11-10 : MGL. Use the new roi field in the struct returned by
;                 momfbd_read to get predictable dimensions that match
;                 the output of red_readdata.
;
;    2017-11-12 : MGL. Update dimensions to match new convention in
;                 struct returned by momfbd_read.
; 
; 
; 
;-
function red_readhead_momfbd, fname, version = version

  compile_opt idl2

  mr = momfbd_read(fname, /names) ; Use /names to avoid reading the data parts

  red_mkhdr, header, 4, [mr.roi[1]-mr.roi[0]+1 $
                         , mr.roi[3]-mr.roi[2]+1] - 2*mr.margin 

  header = header[where(header ne replicate(' ',80))] ; Remove blank lines

  date_avg = mr.date + 'T' + mr.time
  fxaddpar, header, 'DATE-AVG', date_avg, ' ', before='COMMENT'

  ;; Set version to the version of the momfbd (or redux) program that
  ;; wrote the file.
  version = mr.version

  return, header
        
end


head = red_readhead_momfbd('/scratch/mats/2016.09.19/CHROMIS/momfbd_test_pd2/09:28:36/3950/cfg/results/camXXX_2016-09-19T09:28:36_00071_12.00ms_G10.00_3934_3934_-313.momfbd', version = version)

print, version

end
