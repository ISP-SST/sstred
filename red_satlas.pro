; docformat = 'rst'

;+
;
; :history:
; 
;    2013-07-11 : MGL. Renamed from satlas. Modified mechanism for
;                 finding the data file within the path. Should now
;                 pick the version that is in the same directory as
;                 this file, that is, in the crispred repository.
; 
; 
; 
;
;-
pro red_satlas,xstart,xend,outx,outy,nm=nm,nograv=nograv,nocont=nocont

  ;; Find the input data
  this_dir = file_dirname( routine_filepath("red_satlas"), /mark )
  restore, this_dir+'ftsatlas.idlsave'

  if not keyword_set(nocont) then YL_FTS/=CINT_FTS
  if keyword_set(nm) then begin
     xstart=xstart/10.d0
     xend=xend/10.d0
  endif
;
  c=299792458.d0                ; light speed in m/s
;
  ;pos=where(XL_FTS gt xstart AND XL_FTS lt xend)
  ;outx=xl_fts;[pos]
  if not keyword_set(nograv) then xl_fts*=(1.d0-633.d0/c)
  pos=where(XL_FTS ge xstart AND XL_FTS le xend)
  outx=xl_fts[pos]
  outy=YL_FTS[pos]
  return
end
