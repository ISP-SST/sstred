PRO phv_update_img, wid, img, xsz, ysz, xpos, ypos

  COMPILE_OPT hidden
  
  WSET, wid
  dispim = congrid( img, xsz, ysz, cubic = -0.5 )
  TVSCL, red_histo_opt(dispim), xpos, ypos
  
END

PRO phv_timer_callback, id, info

  COMPILE_OPT hidden

  IF (*info).exit THEN return
 
  IF (*info).disp_mapped THEN BEGIN
    phv_update_img, (*info).draw, (*info).ph2, (*info).drawSz[0]/2, (*info).drawSz[1], (*info).drawSz[0]/2, 0
    XYOUTS, (*info).drawSz[0]/2, (*info).drawSz[1]-20, 'ph2', CHARSIZE=1, /DEVICE
    (*info).disp_mapped = 0
  ENDIF ELSE BEGIN
    phv_update_img, (*info).draw, rdx_img_project( (*info).map, (*info).ph1, /pres ), (*info).drawSz[0]/2, (*info).drawSz[1], (*info).drawSz[0]/2, 0
    XYOUTS, (*info).drawSz[0]/2, (*info).drawSz[1]-20, 'map(ph1)', CHARSIZE=1, /DEVICE
    (*info).disp_mapped = 1
  ENDELSE
  (*info).timerid = Timer.set( 1, 'phv_timer_callback', info )
END

PRO phv_update_imgs, info

  COMPILE_OPT hidden
  
  res = Timer.Cancel( (*info).timerid )
  
  phv_update_img, (*info).draw, (*info).ph1, (*info).drawSz[0]/2, (*info).drawSz[1], 0, 0
  XYOUTS, 5, (*info).drawSz[1]-20, 'ph1', CHARSIZE=1, /DEVICE
  phv_update_img, (*info).draw, (*info).ph2, (*info).drawSz[0]/2, (*info).drawSz[1], (*info).drawSz[0]/2, 0
  XYOUTS, (*info).drawSz[0]/2, (*info).drawSz[1]-20, 'ph2', CHARSIZE=1, /DEVICE

 (*info).timerid = Timer.set( 1, 'phv_timer_callback', info )

END

PRO red_phverify_event, ev

  COMPILE_OPT hidden
  WIDGET_CONTROL, ev.top, GET_UVALUE=info
  
  IF (ev.ID eq (*info).m_flipx) OR (ev.ID eq (*info).m_flipy) THEN BEGIN
    IF (ev.ID eq (*info).m_flipx) THEN mm = [ [ -1, 0, (*info).sz2[0] ], [0, 1, 0], [0, 0, 1] ] $
    ELSE mm = [ [ 1, 0, 0], [ 0, -1, (*info).sz2[1] ], [0, 0, 1] ]
    (*info).map = (*info).map # mm
    WIDGET_CONTROL, (*info).table, SET_VALUE=(*info).map
  ENDIF
  
  IF (ev.ID eq (*info).m_calib) THEN BEGIN
    h_init = (*info).map
    nref = 10
    (*info).map = rdx_img_align( (*info).ph1, (*info).ph2, nref=(*info).nref, max_shift=(*info).max_shift, $
      h_init=h_init, verbose=0 )
    WIDGET_CONTROL, (*info).table, SET_VALUE=(*info).map
  ENDIF
  
  IF (ev.ID eq (*info).m_restore) THEN BEGIN
    (*info).map = (*info).omap
    WIDGET_CONTROL, (*info).table, SET_VALUE=(*info).map
  ENDIF

  IF (ev.ID eq (*info).m_unity) THEN BEGIN
    (*info).map = [ [ 1, 0, 0], [0, 1, 0], [0, 0, 1] ]
    WIDGET_CONTROL, (*info).table, SET_VALUE=(*info).map
  ENDIF

  IF (ev.ID eq (*info).b_done) THEN BEGIN
    (*info).exit = 1
    res = Timer.Cancel( (*info).timerid )
    wait,1
    WIDGET_CONTROL, ev.TOP, /DESTROY
  ENDIF

END

FUNCTION red_phverify, ph1, ph2, map, max_shift=max_shift, nref=nref, threshold=threshold

  sz1 = size( ph1, /dim )
  sz2 = size( ph2, /dim )

  if( (n_elements(sz1) ne 2) OR (n_elements(sz2) ne 2) ) then return, map

  if(n_elements(threshold) eq 0) then threshold = 0.0
  if(n_elements(nref) eq 0) then nref = 10
  if(n_elements(max_shift) eq 0) then max_shift = 200
  
  drawSz = max([[sz1],[sz2]],dim=2) / [1,2]

  base = WIDGET_BASE( /COLUMN )
  draw = WIDGET_DRAW( base, XSIZE = drawSz[0], YSIZE = drawSz[1] ) 
  control = WIDGET_BASE( base, /ROW )
  table = WIDGET_TABLE( control, VALUE=map, /NO_HEADERS, /EDITABLE, XSIZE=3, YSIZE=3, FRAME=3, $
    SCR_XSIZE=245, SCR_YSIZE=64, COLUMN_WIDTHS=[80,80,80], ROW_HEIGHTS=[20,20,20] )
  mat_buttons = WIDGET_BASE( control, COLUMN=2 )
  m_flipx = WIDGET_BUTTON( mat_buttons, value='Flip X' )
  m_flipy = WIDGET_BUTTON( mat_buttons, value='Flip Y' )
  m_calib = WIDGET_BUTTON( mat_buttons, value='Calibrate' )
  m_restore = WIDGET_BUTTON( mat_buttons, value='Restore' )
  m_unity = WIDGET_BUTTON( mat_buttons, value='Unity' )
  buttons = WIDGET_BASE( control, /COLUMN )
  b_done = WIDGET_BUTTON( buttons, value='Done' )

  WIDGET_CONTROL, base, /REALIZE
  WIDGET_CONTROL, draw, GET_VALUE = draw_ind
  
  info = { base:base, draw:draw_ind, ph1:ph1, ph2:ph2, map:map, omap:map, max_shift:max_shift, nref:nref, threshold:threshold, $
    sz1:sz1, sz2:sz2, drawSz:drawSz, table:table, $
    m_flipx:m_flipx, m_flipy:m_flipy, m_calib:m_calib, m_restore:m_restore, m_unity:m_unity, b_done:b_done, $
    timerid:0L, disp_mapped:0, exit:0 }
    
  info = PTR_NEW( info, /NO_COPY )
  WIDGET_CONTROL, base, SET_UVALUE=info
  
  phv_update_imgs, info

  XMANAGER, 'red_phverify', base
  
  return, (*info).map
  
END

PRO red_phverify_test
  cdir = '/nadir-scratch/tomas/alex/2017-10-19/CHROMIS/calib/'
  cfile = cdir + 'alignments.sav'
  restore, cfile
  nAl = n_elements(alignments)
  last_pref = ''
  for i=0, nAl-1 do begin
    if alignments[i].state2.prefilter ne last_pref then begin
      print, 'Verifying map: ',alignments[i].state1.fullstate,' -> ',alignments[i].state2.fullstate
      ph1 = red_readdata( alignments[i].state1.filename )
      ph2 = red_readdata( alignments[i].state2.filename )
      alignments[i].map = red_phverify( ph1, ph2, alignments[i].map, nref=5 )
      last_pref = alignments[i].state2.prefilter
    endif
  endfor
  
END
