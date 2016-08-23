; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :author:
; 
; 
; 
; 
; :returns:
; 
; 
; :Params:
; 
; 
; :Keywords:
; 
;    wb_states  : 
;   
;   
;   
;    outformat  : 
;   
;   
;   
;    numpoints  : 
;   
;   
;   
;    modes  : 
;   
;   
;   
;    date_obs  : 
;   
;   
;   
;    state  : 
;   
;   
;   
;    no_descatter  : 
;   
;   
;   
;    global_keywords  : 
;   
;   
;   
;    unpol  : 
;   
;   
;   
;    skip  : 
;   
;   
;   
;    pref  : 
;   
;   
;   
;    escan  : 
;   
;   
;   
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
;
;   2016-02-15 : MGL. Use red_loadbackscatter. Remove keyword descatter,
;                new keyword no_descatter.
;
;   2016-08-23 : THI. Rename camtag to detector and channel to camera,
;                so the names match those of the corresponding SolarNet
;                keywords.
; 
;-
pro red::prepmomfbd2, wb_states = wb_states, outformat = outformat, numpoints = numpoints, $
                      modes = modes, date_obs = date_obs, state = state, no_descatter = no_descatter, $
                      global_keywords = global_keywords, unpol = unpol, skip = skip, $
                      pref = pref, escan = escan
                                ;
  inam = 'red::prepmomfbd2 : '
                                ;
                                ; get keywords
  if(~keyword_set(date_obs)) then begin
     date_obs = ' '
     read, date_obs, prompt = inam+'type date_obs (YYYY-MM-DD): '
  endif
  if(~keyword_set(numpoints)) then numpoints = '88'
  if(~keyword_set(outformat)) then outformat = 'MOMFBD'
  if(~keyword_set(modes)) then modes = '2-29,31-36,44,45'
  
  if ~ptr_valid(self.data_dirs) then begin
     print, inam+' : ERROR : undefined data_dir'
     return
  endif

  Ndirs = n_elements(*self.data_dirs)
                                ;
                                ; get states from the data folder
                                ;
  for fff = 0, Ndirs - 1 do begin
     data_dir = self.data_dirs[fff]
     spawn, 'find ' + (data_dir) + '/' + self.camt + '/ | grep im.ex', files
     folder_tag = strsplit(data_dir,'/',/extract)
     nn = n_elements(folder_tag) - 1
     folder_tag = folder_tag[nn]

     nf = n_elements(files)
     if(files[0] eq '') then begin
        print, inam + 'ERROR -> no frames found in '+(data_dir)
        return
     endif
                                ;
     files = red_sortfiles(temporary(files))
                                ;
                                ; get image unique states
     stat = red_getstates(files)
                                ;red_flagtuning, stat
     states = stat.hscan+'.'+stat.state
     pos = uniq(states, sort(states))
     ustat = stat.state[pos]
     ustatp = stat.pref[pos]
                                ;ustats = stat.scan[pos]
     ntt = n_elements(ustat)
     hscans = stat.hscan[pos]
                                ;
                                ; Get unique prefilters
     upref = stat.pref[uniq(stat.pref, sort(stat.pref))]
     np = n_elements(upref)
                                ;
                                ; Get scan numbers
     uscan = stat.rscan[uniq(stat.rscan, sort(stat.rscan))]
     ns = n_elements(uscan)
                                ;
                                ; Create a reduc file per prefilter and scan number?
     outdir = self.out_dir + '/momfbd/'+folder_tag
     file_mkdir, outdir
                                ;
                                ; self -> getdetectors, dir = self.data_dir
     self.camttag = 'camXXV'
     self.camrtag = 'camXVIII'
     self.camwbtag = 'camXIX'
                                ;
                                ; Print cams
                                ;
     print, inam + 'cameras found:'
     print, ' WB   -> '+self.camwbtag
     print, ' NB_T -> '+self.camttag
     print, ' NB_R -> '+self.camrtag

                                ;
                                ; Extensions
                                ;
     case outformat of
        'ANA': exten = '.f0'
        'MOMFBD': exten = '.momfbd'
        ELSE: begin
           print, inam+'WARNING -> could not determine a file type for the output'
           exten = ''
        end
        
     endcase
                                ;
                                ; choose offset state
                                ;
     for ss = 0L, ns - 1 do begin
        if(keyword_set(escan)) then if(ss ne escan) then continue
        
        scan = uscan[ss]
                                ;
        for pp = 0L, np - 1 do begin
           if(keyword_set(pref)) then begin
              if(upref[pp] NE pref) then begin
                 print, inam + 'Skipping prefilter -> ' + upref[pp]
                 continue
              endif
           endif
                                ;
                                ;
                                ; load align clips
                                ;

           clipfile = self.out_dir + '/calib/align_clips.'+upref[pp]+'.sav'
           IF(~file_test(clipfile)) THEN BEGIN
              print, inam + 'ERROR -> align_clip file not found'
              print, inam + '-> you must run red::getalignclips first!'
              return
           endif
           restore, clipfile
           wclip = acl[0]
           tclip = acl[1]
           rclip = acl[2]

           lam = strmid(string(float(upref[pp]) * 1.e-10), 2)
                                ;
           cfg_file = 'momfbd.reduc.'+upref[pp]+'.'+scan+'.cfg'
           outdir = self.out_dir + '/momfbd/'+upref[pp]+'/cfg/'
           file_mkdir, outdir
           rdir = self.out_dir + '/momfbd/'+upref[pp]+'/cfg/results/'
           file_mkdir, rdir
           openw, lun, outdir + cfg_file, /get_lun
                                ;
                                ; Image numbers
                                ;
           numpos = where((stat.rscan eq uscan[ss]) AND (stat.star eq 0B) AND (stat.pref eq upref[pp]), ncount)
           if(ncount eq 0) then continue
           n0 = stat.nums[numpos[0]]
           n1 = stat.nums[numpos[ncount-1]]
           print, inam+'Prefilter = '+upref[pp]+' -> scan = '+uscan[ss]+' -> image range = ['+n0+' - '+n1+']'
                                ;
                                ; WB anchor channel
                                ;
           printf, lun, 'object{'
           printf, lun, '  WAVELENGTH=' + lam
           printf, lun, '  OUTPUT_FILE=results/'+self.camwbtag+'.'+scan+'.'+upref[pp]
           printf, lun, '  channel{'
           printf, lun, '    IMAGE_DATA_DIR='+self.out_dir+'/data/'+folder_tag+ '/' +self.camwbtag+'_nostate/'
           printf, lun, '    FILENAME_TEMPLATE='+self.camwbtag+'.'+scan+'.'+upref[pp]+'.%07d'
           printf, lun, '    DIVERSITY=0.0 mm'
           printf, lun, '    GAIN_FILE=' + file_search(self.out_dir+'gaintables/'+self.camwbtag + $
                                                       '.' + upref[pp]+'*.gain')
           printf, lun, '    DARK_TEMPLATE='+self.out_dir+'darks/'+self.camwbtag+'.summed.0000001'
           printf, lun, '    DARK_NUM=0000001'
           printf, lun, '    ' + wclip
           if upref[pp] EQ '8542' AND ~keyword_set(no_descatter) then begin
              self -> loadbackscatter, self.camwbtag, upref[pp], bg, psf, bgfile = bgf, bpfile = psff
              printf, lun, '    PSF='+psff
              printf, lun, '    BACK_GAIN='+bgf
           endif 
                                ;
           xofile = (file_search(self.out_dir+'/calib/'+self.camwbtag+'.*.xoffs'))[0]
           yofile = (file_search(self.out_dir+'/calib/'+self.camwbtag+'.*.yoffs'))[0]
           
           if(file_test(xofile)) then printf, lun, '    XOFFSET='+xofile
           if(file_test(yofile)) then printf, lun, '    YOFFSET='+yofile
                                ;
           printf, lun, '    INCOMPLETE'
           printf, lun, '  }'
           printf, lun, '}'
                                ;
                                ; Loop all wavelengths
                                ;
           pos1 = where((ustatp eq upref[pp]), count)
           if(count eq 0) then continue
           ustat1 = ustat[pos1]
                                ;
           for ii = 0L, count - 1 do BEGIN
              
                                ;
                                ; external states?
                                ;
              if(keyword_set(state)) then begin
                 dum = where(state eq ustat1[ii], cstate)
                 if(cstate eq 0) then continue
                 print, inam+'found '+state+' -> scan = '+uscan[ss]
              endif
                                ;
              self -> whichoffset, ustat1[ii], xoff = xoff, yoff = yoff
                                ;
                                ; Trans. camera
                                ;
              istate = red_encode_scan(hscans[pos1[ii]], scan)+'.'+ustat1[ii]
              
                                ;
                                ; lc4?
                                ;
              tmp = strsplit(istate,'.', /extract)
              ntmp = n_elements(tmp)
                                ;
              idx = strsplit(ustat1[ii],'.')
              nidx = n_elements(idx)
              iwavt = strmid(ustat1[ii], idx[0], idx[nidx-1]-1)
                                ;
              if(keyword_set(skip)) then begin
                 dum = where(iwavt eq skip, ccout)
                 if ccout ne 0 then begin
                    print, inam+'skipping state -> '+ustat1[ii]
                    continue
                 endif
              endif
                                ;
                                ;
              printf, lun, 'object{'
              printf, lun, '  WAVELENGTH=' + lam
              printf, lun, '  OUTPUT_FILE=results/'+self.camttag+'.'+istate 
              printf, lun, '  channel{'
              printf, lun, '    IMAGE_DATA_DIR='+self.out_dir+'/data/'+folder_tag+ '/' +self.camttag+'/'
              printf, lun, '    FILENAME_TEMPLATE='+self.camttag+'.'+istate+'.%07d'
              printf, lun, '    DIVERSITY=0.0 mm'
              if(~keyword_set(unpol)) then begin
                                ;search = self.out_dir+'/gaintables/'+self.camttag + '.' + ustat1[ii] + '*.gain'
                 search = self.out_dir+'/gaintables/'+self.camttag + '.' + istate+'.gain'
              endif Else begin
                                ;

                 iistate = strjoin((strsplit(istate,'.',/extract))[0:2], '.')
                 
                 search = self.out_dir+'/gaintables/'+self.camttag + $
                          '.' + iistate+ '*unpol.gain'
                                ;if tmp[ntmp-1] eq 'lc4' then search = self.out_dir+'/gaintables/'+$
                                ;                                      self.camttag + '.' + ustat[pos[ii]] + $
                                ;                                      '*.gain'
                                ;
              endelse
              
              printf, lun, '    GAIN_FILE=' + file_search(search)
              printf, lun, '    DARK_TEMPLATE='+self.out_dir+'/darks/'+self.camttag+'.summed.0000001'
              printf, lun, '    DARK_NUM=0000001'
              printf, lun, '    ' + tclip
                                ;
              xofile = self.out_dir+'/calib/'+self.camttag+'.'+xoff
              yofile = self.out_dir+'/calib/'+self.camttag+'.'+yoff
              if(file_test(xofile)) then printf, lun, '    XOFFSET='+xofile
              if(file_test(yofile)) then printf, lun, '    YOFFSET='+yofile
                                ;
              if upref[pp] EQ '8542' AND ~keyword_set(no_descatter) then begin
                 self -> loadbackscatter, self.camttag, upref[pp], bg, psf, bgfile = bgf, bpfile = psff
                 printf, lun, '    PSF='+psff
                 printf, lun, '    BACK_GAIN='+bgf
              endif 
              printf, lun, '    INCOMPLETE'
              printf, lun, '  }'
              printf, lun, '}'  
                                ;
                                ; Reflec. camera
                                ;
              printf, lun, 'object{'
              printf, lun, '  WAVELENGTH=' + lam
              printf, lun, '  OUTPUT_FILE=results/'+self.camrtag+'.'+istate 
              printf, lun, '  channel{'
              printf, lun, '    IMAGE_DATA_DIR='+self.out_dir+'/data/'+folder_tag+ '/' +self.camrtag+'/'
              printf, lun, '    FILENAME_TEMPLATE='+self.camrtag+'.'+istate+'.%07d'
              printf, lun, '    DIVERSITY=0.0 mm' 
              if(~keyword_set(unpol)) then begin
;              search = self.out_dir+'/gaintables/'+self.camrtag + '.' + ustat1[ii] + '*.gain'
                 search = self.out_dir+'/gaintables/'+self.camrtag + '.' + istate+'.gain'

              endif Else begin
                 
                                ;             search = self.out_dir+'/gaintables/'+self.camrtag + '.' + istate+'.gain'

                 iistate = strjoin((strsplit(istate,'.',/extract))[0:2], '.')
                 
                 search = self.out_dir+'/gaintables/'+self.camrtag + $
                          '.' + iistate+ '*unpol.gain'

                                ;if tmp[ntmp-1] eq 'lc4' then search = self.out_dir+'/gaintables/'+$
                                ;                                      self.camrtag + '.' + ustat[pos[ii]] + $
                                ;                                      '*.gain'
              endelse
              printf, lun, '    GAIN_FILE=' + file_search(search)
              printf, lun, '    DARK_TEMPLATE='+self.out_dir+'/darks/'+self.camrtag+'.summed.0000001'
              printf, lun, '    DARK_NUM=0000001'
              printf, lun, '    ' + rclip
              xofile = self.out_dir+'/calib/'+self.camrtag+'.'+xoff
              yofile = self.out_dir+'/calib/'+self.camrtag+'.'+yoff
              if(file_test(xofile)) then printf, lun, '    XOFFSET='+xofile
              if(file_test(yofile)) then printf, lun, '    YOFFSET='+yofile
                                ;
              if upref[pp] EQ '8542' AND ~keyword_set(no_descatter) then begin
                 self -> loadbackscatter, self.camrtag, upref[pp], bg, psf, bgfile = bgf, bpfile = psff
                 printf, lun, '    PSF='+psff
                 printf, lun, '    BACK_GAIN='+bgf
              endif 
              printf, lun, '    INCOMPLETE'
              printf, lun, '  }'
              printf, lun, '}'
                                ;
                                ; WB with states (for de-warping to the anchor, 
                                ; only to remove rubbersheet when differential seeing is strong)
                                ;
              if(keyword_set(wb_states)) then begin
                 printf, lun, 'object{'
                 printf, lun, '  WAVELENGTH=' + lam
                 printf, lun, '  WEIGHT=0.00'
                 printf, lun, '  OUTPUT_FILE=results/'+self.camwbtag+'.'+istate 
                 printf, lun, '  channel{'
                 printf, lun, '    IMAGE_DATA_DIR='+self.out_dir+'/data/'+folder_tag+ '/' +self.camwbtag+'/'
                 printf, lun, '    FILENAME_TEMPLATE='+self.camwbtag+'.'+istate+'.%07d'
                 printf, lun, '    DIVERSITY=0.0 mm'
                 printf, lun, '    GAIN_FILE=' + file_search(self.out_dir+'/gaintables/'+self.camwbtag + $
                                                             '.' + upref[pp] + '*.gain')
                 printf, lun, '    DARK_TEMPLATE='+self.out_dir+'/darks/'+self.camwbtag+'.summed.0000001'
                 printf, lun, '    DARK_NUM=0000001'
                 printf, lun, '    ' + wclip
                 if upref[pp] EQ '8542' AND ~keyword_set(no_descatter) then begin
                    self -> loadbackscatter, self.camwbtag, upref[pp], bg, psf, bgfile = bgf, bpfile = psff
                    printf, lun, '    PSF='+psff
                    printf, lun, '    BACK_GAIN='+bgf
                 endif 
                 xofile = self.out_dir+'/calib/'+self.camwbtag+'.'+xoff
                 yofile = self.out_dir+'/calib/'+self.camwbtag+'.'+yoff
                 if(file_test(xofile)) then printf, lun, '    XOFFSET='+xofile
                 if(file_test(yofile)) then printf, lun, '    YOFFSET='+yofile

                 printf, lun, '    INCOMPLETE'
                 printf, lun, '  }'
                 printf, lun, '}'
              endif          
                                ;
           endfor
                                ;
                                ; Global keywords
                                ;
           printf, lun, 'PROG_DATA_DIR=./data/'
           printf, lun, 'DATE_OBS='+date_obs
           printf, lun, 'IMAGE_NUMS='+n0+'-'+n1
           printf, lun, 'BASIS=Karhunen-Loeve'
           printf, lun, 'MODES='+modes
           printf, lun, 'NUM_POINTS='+numpoints
           printf, lun, 'TELESCOPE_D=0.97'
           printf, lun, 'ARCSECPERPIX='+self.image_scale
           printf, lun, 'PIXELSIZE=16.0E-6'
           printf, lun, 'GETSTEP=getstep_conjugate_gradient'
           printf, lun, 'GRADIENT=gradient_diff'
           printf, lun, 'MAX_LOCAL_SHIFT=30'
           printf, lun, 'NEW_CONSTRAINTS'
           printf, lun, 'FILE_TYPE='+outformat
           printf, lun, 'FAST_QR'
           if(outformat eq 'MOMFBD') then printf, lun, 'GET_PSF'
           if(outformat eq 'MOMFBD') then printf, lun, 'GET_PSF_AVG'
                                ;
                                ; External keywords?
                                ;
           if(keyword_set(global_keywords)) then begin
              nk = n_elements(global_keywords)
              for ki = 0L, nk -1 do printf, lun, global_keywords[ki]
           endif
                                ;
           free_lun, lun
        endfor
     endfor
  endfor
                                ;
  print, inam+'done!'
  return
end
