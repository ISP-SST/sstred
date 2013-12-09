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
;    descatter  : 
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
;    div  : 
;   
;   
;   
;    nremove : 
;   
;   
;   
;    newgains :
; 
;
;
;
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
;
;   2013-06-13, JdlCR : added support for scan-dependent gains -> 
;                       using keyword "/newgains".
;
;   2013-06-28, JdlCR : added NF (object) option 
; 
;   2013-08-27 : MGL. Added support for logging. Let the subprogram
;                find out its own name.
;
;-
pro red::prepmomfbd, wb_states = wb_states, outformat = outformat, numpoints = numpoints, $
                     modes = modes, date_obs = date_obs, state = state, descatter = descatter, $
                     global_keywords = global_keywords, unpol = unpol, skip = skip, $
                     pref = pref, escan = escan, div = div, nremove=nremove, $
                     newgains=newgains, nf = nfac, weight = weight

  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  ;; Logging
  help, /obj, self, output = selfinfo 
  red_writelog, selfinfo = selfinfo

  ;; Get keywords
  if(~keyword_set(date_obs)) then begin
     date_obs = ' '
     read, date_obs, prompt = inam+' : type date_obs (YYYY-MM-DD): '
  endif
  if(~keyword_set(numpoints)) then numpoints = '88'
  if(~keyword_set(outformat)) then outformat = 'MOMFBD'
  if(~keyword_set(modes)) then modes = '2-29,31-36,44,45'
  if(n_elements(nremove) eq 0) then nremove=1
  if(n_elements(nfac) gt 0) then begin
     if(n_elements(nfac) eq 1) then nfac = replicate(nfac,3)
  endif
  
  
  ;; Get states from the data folder
  for fff = 0, self.ndir - 1 do begin
     data_dir = self.data_list[fff]
     spawn, 'find ' + data_dir + '/' + self.camt + '/ | grep im.ex | grep -v ".lcd."', files
     folder_tag = strsplit(data_dir,'/',/extract)
     nn = n_elements(folder_tag) - 1
     folder_tag = folder_tag[nn]

     nf = n_elements(files)
     if(files[0] eq '') then begin
        print, inam + ' : ERROR -> no frames found in '+data_dir
        return
     endif

     files = red_sortfiles(temporary(files))

     ;; Get image unique states
     stat = red_getstates(files)
     red_flagtuning, stat, nremove
     states = stat.hscan+'.'+stat.state
     pos = uniq(states, sort(states))
     ustat = stat.state[pos]
     ustatp = stat.pref[pos]
                                ;ustats = stat.scan[pos]

     ntt = n_elements(ustat)
     hscans = stat.hscan[pos]

     ;; Get unique prefilters
     upref = stat.pref[uniq(stat.pref, sort(stat.pref))]
     np = n_elements(upref)

     ;; Get scan numbers
     uscan = stat.rscan[uniq(stat.rscan, sort(stat.rscan))]
     ns = n_elements(uscan)

     ;; Create a reduc file per prefilter and scan number?
     outdir = self.out_dir + '/momfbd/' + folder_tag
     file_mkdir, outdir

     self -> getcamtags, dir = data_dir

     ;; Print cams
     print, inam + ' : cameras found:'
     print, ' WB   -> '+self.camwbtag
     print, ' NB_T -> '+self.camttag
     print, ' NB_R -> '+self.camrtag

     ;; Extensions
     case outformat of
        'ANA': exten = '.f0'
        'MOMFBD': exten = '.momfbd'
        ELSE: begin
           print, inam+' : WARNING -> could not determine a file type for the output'
           exten = ''
        end
     endcase

     ;; Choose offset state
     for ss = 0L, ns - 1 do begin
        if(keyword_set(escan)) then if(ss ne escan) then continue
        
        scan = uscan[ss]

        for pp = 0L, np - 1 do begin
           if(keyword_set(pref)) then begin
              if(upref[pp] NE pref) then begin
                 print, inam + ' : Skipping prefilter -> ' + upref[pp]
                 continue
              endif
           endif

           ;; Load align clips
           clipfile = self.out_dir + '/calib/align_clips.'+upref[pp]+'.sav'
           IF(~file_test(clipfile)) THEN BEGIN
              print, inam + ' : ERROR -> align_clip file not found'
              print, inam + ' : -> you must run red::getalignclps first!'
              continue
           endif
           restore, clipfile
           wclip = acl[0]
           tclip = acl[1]
           rclip = acl[2]

           lam = strmid(string(float(upref[pp]) * 1.e-10), 2)

           cfg_file = 'momfbd.reduc.'+upref[pp]+'.'+scan+'.cfg'
           outdir = self.out_dir + '/momfbd/'+folder_tag+'/'+upref[pp]+'/cfg/'
           file_mkdir, outdir
           rdir = self.out_dir + '/momfbd/'+folder_tag+'/'+upref[pp]+'/cfg/results/'
           file_mkdir, rdir
           ddir = self.out_dir + '/momfbd/'+folder_tag+'/'+upref[pp]+'/cfg/data/'
           file_mkdir, ddir
           if(n_elements(lun) gt 0) then free_lun, lun
           openw, lun, outdir + cfg_file, /get_lun, width=2500

           ;; Image numbers
           numpos = where((stat.rscan eq uscan[ss]) AND (stat.star eq 0B) AND (stat.pref eq upref[pp]), ncount)
           if(ncount eq 0) then continue
           n0 = stat.nums[numpos[0]]
           n1 = stat.nums[numpos[ncount-1]]
           nall = strjoin(stat.nums[numpos],',')
           print, inam+' : Prefilter = '+upref[pp]+' -> scan = '+uscan[ss]+' -> image range = ['+n0+' - '+n1+']'

           ;; WB anchor channel
           printf, lun, 'object{'
           printf, lun, '  WAVELENGTH=' + lam
           printf, lun, '  OUTPUT_FILE=results/'+self.camwbtag+'.'+scan+'.'+upref[pp]
           if(n_elements(weight) eq 3) then printf, lun, '  WEIGHT='+string(weight[0])
           printf, lun, '  channel{'
           printf, lun, '    IMAGE_DATA_DIR='+self.out_dir+'/data/'+folder_tag+ '/' +self.camwbtag+'_nostate/'
           printf, lun, '    FILENAME_TEMPLATE='+self.camwbtag+'.'+scan+'.'+upref[pp]+'.%07d'
                                ; printf, lun, '    DIVERSITY=0.0 mm'
           printf, lun, '    GAIN_FILE=' + file_search(self.out_dir+'gaintables/'+self.camwbtag + $
                                                       '.' + upref[pp]+'*.gain')
           printf, lun, '    DARK_TEMPLATE='+self.out_dir+'darks/'+self.camwbtag+'.summed.0000001'
           printf, lun, '    DARK_NUM=0000001'
           printf, lun, '    ' + wclip
           if((upref[pp] EQ '8542' OR upref[pp] EQ '7772' ) AND (keyword_set(descatter))) then begin
              psff = self.descatter_dir+'/'+self.camwbtag+'.psf.f0'
              bgf = self.descatter_dir+'/'+self.camwbtag+'.backgain.f0'
              if(file_test(psff) AND file_test(bgf)) then begin
                 printf, lun, '    PSF='+psff
                 printf, lun, '    BACK_GAIN='+bgf
              endif
           endif 

           if(keyword_set(div)) then begin
              printf, lun, '    DIVERSITY='+string(div[0])+' mm'
           endif

           xofile = (file_search(self.out_dir+'/calib/'+self.camwbtag+'.*.xoffs'))[0]
           yofile = (file_search(self.out_dir+'/calib/'+self.camwbtag+'.*.yoffs'))[0]
           
           if(file_test(xofile)) then printf, lun, '    XOFFSET='+xofile
           if(file_test(yofile)) then printf, lun, '    YOFFSET='+yofile
           
                                ; printf, lun, '    INCOMPLETE'
           if(n_elements(nfac) gt 0) then printf,lun,'    NF=',red_stri(nfac[0])
           printf, lun, '  }'
           printf, lun, '}'

           ;; Loop all wavelengths
           pos1 = where((ustatp eq upref[pp]), count)
           if(count eq 0) then continue
           ustat1 = ustat[pos1]

           for ii = 0L, count - 1 do BEGIN
              
              ;; External states?
              if(keyword_set(state)) then begin
                 dum = where(state eq ustat1[ii], cstate)
                 if(cstate eq 0) then continue
                 print, inam+' : found '+state+' -> scan = '+uscan[ss]
              endif

              self -> whichoffset, ustat1[ii], xoff = xoff, yoff = yoff

              ;; Trans. camera
              istate = red_encode_scan(hscans[pos1[ii]], scan)+'.'+ustat1[ii]

              ;; lc4?
              tmp = strsplit(istate,'.', /extract)
              ntmp = n_elements(tmp)

              idx = strsplit(ustat1[ii],'.')
              nidx = n_elements(idx)
              iwavt = strmid(ustat1[ii], idx[0], idx[nidx-1]-1)

              if(keyword_set(skip)) then begin
                 dum = where(iwavt eq skip, ccout)
                 if ccout ne 0 then begin
                    print, inam+' : skipping state -> '+ustat1[ii]
                    continue
                 endif
              endif

              printf, lun, 'object{'
              printf, lun, '  WAVELENGTH=' + lam
              printf, lun, '  OUTPUT_FILE=results/'+self.camttag+'.'+istate 
              if(n_elements(weight) eq 3) then printf, lun, '  WEIGHT='+string(weight[1])
              printf, lun, '  channel{'
              printf, lun, '    IMAGE_DATA_DIR='+self.out_dir+'/data/'+folder_tag+ '/' +self.camttag+'/'
              printf, lun, '    FILENAME_TEMPLATE='+self.camttag+'.'+istate+'.%07d'
                                ;  printf, lun, '    DIVERSITY=0.0 mm'
              if(~keyword_set(unpol)) then begin
                 if(keyword_set(newgains)) then begin
                    search = self.out_dir+'/gaintables/'+folder_tag+'/'+self.camttag + '.' + istate+'.gain'
                 endif else begin
                    search = self.out_dir+'/gaintables/'+self.camttag + '.' + ustat1[ii] + '*.gain'
                 endelse
              endif Else begin

                 search = self.out_dir+'/gaintables/'+self.camttag + $
                          '.' + strmid(ustat1[ii], idx[0], $
                                       idx[nidx-1])+ '*unpol.gain'
                                ;if tmp[ntmp-1] eq 'lc4' then search = self.out_dir+'/gaintables/'+$
                                ;                                      self.camttag + '.' + ustat[pos[ii]] + $
                                ;                                      '*.gain'
                                ;
              endelse
              printf, lun, '    GAIN_FILE=' + file_search(search)
              printf, lun, '    DARK_TEMPLATE='+self.out_dir+'/darks/'+self.camttag+'.summed.0000001'
              printf, lun, '    DARK_NUM=0000001'
              printf, lun, '    ' + tclip

              xofile = self.out_dir+'/calib/'+self.camttag+'.'+xoff
              yofile = self.out_dir+'/calib/'+self.camttag+'.'+yoff
              if(file_test(xofile)) then printf, lun, '    XOFFSET='+xofile
              if(file_test(yofile)) then printf, lun, '    YOFFSET='+yofile

              if((upref[pp] EQ '8542' OR upref[pp] EQ '7772' ) AND (keyword_set(descatter))) then begin
                 psff = self.descatter_dir+'/'+self.camttag+'.psf.f0'
                 bgf = self.descatter_dir+'/'+self.camttag+'.backgain.f0'
                 if(file_test(psff) AND file_test(bgf)) then begin
                    printf, lun, '    PSF='+psff
                    printf, lun, '    BACK_GAIN='+bgf
                 endif
              endif 

              if(keyword_set(div)) then begin
                 printf, lun, '    DIVERSITY='+string(div[1])+' mm'
              endif
              if(n_elements(nfac) gt 0) then printf,lun,'    NF=',red_stri(nfac[1])
           
              printf, lun, '    INCOMPLETE'
              printf, lun, '  }'
              printf, lun, '}'  

              ;; Reflected camera
              printf, lun, 'object{'
              printf, lun, '  WAVELENGTH=' + lam
              printf, lun, '  OUTPUT_FILE=results/'+self.camrtag+'.'+istate 
              if(n_elements(weight) eq 3) then printf, lun, '  WEIGHT='+string(weight[2])
              printf, lun, '  channel{'
              printf, lun, '    IMAGE_DATA_DIR='+self.out_dir+'/data/'+folder_tag+ '/' +self.camrtag+'/'
              printf, lun, '    FILENAME_TEMPLATE='+self.camrtag+'.'+istate+'.%07d'
                                ;   printf, lun, '    DIVERSITY=0.0 mm' 
              if(~keyword_set(unpol)) then begin
                 if(keyword_set(newgains)) then begin
                    search = self.out_dir+'/gaintables/'+folder_tag+'/'+self.camrtag + '.' + istate+'.gain'
                 endif else begin
                    search = self.out_dir+'/gaintables/'+self.camrtag + '.' + ustat1[ii] + '*.gain'
                 endelse
              endif Else begin
                 idx = strsplit(ustat1[ii],'.')
                 nidx = n_elements(idx)
                 search = file_search(self.out_dir+'/gaintables/'+self.camrtag + $
                                      '.' + strmid(ustat1[ii], idx[0], $
                                                   idx[nidx-1])+ '*unpol.gain')
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
              if((upref[pp] EQ '8542' OR upref[pp] EQ '7772' ) AND (keyword_set(descatter))) then begin
                 psff = self.descatter_dir+'/'+self.camrtag+'.psf.f0'
                 bgf = self.descatter_dir+'/'+self.camrtag+'.backgain.f0'
                 if(file_test(psff) AND file_test(bgf)) then begin
                    printf, lun, '    PSF='+psff
                    printf, lun, '    BACK_GAIN='+bgf
                 endif
              endif 

              if(keyword_set(div)) then begin
                 printf, lun, '    DIVERSITY='+string(div[2])+' mm'
              endif
              if(n_elements(nfac) gt 0) then printf,lun,'    NF=',red_stri(nfac[2])

              printf, lun, '    INCOMPLETE'
              printf, lun, '  }'
              printf, lun, '}'

              ;; WB with states (for de-warping to the anchor, only to
              ;; remove rubbersheet when differential seeing is
              ;; strong)
              if(keyword_set(wb_states)) then begin
                 printf, lun, 'object{'
                 printf, lun, '  WAVELENGTH=' + lam
                 printf, lun, '  WEIGHT=0.00'
                 printf, lun, '  OUTPUT_FILE=results/'+self.camwbtag+'.'+istate 
                 printf, lun, '  channel{'
                 printf, lun, '    IMAGE_DATA_DIR='+self.out_dir+'/data/'+folder_tag+ '/' +self.camwbtag+'/'
                 printf, lun, '    FILENAME_TEMPLATE='+self.camwbtag+'.'+istate+'.%07d'
                                ; printf, lun, '    DIVERSITY=0.0 mm'
                 printf, lun, '    GAIN_FILE=' + file_search(self.out_dir+'/gaintables/'+self.camwbtag + $
                                                             '.' + upref[pp] + '*.gain')
                 printf, lun, '    DARK_TEMPLATE='+self.out_dir+'/darks/'+self.camwbtag+'.summed.0000001'
                 printf, lun, '    DARK_NUM=0000001'
                 printf, lun, '    ' + wclip
                 
                 if((upref[pp] EQ '8542' OR upref[pp] EQ '7772' ) AND (keyword_set(descatter))) then begin
                    psff = self.descatter_dir+'/'+self.camwbtag+'.psf.f0'
                    bgf = self.descatter_dir+'/'+self.camwbtag+'.backgain.f0'
                    if(file_test(psff) AND file_test(bgf)) then begin
                       printf, lun, '    PSF='+psff
                       printf, lun, '    BACK_GAIN='+bgf
                    endif
                 endif 

                 if(keyword_set(div)) then begin
                    printf, lun, '    DIVERSITY='+string(div[0])+' mm'
                 endif
                 xofile = self.out_dir+'/calib/'+self.camwbtag+'.'+xoff
                 yofile = self.out_dir+'/calib/'+self.camwbtag+'.'+yoff
                 if(file_test(xofile)) then printf, lun, '    XOFFSET='+xofile
                 if(file_test(yofile)) then printf, lun, '    YOFFSET='+yofile
                 if(n_elements(nfac) gt 0) then printf,lun,'    NF=',red_stri(nfac[0])

                 printf, lun, '    INCOMPLETE'
                 printf, lun, '  }'
                 printf, lun, '}'
              endif          
           endfor

           ;; Global keywords
           printf, lun, 'PROG_DATA_DIR=./data/'
           printf, lun, 'DATE_OBS='+date_obs
           printf, lun, 'IMAGE_NUMS='+nall ;n0+'-'+n1
           printf, lun, 'BASIS=Karhunen-Loeve'
           printf, lun, 'MODES='+modes
           printf, lun, 'NUM_POINTS='+numpoints
           printf, lun, 'TELESCOPE_D=0.97'
           printf, lun, 'ARCSECPERPIX=0.0592'
           printf, lun, 'PIXELSIZE=16.0E-6'
           printf, lun, 'GETSTEP=getstep_conjugate_gradient'
           printf, lun, 'GRADIENT=gradient_diff'
           printf, lun, 'MAX_LOCAL_SHIFT=30'
           printf, lun, 'NEW_CONSTRAINTS'
           printf, lun, 'FILE_TYPE='+outformat
           printf, lun, 'FAST_QR'
           if(outformat eq 'MOMFBD') then printf, lun, 'GET_PSF'
           if(outformat eq 'MOMFBD') then printf, lun, 'GET_PSF_AVG'

           ;; External keywords?
           if(keyword_set(global_keywords)) then begin
              nk = n_elements(global_keywords)
              for ki = 0L, nk -1 do printf, lun, global_keywords[ki]
           endif

           free_lun, lun
        endfor
     endfor
  endfor

  print, inam+' : done!'
  return
end
