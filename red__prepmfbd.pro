; docformat = 'rst'

;+
; Prepare for MFBD processing of WB data only using MFBD without the
; MO part.
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :author:
; 
; 
; :Keywords:
; 
;    wb_states  :  in, optional, 
;   
;   
;   
;    outformat  :  in, optional, 
;   
;   
;   
;    numpoints  :  in, optional, 
;   
;   
;   
;    modes  :  in, optional, 
;   
;   
;   
;    date_obs  :  in, optional, 
;   
;   
;   
;    no_descatter  :  in, optional, 
;   
;   
;   
;    global_keywords  :  in, optional, 
;   
;   
;   
;    skip  :  in, optional, 
;   
;   
;   
;    nremove :  in, optional, 
;   
;   
;   
;    mmfbddir :  in, optional, type=string, default='mfbd'
;   
;       Top directory of output tree.
;   
;    nimages : in, optional, type=integer, default="One scan"
;
;       The number of images per mfbd data set.
;
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
;
;   2013-06-13 : JdlCR. Added support for scan-dependent gains ->
;                using keyword "/newgains".
;
;   2013-06-28 : JdlCR. Added NF (object) option 
; 
;   2013-08-27 : MGL. Added support for logging. Let the subprogram
;                find out its own name.
;
;   2013-09-05 : MGL. Split off from red::prepmomfbd, stripped parts
;                that have to do with the narrowband cameras. New
;                keyword nimages. Change subdirectory to "mfbd".
;
;   2014-01-10 : MGL. Follow Pit's lead and remove keyword
;                outformat, use self.filetype instead. 
;
;   2014-04-02 : MGL. Fixed bug in IMAGE_DATA_DIR.
;
;   2014-05-13 : MGL. Remove keyword newgains.
;
;   2015-10-13 : MGL. Added keyword mfbddir. Allow keyword num_points
;                to not be a string.
;
;   2016-02-15 : MGL. Use red_loadbackscatter. Remove keyword descatter,
;                new keyword no_descatter.
;
;-
pro red::prepmfbd, numpoints = numpoints, $
                   modes = modes, date_obs = date_obs, no_descatter = no_descatter, $
                   global_keywords = global_keywords, skip = skip, $
                   pref = pref, $
                   nf = nfac, nimages = nimages, $
                   mfbddir = mfbddir

  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  ;; Logging
  help, /obj, self, output = selfinfo 
  red_writelog, selfinfo = selfinfo

  ;; Get keywords
  if n_elements(date_obs) eq 0 then begin
     date_obs = ' '
     read, date_obs, prompt = inam+' : type date_obs (YYYY-MM-DD): '
  endif
  if n_elements(numpoints) eq 0 then numpoints = '88'
  if n_elements(modes) eq 0 then modes = '2-29,31-36,44,45'
  if n_elements(nfac) eq 1 then nfac = replicate(nfac,3)
  if n_elements(mfbddir) eq 0 then mfbddir = 'mfbd' 
     
  for fff = 0, self.ndir - 1 do begin
     data_dir = (*self.data_list)[fff]
     spawn, 'find ' + data_dir + '/' + self.camwb+ '/ | grep im.ex | grep -v ".lcd."', files
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

     ntt = n_elements(ustat)
     hscans = stat.hscan[pos]

     ;; Get unique prefilters
     upref = stat.pref[uniq(stat.pref, sort(stat.pref))]
     np = n_elements(upref)

     ;; Get scan numbers
     uscan = stat.rscan[uniq(stat.rscan, sort(stat.rscan))]
     ns = n_elements(uscan)

     ;; Create a reduc file per prefilter and scan number?
     outdir = self.out_dir + '/' + mfbddir + '/' + folder_tag
     file_mkdir, outdir

     self -> getcamtags, dir = data_dir

     ;; Print cams
     print, ' WB  -> '+self.camwbtag

     ;; Loop over scans
     for ss = 0L, ns - 1 do begin
        
        scan = uscan[ss]

        ;; Loop over prefilters
        for pp = 0L, np - 1 do begin
           if(keyword_set(pref)) then begin
              if(upref[pp] NE pref) then begin
                 print, inam + ' : Skipping prefilter -> ' + upref[pp]
                 continue
              endif
           endif

           lam = strmid(string(float(upref[pp]) * 1.e-10), 2)

           outdir = self.out_dir + '/' + mfbddir + '/' + folder_tag + '/' + upref[pp] + '/cfg/'
           file_mkdir, outdir
           rdir = self.out_dir + '/' + mfbddir + '/' + folder_tag + '/' + upref[pp] + '/cfg/results/'
           file_mkdir, rdir
           ddir = self.out_dir + '/' + mfbddir + '/' + folder_tag + '/' + upref[pp] + '/cfg/data/'
           file_mkdir, ddir

           ;; Image numbers for the current scan
           numpos = where((stat.rscan eq uscan[ss]) AND (stat.star eq 0B) AND (stat.pref eq upref[pp]), ncount)
           if(ncount eq 0) then continue
    
           ;; Split scan into subsets?
           if n_elements(nimages) eq 0 then begin

              Nsubscans = 1

           endif else begin
              
              Nsubscans = round(ncount/float(nimages))
              
           endelse

           ;; Loop over subsets of the current scan
           for isub = 0, Nsubscans-1 do begin

              if Nsubscans eq 1 then begin
                 cfg_file = 'momfbd.reduc.'+upref[pp]+'.'+scan+'.cfg'
                 n0 = stat.nums[numpos[0]]
                 n1 = stat.nums[numpos[ncount-1]]
;                 nall = strjoin(stat.nums[numpos],',')
                 nall = strjoin([n0,n1],'-')
                 print, inam+' : Prefilter = '+upref[pp]+' -> scan = '+uscan[ss]+' -> image range = ['+n0+' - '+n1+']'
                 oname = self.camwbtag+'.'+scan+'.'+upref[pp]
             endif else begin
                 subscan = red_stri(isub, ni = '(i03)')
                 cfg_file = 'momfbd.reduc.'+upref[pp]+'.'+scan+':'+subscan+'.cfg'
                 n0 = stat.nums[numpos[isub*ncount/Nsubscans]]
                 n1 = stat.nums[numpos[((isub+1)*ncount/Nsubscans-1) <ncount]]
                 nall = strjoin([n0,n1],'-')                 
                 print, inam+' : Prefilter = '+upref[pp]+' -> scan = '+scan+':'+subscan+' -> image range = ['+n0+' - '+n1+']'
                 oname = self.camwbtag+'.'+scan+':'+subscan+'.'+upref[pp]
              endelse
              
              ;; Open config file for writing
              openw, lun, outdir + cfg_file, /get_lun, width=2500

              ;; WB only
              printf, lun, 'object{'
              printf, lun, '  WAVELENGTH=' + lam
              printf, lun, '  OUTPUT_FILE=results/'+oname
              printf, lun, '  channel{'
              printf, lun, '    IMAGE_DATA_DIR='+self.out_dir+'/data/'+folder_tag+ '/' +self.camwb+'_nostate/'
              ;printf, lun, '    IMAGE_DATA_DIR='+self.out_dir+'/data/'+folder_tag+ '/' +self.camwbtag+'_nostate/'
              printf, lun, '    FILENAME_TEMPLATE='+self.camwbtag+'.'+scan+'.'+upref[pp]+'.%07d'
                                ; printf, lun, '    DIVERSITY=0.0 mm'
              printf, lun, '    GAIN_FILE=' + file_search(self.out_dir+'gaintables/'+self.camwbtag + $
                                                          '.' + upref[pp]+'*.gain')
              printf, lun, '    DARK_TEMPLATE='+self.out_dir+'darks/'+self.camwbtag+'.summed.0000001'
              printf, lun, '    DARK_NUM=0000001'
              if (upref[pp] EQ '8542' OR upref[pp] EQ '7772' ) AND ~keyword_set(no_descatter) then begin
                 self -> loadbackscatter, self.camwbtag, upref[pp], bg, psf, bgfile = bgf, bpfile = psff
;                 psff = self.descatter_dir+'/'+self.camwbtag+'.psf.f0'
;                 bgf = self.descatter_dir+'/'+self.camwbtag+'.backgain.f0'
;                 if(file_test(psff) AND file_test(bgf)) then begin
                 printf, lun, '    PSF='+psff
                 printf, lun, '    BACK_GAIN='+bgf
;                 endif
              endif 

              printf, lun, '    INCOMPLETE'

              if(n_elements(nfac) gt 0) then printf,lun,'    NF=',red_stri(nfac[0])
              printf, lun, '  }'
              printf, lun, '}'

              ;; Global keywords
              printf, lun, 'PROG_DATA_DIR=./data/'
              printf, lun, 'DATE_OBS='+date_obs
              printf, lun, 'IMAGE_NUMS='+nall ;n0+'-'+n1
              printf, lun, 'BASIS=Karhunen-Loeve'
              printf, lun, 'MODES='+modes
              printf, lun, 'NUM_POINTS='+strtrim(numpoints, 2)
              printf, lun, 'TELESCOPE_D=0.97'
              printf, lun, 'ARCSECPERPIX='+self.image_scale
              printf, lun, 'PIXELSIZE=16.0E-6'
              printf, lun, 'GETSTEP=getstep_conjugate_gradient'
              printf, lun, 'GRADIENT=gradient_diff'
              printf, lun, 'MAX_LOCAL_SHIFT=30'
              printf, lun, 'NEW_CONSTRAINTS'
              printf, lun, 'FILE_TYPE='+self.filetype
              printf, lun, 'FAST_QR'
              if self.filetype eq 'MOMFBD' then begin
                 printf, lun, 'GET_PSF'
                 printf, lun, 'GET_PSF_AVG'
              endif
              ;; External keywords?
              if(keyword_set(global_keywords)) then begin
                 nk = n_elements(global_keywords)
                 for ki = 0L, nk -1 do printf, lun, global_keywords[ki]
              endif

              free_lun, lun
           endfor               ; isub
        endfor                  ; pp
     endfor                     ; ss
  endfor                        ; fff

  print, inam+' : done!'
  return
end
