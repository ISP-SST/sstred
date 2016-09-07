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
; :Author:
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
;    nostate  :  in, optional, type=boolean
;   
;       For wideband data, if not set, the mfbd sets will be set up as
;       if the data had states.
;   
;    numpoints  :  in, optional, 
;   
;   
;   
;    modes  :  in, optional, 
;   
;   
;   
;    date_obs  :  in, optional, type=string
;   
;      The date of observations in ISO (YYYY-MM-DD) format. If not
;      given, taken from class object.
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
;    pref : in, optional
;   
;        The prefilter for which to generate mfbd config files. 
;   
;    nremove :  in, optional, 
;   
;    camera : in, optional, type=string, default="Crisp-W or Chromis-W"
;
;       The camera to run mfbd on.
;   
;    mfbddir :  in, optional, type=string, default='mfbd'
;   
;       Top directory of output tree.
;   
;    nimages : in, optional, type=integer, default="One scan"
;
;       The number of images per mfbd data set.
;
; :History:
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
;   2016-05-31 : MGL. Default date from class object. Better loop
;                indices. Added dirs keyword.
;
;   2016-06-01 : MGL. Replace CRISP specific code.
;
;   2016-06-02 : MGL. Use new class methods and camerainfo.
;
;   2016-06-05 : MGL. New keyword camera. Default cam that works
;                for both CRISP and CHROMIS.
;
;   2016-06-06 : MGL. Minor improvements from prepmomfbd.
;
;   2016-06-08 : MGL. Renamed keyword nf to nfac so scope_varfetch
;                does not get confused in writelog.
;
;   2016-06-09 : MGL. Bugfix in the search string.
;
;   2016-08-23 : THI. Rename camtag to detector and channel to camera,
;                so the names match those of the corresponding SolarNet
;                keywords.
;
;   2016-09-01 : MGL. Make it work in terms of state rather than
;                prefilter. 
;
;   2016-09-07 : MGL. New keyword nostate.
;
;-
pro red::prepmfbd, numpoints = numpoints, $
                   modes = modes, $
                   date_obs = date_obs, $
                   dirs = dirs, $
                   no_descatter = no_descatter, $
                   global_keywords = global_keywords, $
                   pref = pref, $
                   nfac = nfac, $
                   nimages = nimages, $
                   mfbddir = mfbddir, $
                   camera = camera, $
                   nostate = nostate

  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  ;; Logging
  help, /obj, self, output = selfinfo 
  red_writelog, selfinfo = selfinfo

  ;; Get keywords
  if n_elements(date_obs) eq 0 then date_obs = self.isodate
  if n_elements(modes) eq 0 then modes = '2-29,31-36,44,45'
  if n_elements(mfbddir) eq 0 then mfbddir = 'mfbd' 

  if n_elements(dirs) gt 0 then begin
     dirs = [dirs] 
  endif else begin
     if ~ptr_valid(self.data_dirs) then begin
        print, inam+' : ERROR : undefined data_dir'
        return
     endif
     dirs = *self.data_dirs
  endelse

  Ndirs = n_elements(dirs)
  if Ndirs eq 0 then begin
     print, inam+' : ERROR : no directories defined'
     return
  endif else begin
     if Ndirs gt 1 then dirstr = '['+ strjoin(dirs,';') + ']' $
     else dirstr = dirs[0]
  endelse

  if n_elements(camera) ne 0 then begin
     ;; Check that the specified camera exists
     campos = where(strmatch(*self.cameras,camera), nn)
     if nn eq 1 then begin
        cam = camera
     endif else begin
        print, inam + ' : ERROR -> camera not defined '+camera
        retall
     endelse
  endif else begin
     ;; Select as default the wideband camera
     campos = where(strmatch(*self.cameras,'*-W'), nn)
     if nn eq 1 then begin
        cam = (*self.cameras)[campos]
     endif else begin
        print, inam + ' : ERROR -> camera not specified and no wideband camera defined.'
        retall        
     endelse
  endelse

  detector = self -> getdetector(cam)
  
  if n_elements(numpoints) eq 0 then begin
     ;; About the same subfield size in arcsec as CRISP:
     numpoints = strtrim(round(88*0.0590/self.image_scale/2)*2, 2)
  endif else begin
     ;; Convert strings, just to avoid breaking existing codes.
     if( size(numpoints, /type) eq 7 ) then numpoints = fix(numpoints) 
  endelse

  for idir = 0L, Ndirs - 1 do begin

     data_dir = dirs[idir]
     folder_tag = file_basename(data_dir)

     ;; Directory
     if (strmatch(cam,'*-W') or strmatch(cam,'*-D')) $
        and keyword_set(nostate) then begin
        ;; Wideband data
        camdir = data_dir + '/' + cam + '_nostate/'
     endif else begin
        ;; Narrowband data
        camdir = data_dir + '/' + cam + '/'
     endelse

     search = camdir + detector
     files = file_search(search+'*', count = Nfiles) 

     IF Nfiles EQ 0 THEN BEGIN
         print, inam + ' : ERROR -> no frames found in '+search
         print, inam + '   Did you run link_data?'
         return
     ENDIF 

     files = red_sortfiles(temporary(files))
     
     ;; Get image unique states
     self -> extractstates, files, states

     ;; Get unique prefilters
     upref = states[uniq(states.prefilter, sort(states.prefilter))].prefilter
     Npref = n_elements(upref)

     ;; Get fullstates
     ustate = states[uniq(states.fullstate, sort(states.fullstate))].fullstate
     Nstates = n_elements(ustate)

     ;; Get scan numbers
     uscan = states[uniq(states.scannumber, sort(states.scannumber))].scannumber
     Nscans = n_elements(uscan)

     ;; Create a config file per state and scan number

  
     ;; Loop over scans
     for iscan = 0L, Nscans-1 do begin
        
        scan = string(uscan[iscan], format = '(i05)')

        ;; Loop over prefilters
        for istate = 0L, Nstates - 1 do begin
           
;           if n_elements(pref) ne 0 then begin
;              if upref[ipref] NE pref then begin
;                 print, inam + ' : Skipping prefilter -> ' + upref[ipref]
;                 continue
;              endif
;           endif


           ;; Where to put the config files:
           outdir = self.out_dir + '/' + mfbddir + '/' + folder_tag $
                    + '/' + ustate[istate] + '/cfg/'

           ;; Needed by momfbd:
           file_mkdir, outdir + 'results'
           file_mkdir, outdir + 'data'

           ;; Image numbers for the current scan:
           numpos = where((states.scannumber eq scan) $
;                          AND (stat.star eq 0B) $
                          AND (states.fullstate eq ustate[istate]) $
                          , Ncount)
           if(Ncount eq 0) then continue
    
           ;; Split scan into subsets?
           if n_elements(Nimages) eq 0 then begin
              Nsubscans = 1
           endif else begin
              Nsubscans = round(Ncount/float(Nimages))
           endelse

           ;; Loop over subsets of the current scan
           for isub = 0L, Nsubscans-1 do begin

              if Nsubscans eq 1 then begin
                 cfg_file = 'momfbd.reduc.' + ustate[istate] + '.' + scan + '.cfg'
                 n0 = strtrim(states[numpos[0]].framenumber, 2)
                 n1 = strtrim(states[numpos[ncount-1]].framenumber, 2)
                 nall = strjoin([n0,n1],'-')
                 print, inam+' : State = ' + ustate[istate] + ' -> scan = ' $
                        + scan + ' -> image range = [' + nall + ']'
                 oname = detector + '_' + scan + '_' + ustate[istate]
             endif else begin
                 subscan = red_stri(isub, ni = '(i03)')
                 cfg_file = 'momfbd.reduc.' + ustate[istate] + '.' + scan + ':' $
                            + subscan + '.cfg'
                 n0 = strtrim(states[numpos[isub*ncount/Nsubscans]].framenumber, 2)
                 n1 = strtrim(states[numpos[((isub+1)*ncount/Nsubscans-1) $
                                            <ncount]].framenumber, 2)
                 nall = strjoin([n0,n1],'-')                 
                 print, inam+' : State = ' + ustate[istate] + ' -> scan = ' $
                        + scan + ':' + subscan + ' -> image range = [' + nall + ']'
                 oname = detector + '_' + scan + ':' + subscan + '_' + ustate[istate]
              endelse

              ;; Open config file for writing
              openw, lun, outdir + cfg_file, /get_lun, width=2500

              self -> get_calib, states[numpos[0]], darkname = darkname, gainname = gainname
              caminfo = red_camerainfo(detector)
              
              printf, lun, 'object{'
              printf, lun, '  WAVELENGTH=' + strtrim(states[0].pf_wavelength, 2)
              printf, lun, '  OUTPUT_FILE=results/'+oname
              printf, lun, '  channel{'
              printf, lun, '    IMAGE_DATA_DIR=' + self.out_dir + camdir
              if (strmatch(cam,'*-W') or strmatch(cam,'*-D')) $
                 and keyword_set(nostate)  then begin
                 ;; Wideband, no state
                 template = detector + '_' + scan + '_' + upref[ipref] + '_%07d.fits'
              endif else begin
                 ;; With state
                 template = detector + '_' + scan + '_' + states[numpos[0]].fullstate $
                            + '_%07d.fits'
              endelse
              printf, lun, '    FILENAME_TEMPLATE=' + template
              printf, lun, '    GAIN_FILE=' + gainname
              printf, lun, '    DARK_TEMPLATE=' + darkname
              printf, lun, '    DARK_NUM=0000001'

;              if upref[ipref] EQ '8542' OR upref[ipref] EQ '7772' $
;                 AND ~keyword_set(no_descatter) then begin
;                 self -> loadbackscatter, cam, upref[ipref], bg, psf $
;                                          , bgfile = bgf, bpfile = psff
;                 printf, lun, '    PSF='+psff
;                 printf, lun, '    BACK_GAIN='+bgf
;              endif 

              printf, lun, '    INCOMPLETE'

              if(n_elements(nfac) gt 0) then printf,lun,'    NF=',red_stri(nfac[0])

              printf, lun, '  }'
              printf, lun, '}'

              ;; Global keywords
              printf, lun, 'PROG_DATA_DIR=./data/'
              printf, lun, 'DATE_OBS=' + date_obs
              printf, lun, 'IMAGE_NUMS=' + nall 
              printf, lun, 'BASIS=Karhunen-Loeve'
              printf, lun, 'MODES=' + modes
              printf, lun, 'NUM_POINTS=' + strtrim(numpoints, 2)
              printf, lun, 'TELESCOPE_D=0.97'
              printf, lun, 'ARCSECPERPIX=' + self.image_scale
              printf, lun, 'PIXELSIZE=' + strtrim(caminfo.pixelsize, 2)
              printf, lun, 'GETSTEP=getstep_conjugate_gradient'
              printf, lun, 'GRADIENT=gradient_diff'
              printf, lun, 'MAX_LOCAL_SHIFT=30'
              printf, lun, 'NEW_CONSTRAINTS'
              printf, lun, 'FILE_TYPE=' + self.filetype
              printf, lun, 'FAST_QR'
              printf, lun, 'FIT_PLANE'
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
        endfor                  ; istate

     endfor                     ; iscan
  endfor                        ; idir

  print, inam+' : done!'
  return
end
