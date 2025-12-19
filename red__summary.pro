;+
; Make a document with dataset summary and quicklook info.
; 
; Needs pdflatex to be installed, as well as one of pdfunite and
; pdftk.
; 
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
;   force : in, optional, type=boolean
;   
;     Normally quicklook/quicklook_mosaic is run only if there are no
;     quicklook images. This keyword forces a rerun. Option
;     /core_and_wings is used by default, unless there are extra
;     keywords to be passed on to quicklook/quicklook_mosaic.
; 
; 
; :History:
;   
;   2025-12-07 : MGL. First version.
; 
;-
pro red::summary, _extra = extra $
                  , force = force $
                  , no_mosaics = no_mosaics $
                  , no_regular = no_regular 


  pwd = getenv('PWD')

  ;; Check that we have pdflatex and pdfunite
  spawn, 'pdflatex -v', stdout, stderr, exit_status = status_pdflatex  
  spawn, 'pdfunite -v', stdout, stderr, exit_status = status_pdfunite
  spawn, 'pdftk -v', stdout, stderr, exit_status = status_pdftk
  if status_pdflatex ne 0 then red_message, 'pdflatex does not seem to be available.'
  if status_pdfunite ne 0 then red_message, 'pdfunite does not seem to be available.'
  if status_pdftk    ne 0 then red_message, 'pdftk does not seem to be available.'
  if status_pdflatex ne 0 then begin
    red_message, 'We need to run pdflatex.'
    return
  endif
  

  latex_options = ''
  self -> cameras $
     , instrument = instrument 

  data_dirs = *self.data_dirs
  Ndirs = n_elements(data_dirs)

  timestamps = file_basename(data_dirs)

  isMosaic = strmatch(data_dirs, '*mosaic/*')

  for idir = 0, Ndirs-1 do begin
  
    ;; First run quicklook
    if isMosaic(idir) then begin
      if keyword_set(no_mosaics) then continue
      qdir = 'quicklook_mosaic/'+timestamps[idir]+'/'
    endif else begin
      if keyword_set(no_regular) then continue
      qdir = 'quicklook/'+timestamps[idir]+'/'
    endelse
    
    ;; Also summarize 
    self -> summarize_datadir, timestamps[idir]

    print, data_dirs[idir]
    
    texfile = timestamps[idir]+'.tex'
    openw, lun, /get_lun, 'summaries/'+texfile

    printf, lun, '\pdfminorversion=7'
    printf, lun, '\documentclass['+latex_options+']{article}'
    printf, lun, '\usepackage{graphicx} '
    printf, lun, '\usepackage{time}'
    printf, lun, '\usepackage{fancyvrb}'
    printf, lun, '\usepackage{amsmath}'
    printf, lun, '\usepackage[labelformat=simple]{subcaption}'
    printf, lun, '\renewcommand\thesubfigure{(\alph{subfigure})}'
    printf, lun, '\usepackage[a4paper]{geometry}'
    printf, lun, ''
    printf, lun, '\begin{document}'
    printf, lun, ''

    printf, lun, '\section*{'+instrument+' '+self.isodate+'T'+timestamps[idir]+'}'

    ;; Summary tables

    self -> summarize_datadir, timestamps[idir]
    
    txtfile = file_search('summaries/'+timestamps[idir]+'.txt', count = Ntxt)
    savfile = file_search('summaries/'+timestamps[idir]+'.sav', count = Nsav)

    case 1 of
      Nsav : begin
        restore, savfile
        printf, lun, '\noindent'
        printf, lun, '\begin{table}[!th]'
        printf, lun, sum_struct.datadir+'\\'
        printf, lun, 'Beginning and end times: ' + sum_struct.started + ' ' + sum_struct.ended+'\\'
        printf, lun, 'Cameras: ' + strjoin(reverse(sum_struct.cameras), ', ')+'\\'
        if isMosaic(idir) then begin
          printf, lun, 'Number of tiles: ' + sum_struct.Ntiles+'\\'
          ttot = red_time2double(sum_struct.ended) - red_time2double(sum_struct.started) ;  [s]
          printf, lun, 'Total time: ' + red_time2string(ttot)+'\\'
        endif else begin
          printf, lun, 'Number of scans: ' + sum_struct.Nscans+'\\'
          printf, lun, 'Number of frames per scan: ' + sum_struct.Nframes_per_scan+'\\'
          tscan = ( red_time2double(sum_struct.ended) - red_time2double(sum_struct.started) ) / long(sum_struct.Nscans)
          case 1 of
            tscan lt 2 : format = '(f8.2)'
            tscan lt 10 : format = '(f8.1)'
            else : format = '(i5)'
          endcase
          printf, lun, 'Approximate time per scan: ' + strtrim(string(tscan, format = format), 2)+ ' s.\\[10mm]'
        endelse

        Npref = n_elements(sum_struct.prefilters)
        
        printf, lun, ''

        for ipref = 0, Npref-1 do begin
          Ntun = n_elements((sum_struct.nb[ipref]).tunings)
          printf, lun, '\begin{minipage}[t]{0.25\textwidth}'
          printf, lun, '\vspace{0pt}\noindent'
          printf, lun, '\textbf{'+sum_struct.prefilters[ipref]+'\hfill}\\[3mm]'
          printf, lun, '\begin{tabular}{rrr}'
          printf, lun, '\hline'
          printf, lun, 'Tuning & LC & $N_{\text{exp}}$ \\'
          printf, lun, '\hline'
          for itun = 0, Ntun-1 do begin
            printf, lun, '$'+red_strreplace((sum_struct.nb[ipref]).tunings[itun],sum_struct.prefilters[ipref]+'_','') + '$'  $
                    + ' & ' $
                    + red_strreplace(rdx_ints2str([(sum_struct.nb[ipref]).lc[itun]]), '-', '--') $
                    + ' & ' $
                    + strtrim((sum_struct.nb[ipref]).nexp[itun], 2) $
                    + ' \\'
          endfor                ; itun
          printf, lun, '\hline'
          
          printf, lun, '\end{tabular}'
          printf, lun, '\end{minipage}\quad'
        endfor                  ; ipref
        printf, lun, '\end{table}'
        printf, lun, '\clearpage'
      end
      Ntxt : begin
        printf, lun, '\VerbatimInput[fontsize=\tiny]{'+pwd+'/'+sumfile[0]+'}'
      end
      else :
    endcase

    
    ;; Images
    if isMosaic(idir) then begin

      imfiles = file_search(qdir + '*png', count = Nims)

      if Nims eq 0 then begin
        if n_elements(extra) gt 0 then begin
          self -> quicklook_mosaic, datasets = timestamps[idir], _strict_extra = extra
        endif else begin
          self -> quicklook_mosaic, datasets = timestamps[idir], /core_and_wing
        endelse
        imfiles = file_search(qdir + '*png', count = Nims)
        if Nim eq 0 then stop
        endif

        
        tuninfo = stregex(file_basename(imfiles), '([0-9][0-9][0-9][0-9])[._]([0-9][0-9][0-9][0-9])_([+-][0-9]*)', /extract, /sub)
        prefs = tuninfo[1,*]
        for iim = 0, Nims-1 do begin
        splt = strsplit(file_basename(imfiles[iim]),'_.',/extract)
        cap = strjoin(splt[[0, 2, 3, 4, 5, 6, 7]], ' ')
        cap = strjoin(tuninfo[1:*, iim], ' ')

        printf, lun, '\noindent'
        printf, lun, '\begin{figure}[!th]'
        printf, lun, '\centering'
        printf, lun, '\includegraphics[width=\linewidth]{'+pwd+'/'+imfiles[iim]+'}'
        printf, lun, '\caption{'+cap+'.}'
        printf, lun, '\end{figure}'
        
      endfor                    ; iim

    endif else begin

      imfiles = file_search(qdir + '*jpg', count = Nims)

      if Nims eq 0 then begin
        if n_elements(extra) gt 0 then begin
          self -> quicklook, datasets = timestamps[idir], min_nscan = 1, _strict_extra = extra
        endif else begin
          self -> quicklook, datasets = timestamps[idir], min_nscan = 1, /core_and_wing
        endelse
        imfiles = file_search(qdir + '*jpg', count = Nims)
        if Nims eq 0 then stop
      endif

      red_extractstates, imfiles, /base, wav = wav, dwav = dwav, pref = pref
      upref = red_uniquify(pref, count = Nprefs)

      for ipref = 0, Nprefs-1 do begin
        indx = where(pref eq upref[ipref], Nim)
        iindx = sort(dwav[indx])
        
        printf, lun, '\noindent'
        printf, lun, '\begin{figure}[!th]'
        printf, lun, '\centering'
        
        for iim = 0, Nim-1 do begin
          printf, lun, '\begin{subfigure}{0.32\linewidth}'
          printf, lun, '\includegraphics[width=\linewidth]{'+pwd+'/'+imfiles[indx[iindx[[iim]]]]+'}'
          cap = red_strreplace(wav[indx[iindx[[iim]]]], '_', '')
          cap = red_strreplace(cap, '+', '$+$')
          cap = red_strreplace(cap, '-', '$-$')
          printf, lun, '\caption{'+cap+'}'
          printf, lun, '\end{subfigure}'
        endfor

;        printf, lun, '\caption{'+upref[ipref]+'.}'
        printf, lun, '\end{figure}'
          
      endfor                    ; ipref
    endelse
    
    printf, lun, ''
    printf, lun, ''

    ;; r0 plot
    r0file = file_search(qdir + '*r0*pdf', count = Nr0)
    if Nr0 ne 0 then begin
      printf, lun, '\noindent'
      printf, lun, '\begin{figure}[!th]'
      printf, lun, '\centering'
      printf, lun, '\includegraphics[width=\linewidth]{'+pwd+'/'+r0file[0]+'}'
      printf, lun, '\caption{'+timestamps[idir]+' line scans $r_0$.}'
      printf, lun, '\end{figure}'
    endif


    printf, lun, '\end{document}'
    free_lun, lun

    cd, 'summaries/'
    spawn, 'pdflatex '+texfile

    cd, pwd

  endfor                        ; idir

  cd, 'summaries/'
  red_message, 'The summary will be combined into a single file, summaries/summary.pdf'
  case 0 of

    status_pdfunite : spawn, 'pdfunite ??:??:??.pdf summary.pdf'

    status_pdftk    : spawn, 'pdftk ??:??:??.pdf cat output summary.pdf'

    else : begin
      red_message, ['We would need pdfunite or pdftk to combine all the summaries/??:??:??.pdf' $
                    , 'files into a single document. The individual pdf files are created, though.']
    end
    
  endcase
  
  cd, pwd
  
end

cd, '/scratch_beegfs/mats/NEW/2025-08-18/CRISP2'
a = crisp2red(/dev, /no)
a -> summary                    ;   , /no_regular

end
