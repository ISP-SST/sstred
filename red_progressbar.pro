; docformat = 'rst'

;+
; Display progress in teminal window.
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, ISP
; 
; 
; :Params:
; 
;    i : in, type="scalar number"
; 
;       The progress so far.
;
;    message : in, optional, type=string, default="'Progress'"
;   
;       A message identifying the task for which progress is being
;       made. 
; 
;    N : in, type="scalar number"
;    
;       Normalize i with this number.
; 
; :Keywords:
; 
;    nobar : in, optional, type=boolean
;
;       If set, do not show a status bar.
;
;    barlength : in, optional, type=integer, default=20.
;
;       The length of the status bar in number of characters.
;
;    predict : in, optional, type=boolean
;
;       Print a prediction for the remaining time. 
; 
; :History:
; 
;     2016-05-25 : MGL. First version, code taken from
;                  red_sumfiles.pro. 
; 
;     2016-10-10 : MGL. Implemented a real bar. 
; 
;     2016-10-25 : MGL. Deal with N=0 and N=1.
; 
;     2016-10-28 : MGL. Don't need to send the message when finishing.
; 
;     2016-10-31 : MGL. Remove keyword finished, just test i vs. N for
;                  completion. Make message a parameter instead of a
;                  keyword. New keyword clock.
; 
;     2016-11-16 : MGL. New keyword predict.
; 
;     2016-11-28 : MGL. Call red_timestring with /interval.
; 
;     2017-05-08 : MGL. Try to make predictions more robust. Adapt
;                  output line to teminal width.
; 
;     2017-10-24 : MGL. Replace keyword clock with a common block. 
; 
;     2018-01-22 : MGL. Allow not calling red_progressbar for every
;                  iteration. 
; 
;-
pro red_progressbar, i, N, message $
                     , barlength = barlength $
                     , nobar = nobar $
                     , predict = predict

  common red_progressbar_common, clock, times, iterations, icnt

  ;; Initialize
  if i eq 0 or n_elements(clock) eq 0 then begin

    clock = tic()
    icnt = i

    times = dblarr(N)
    iterations = lindgen(N)
    
  endif  
  
  if n_elements(message) eq 0 then message = 'Progress'
  bb = string(13B)              ; CR w/o LF

  case N of
    0: return
    1: norm = 100.
    else: norm = 100. / (N - 1.0)  
  endcase

  if i eq N-1 then percentdone = 100. else percentdone = norm * i
  prediction = '' 

  times[i] = toc(clock)
  iterations[i] = i

  outlength = (TERMINAL_SIZE( ))[0] ; Base output length on terminal width
  outline = string(replicate(32B, outlength))

  if n_elements(barlength) eq 0 then barlength = 20
  prediction = ''

  if keyword_set(predict) and i gt 1 then begin

    indx = where(times ne 0.0, Ntimed)
    if Ntimed gt 1 then begin
      
      dt = median( (times[indx[1:*]]-times[indx[0:Ntimed-2]]) / (iterations[indx[1:*]]-iterations[indx[0:Ntimed-2]]) )
      
      
      time_remaining = (N-i) * dt
      
      prediction = ' ('+red_timestring(round(time_remaining), Nsecdec = 0, /interval)+' remaining)'
    endif
    
  endif

  time = 'in ' + red_timestring(round(toc(clock)), Nsecdec = 0, /interval) + prediction

  if keyword_set(nobar) then begin
    
    strput, outline, string( message + ' -> ' $
                             , percentdone, '% ' + time + '   ', FORMAT = '(A,F5.1,A,$)')
    
  endif else begin
    
    time += ': '
    
    elength = floor(percentdone/100*barlength)
    mlength = barlength-elength
    bar = ''
    if elength gt 0 then bar += string(replicate(61B, elength)) ; Replicated '='
    if mlength gt 0 then bar += string(replicate(45B, mlength)) ; Replicated '-'
    strput, outline, string('[' + bar + '] ' $
                            , percentdone, '% ' + time + message + '   ', FORMAT = '(A,F5.1,A,$)')
    
  endelse

  
  print, bb, outline, FORMAT = '(A,A,$)'
  
  if i eq N-1 then begin
    undefine, clock
    print                       ; Finished
  endif

end


N=50000
;for i=0,N-1 do begin red_progressbar,i,N,'Test',clock=clock & wait,.1 & end    
;for i=0,N-1 do begin red_progressbar,i,N,'Test' & wait,.1 & end
;for i=0,N-1 do begin red_progressbar,i,N,'Test',clock=clock, /nobar& wait,.1 & end    
;for i=0,N-1 do begin red_progressbar,i,N,'Test', /nobar & wait,.1 & end

for i=0,N-1 do begin & if i mod 100 eq 0 or i eq N-1 then red_progressbar,i,N,'Test', /predict & wait,.001 & end    



for i=0,N-1 do begin red_progressbar,i,N,'Test', /predict & wait,.001 & end    
;for i=0,N-1 do begin red_progressbar,i,N,'Test',clock=clock, /predict & wait,.1 & end    

end
