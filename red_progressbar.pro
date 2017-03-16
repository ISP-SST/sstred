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
;    clock : in, out, optional, type=struct
;
;       Remember the 'tic/toc' clock in this variable and use it when
;       finishing.
;
;    predict : in, optional, type=boolean
;
;       Print a prediction for the remaining time. Does nothing if the
;       clock keyword is not used.
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
; 
;-
pro red_progressbar, i, N, message $
                     , barlength = barlength $
                     , nobar = nobar $
                     , clock = clock $
                     , predict = predict

  if n_elements(message) eq 0 then message = 'Progress'
  bb = string(13B)              ; CR w/o LF
  
  case N of
    0: return
    1: norm = 100.
    else: norm = 100. / (N - 1.0)  
  endcase

  percentdone = norm * i

  if arg_present(clock) then begin
    if i eq 0 or n_elements(clock) eq 0 then clock = tic()        ; Remember starting time
    time = toc(clock)
;print, time
;stop
    if keyword_set(predict) and i gt 0 then begin
      ;;prediction = ' (total estimated '+strtrim(round(time/(percentdone/100.)), 2)+' s)' else prediction = '' 
      prediction = ' ('+red_timestring(round(time*(1./(percentdone/100.)-1)), Nsecdec = 0, /interval)+' remaining)' 
    endif else prediction = '' 
    time = 'in ' + red_timestring(round(time), Nsecdec = 0, /interval) + prediction
  endif else time = ''

  if n_elements(barlength) eq 0 then barlength = 20
  if keyword_set(nobar) then begin

    print, bb, message + ' -> ' $
           , percentdone, '% ' + time + '   ', FORMAT = '(A,A,F5.1,A,$)'

  endif else begin

    time += ': '

    elength = floor(norm*i/100.*barlength)
    mlength = barlength-elength
    bar = ''
    if elength gt 0 then bar += string(replicate(61B, elength))     ; Replicated '='
    if mlength gt 0 then bar += string(replicate(45B, mlength))     ; Replicated '-'
    print, bb, '[' + bar + '] ' $
           , percentdone, '% ' + time + message + '   ', FORMAT = '(A,A,F5.1,A,$)'

  endelse

  if i eq N-1 then begin
    undefine, clock
    print                       ; Finished
  endif

end

N=50                                                                        
for i=0,N-1 do begin red_progressbar,i,N,'Test',clock=clock & wait,.1 & end    
for i=0,N-1 do begin red_progressbar,i,N,'Test' & wait,.1 & end
for i=0,N-1 do begin red_progressbar,i,N,'Test',clock=clock, /nobar& wait,.1 & end    
for i=0,N-1 do begin red_progressbar,i,N,'Test', /nobar & wait,.1 & end


for i=10,N-1 do begin red_progressbar,i,N,'Test',clock=clock, /predict & wait,.1 & end    

end
