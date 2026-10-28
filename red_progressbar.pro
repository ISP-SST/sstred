; docformat = 'rst'

;+
; Display progress in terminal window. Call near the beginning of a
; loop.
; 
; :Categories:
;
;    SSTRED
; 
; 
; :Author:
; 
;    Mats Löfdahl, ISP
; 
; 
; :Params:
; 
;    i_inner : in, optional, type="scalar number"
;
;       The current inner loop index (0 to N_inner-1).
;
;    i_outer : in, type="scalar number"
; 
;       The current (outer) loop index (0 to N_outer-1).
;
;    message_inner : in, optional, type=string, default="''"
;
;       A message identifying the inner loop task.
;
;    message_outer : in, optional, type=string, default="'Progress'"
;   
;       A message identifying the task for which progress is being
;       made. 
;
;    N_inner : in, optional, type="scalar number"
;
;       Total number of iterations in the inner loop.
; 
;    N_outer : in, type="scalar number"
;    
;       Total number of iterations in the outer loop.
;
; :Keywords:
; 
;    add_loopindex : in, optional, type=boolean
; 
;       Add the loop index value to the progress bar printout.
; 
;    barlength : in, optional, type=integer, default=20
;
;       The length of the status bar in number of characters.
;
;    done : in, optional, type=boolean
; 
;       Call with /done and (optionally) message_done after the loop
;       is finished to print the final 100% completion status.
; 
;    message_done : in, optional, type=string
; 
;       A string to go with the /done option.
; 
;    nobar : in, optional, type=boolean
;
;       If set, do not show a status bar.
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
;     2026-09-22 : MGL. Improve handling of short loops. New keywords
;                  i_inner, N_inner, message_inner
; 
;     2026-09-23 : MGL. New keywords done and message_done.
; 
;-
pro red_progressbar, i_outer, N_outer, message_outer $   ; necessary parameters
                     , i_inner, N_inner, message_inner $ ; optional parameters
                     , add_loopindex = add_loopindex $
                     , barlength = barlength $
                     , done = done $
                     , message_done = message_done $
                     , nobar = nobar $
                     , predict = predict

  common red_progressbar_common, clock, times, iterations, icnt $
     , last_i_outer, last_i_inner, last_N_outer, last_N_inner, last_message_outer

  if keyword_set(done) then begin
    i_outer = last_N_outer
    N_outer = last_N_outer
    if n_elements(last_N_inner) then begin
      i_inner = last_N_inner
      N_inner = last_N_inner
    endif
    if n_elements(message_done) gt 0 then begin
      message_outer = message_done
    endif else begin
      message_outer = last_message_outer
    endelse
  endif
  
  if N_outer eq 0 then return

  ;; 1. Initialization (Run ONLY once at the absolute start of the loop)
  if n_elements(clock) eq 0 then begin
    clock = tic()
    current_toc = 0.0
    icnt = 0L
    last_i_outer = -1L
    last_i_inner = -1L
    times = dblarr(100)
    iterations = lonarr(100)
  endif else begin
    ; If already initialized, capture the exact elapsed time in real-time
    current_toc = double(toc(clock))
  endelse

  ;; Reset if C-c earlier
  is_restart = (i_outer eq 0 && last_i_outer gt 0) || $
               (i_outer eq 0 && n_elements(i_inner) ne 0 && last_i_inner gt 0 && i_inner lt last_i_inner)
               
  if is_restart then begin
    clock = tic()
    current_toc = 0.0
    icnt = 0L
    last_i_outer = -1L
    last_i_inner = -1L
    times = dblarr(100)
    iterations = lonarr(100)
  endif
  
  if n_elements(message_outer) eq 0 then message_outer = 'Progress'
  if n_elements(message_inner) eq 0 then message_inner = ''
  bb = string(13B)              ; CR w/o LF

  norm = 100. / double(N_outer)  
  
  ;; Outer loop percentage
  percentdone = (double(i_outer) * norm) < 100.0
  prediction = '' 

  ;; Only record time history when the OUTER loop actually progresses
  ;; to a new index
  if i_outer ne last_i_outer and i_outer lt N_outer then begin
    if icnt eq n_elements(times) then begin
      times = [times, dblarr(100)]
      iterations = [iterations, lonarr(100)]
    endif
    times[icnt] = current_toc
    iterations[icnt] = i_outer
    icnt++
    last_i_outer = i_outer
    last_N_outer = N_outer
    last_message_outer = message_outer 
  endif

  ;; Track the inner loop index to detect mid-way restarts
  if n_elements(i_inner) ne 0 then begin
    last_i_inner = i_inner
    last_N_inner = N_inner
  endif

  ;; Determine how long the progress line can be
  tsize = (TERMINAL_SIZE())[0]
  outlength = (tsize - 2) > 40
  outline = string(replicate(32B, outlength))

  if n_elements(barlength) eq 0 then begin
    if n_elements(i_inner) eq 0 then begin
      ;; No nested loop
      barlength = 20
    endif else begin
      ;; Nested loop
      barlength = 10
    endelse
  endif
  
  ;; 2. Predict remaining time based on outer loop history
  if keyword_set(predict) and icnt gt 1 and i_outer lt N_outer then begin
    diff_times = times[1:icnt-1] - times[0:icnt-2]
    diff_iters = iterations[1:icnt-1] - iterations[0:icnt-2]
    valid = where(diff_iters gt 0, n_valid)
    
    if n_valid gt 0 then begin
      dt = median(diff_times[valid] / double(diff_iters[valid]))
      time_remaining = (double(N_outer) - double(i_outer)) * dt
      if time_remaining lt 0 then time_remaining = 0.0
      prediction = ' (' + strtrim(red_timestring(double(round(time_remaining)), Nsecdec = 0, /interval), 2) + ' rem)'
    endif
  endif

  elapsed_seconds = double(round(current_toc))
  time_str = 'Time: ' + strtrim(red_timestring(elapsed_seconds, Nsecdec = 0, /interval), 2) + prediction + ' | '

  ;; 3. Construct the Outer Loop bar
  outer_str = ''
  pct_outer_str = string(percentdone, FORMAT='(F5.1)') + '%'
  if keyword_set(nobar) then begin
    outer_str = 'Outer: ' + pct_outer_str + ' '
  end else begin
    elength = floor(percentdone/100*barlength)
    mlength = barlength-elength
    bar = ''
    if elength gt 0 then bar += string(replicate(61B, elength)) ; '='
    if mlength gt 0 then bar += string(replicate(45B, mlength)) ; '-'
    outer_str = '[' + bar + '] ' + pct_outer_str + ' '
  endelse

  ;; 4. Construct the Inner Loop bar
  nested_str = ''
  if n_elements(i_inner) ne 0 && n_elements(N_inner) ne 0 then begin
    norm_inner = 100. / (double(N_inner) > 1.0)
    pct_inner = (double(i_inner) * norm_inner) < 100.0
    pct_inner_str = string(pct_inner, FORMAT='(F5.1)') + '%'
    
    if keyword_set(nobar) then begin
      nested_str = '| Inner: ' + pct_inner_str + ' '
    end else begin
      e_inner = floor(pct_inner/100 * barlength)
      m_inner = barlength - e_inner
      bar_inner = ''
      if e_inner gt 0 then bar_inner += string(replicate(61B, e_inner))
      if m_inner gt 0 then bar_inner += string(replicate(45B, m_inner))
      nested_str = '| [' + bar_inner + '] ' + pct_inner_str + ' '
    endelse
  endif

  ;; 5. Construct the messages tail
  msg_str = '| ' + strtrim(message_outer, 2)
  if keyword_set(add_loopindex) then msg_str += ' '+strtrim(i_outer, 2)
  if strtrim(message_inner, 2) ne '' then msg_str += ' -> ' + strtrim(message_inner, 2)
  if keyword_set(add_loopindex) then msg_str += ' '+strtrim(i_inner, 2)

  ;; Combine and truncate safely
  full_line = time_str + outer_str + nested_str + msg_str
  if strlen(full_line) gt outlength then full_line = strmid(full_line, 0, outlength)
  
  strput, outline, full_line

  print, bb, outline, FORMAT = '(A,A,$)'
  
  ;; Final cleanup when outer loop hits completion target
  if i_outer ge N_outer then begin
    undefine, clock
    undefine, times
    undefine, iterations
    undefine, icnt
    undefine, last_i_outer
    undefine, last_i_inner
    undefine, last_N_outer
    undefine, last_N_inner
    undefine, last_message_outer
    print
  endif

end

N = 3
M = 4
messages = ['short.fits', 'very_long_filename_record_2026.fits', 'med.fits']

for i = 0, N-1 do begin
  red_progressbar, i, N, messages[i mod N]
  wait, 1
endfor
red_progressbar, /done, message_done = 'files'
;red_progressbar, N, N, 'Done'

for i = 0, N-1 do begin 
  for j = 0, M-1 do begin 
    red_progressbar, i, N, messages[i], j, M, 'SubStep', /predict, /add
     wait, 1
  endfor 
endfor
;red_progressbar, N, N, 'Done', M, M, 'SubSteps'
red_progressbar, /done, message_done = 'files'

stop

N = 3
M = 4
messages = ['short.fits', 'very_long_filename_record_2026.fits', 'med.fits']
for i = 0, N-1 do begin 
  for j = 0, M-1 do begin 
    red_progressbar, i, N, messages[i], j, M, 'SubStep_' + strtrim(j,2), /predict 
     wait, 2.4 
  endfor 
endfor
red_progressbar, N, N, 'Done', M, M, 'SubSteps'


stop

N = 3
M = 4
messages = ['short.fits', 'very_long_filename_record_2026.fits', 'med.fits']
for i = 0, N-1 do begin 
  for j = 0, M-1 do begin 
    red_progressbar, i, N, messages[i], j, M, 'Sub', /predict 
    wait, 0.4 
  endfor 
endfor
red_progressbar, N, N, 'Done'

stop
  


N = 3
M = 5
for i = 0, N-1 do begin 
  for j = 0, M-1 do begin 
    red_progressbar, i, N, 'Nested Test', j, M, /predict
    wait, 1.0 
  endfor 
endfor
red_progressbar, N, N, 'Nested Test'


N = 3
for i = 0, N-1 do begin & red_progressbar,i,N,'Test', /predict & wait,5.31 & endfor
red_progressbar,N,N,'Test'


stop

N=50000
;for i=0,N-1 do begin red_progressbar,i,N,'Test',clock=clock & wait,.1 & end    
;for i=0,N-1 do begin red_progressbar,i,N,'Test' & wait,.1 & end
;for i=0,N-1 do begin red_progressbar,i,N,'Test',clock=clock, /nobar& wait,.1 & end    
;for i=0,N-1 do begin red_progressbar,i,N,'Test', /nobar & wait,.1 & end


for i=0,N-1 do begin & if i mod 100 eq 0 or i eq N-1 then red_progressbar,i,N,'Test', /predict & wait,.001 & end  
red_progressbar,N,N,'Test'  


for i=0,N-1 do begin red_progressbar,i,N,'Test', /predict & wait,.001 & end    
;for i=0,N-1 do begin red_progressbar,i,N,'Test',clock=clock, /predict & wait,.1 & end    

end
