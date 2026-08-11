; docformat = 'rst'

;+
; Print out some statistics for an array.
;
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, Institute for Solar Physics, 1990 (Original ANA version.)
; 
; 
;
; :Params:
;
;    x : in, type=array
;
;      The array for which statistics is to be printed.
;
;
; :Keywords:
; 
;    name : in, optional, type=string
; 
;       The name of the variable or some other label that identifies
;       the statistics.
; 
;    nan : in, optional, type=boolean
;
;       Passed when calling min(), max(), etc.
;
;
; :History:
; 
;      2011 : MGL. Ported from ANA.
;
;      2013-09-11 : MGL. Renamed red_stats for inclusion in crispred
;                   pipeline. 
;
;      2013-09-12 : MGL. Substituted "print" for ANA-ism "type". 
;
;      2021-12-01 : MGL. New keywords out and heading.
;
;      2022-06-07 : MGL. New keyword nan.
;
;-
pro red_stats, x, name=name, out = out, heading = heading, nan = nan ; Show statistics about x

  DivisorString = '------------------------------------------------------------------'

  MaxString  = strtrim(Max(float(x), nan = nan),2)
  MinString  = strtrim(Min(float(x), nan = nan),2)
  MeanString = strtrim(Mean(float(x), nan = nan),2)
  DevString  = strtrim(StDdev(float(x), nan = nan),2)

  heading = 'Min, Mean, Max, StDev:'
  
  if arg_present(out) then begin
    out = MinString+' '+MeanString+' '+MaxString+' '+DevString
    return
  endif
  
  IF n_elements(name) gt 0 THEN BEGIN
     print,DivisorString
     print,'Min, Mean, Max, StDev of '+name+':'  
  end else begin
     print,DivisorString
     print,'Min, Mean, Max, StDev:'  
  end
  print,MinString+' '+MeanString+' '+MaxString+' '+DevString
  print,DivisorString

  return

; Port the rest of the ana subroutine to get more nicely formatted output:

  MaxBlanks  = strlen(MaxString)  - 3
  MinBlanks  = strlen(MinString)  - 3
  MeanBlanks = strlen(MeanString) - 4 
  DevBlanks  = strlen(DevString)  - 4
print,maxblanks,minblanks,meanblanks,devblanks
  DataString  = '| ' + MinString  + Blanks(-MinBlanks)
  DataString += ' | ' + MeanString + Blanks(-MeanBlanks)
  DataString += ' | ' + MaxString  + Blanks(-MaxBlanks)
  DataString += ' | ' + DevString  + Blanks(-DevBlanks) + ' |'
print,datastring  
  TitleString  = '| Min'  + Blanks(MaxBlanks)
  TitleString += ' | Mean' + Blanks(MeanBlanks)
  TitleString += ' | Max'  + Blanks(MinBlanks)
  TitleString += ' | SDev' + Blanks(DevBlanks) + ' |' 
  
  DataLength = strlen(DataString)
  
  DivisorString = '------------------------------------------------------------------'

  IF n_elements(name) gt 0 THEN BEGIN
    NameString = '| Statistics for '+name+':'
    NameLength = strlen(NameString)
    NameString = NameString + Blanks(DataLength-NameLength-1) + '|'
    NameLength = strlen(NameString)
    print, strmid(DivisorString,0,NameLength-1)
    print,NameString
    print, strmid(DivisorString,0,NameLength-1)
  END ELSE print, strmid(DivisorString,0,DataLength-1)
  
  print, TitleString
  print, DataString
  print, strmid(DivisorString,0,DataLength-1)


END; Stats
