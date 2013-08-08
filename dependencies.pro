;; Use this program to locate subprograms that are not in the local
;; directory, not in the IDL distribution, not in astrolib, and not in
;; the coyote library.
;;
;; First compile the main program and then run the program:
;;
;; IDL> .comp crispred
;; IDL> .run dependencies
;;
;; Two output files are generated. 
;;
;;  * dependencies.out : lists the paths of the subprograms that are found.
;;
;;  * dependencies.grep : lists the local files where the found
;;                        subprograms are called.
;;
;; Before making local copies of the found subprograms and changing
;; their names to the red_* namespace, make sure this isn't
;; already done so you only have to change how they are called.
;;
;; Finding everything is a recursive process. Run again after fixing
;; what was found.
;;
;; Mats Löfdahl 2013-07-24.

pro prnt, text, lun

  print, text
  printf, lun, text

end

pwdstring = getenv('PWD')+'*'
idlstring = '*common/pkg/itt/idl*'
coyotestring = '*Fanning/coyote*'
astrolibstring = '*astronomy_users_library*'

okstrings = [pwdstring, idlstring, coyotestring, astrolibstring]

spawn, 'rm -f dependencies.grep'
openw, lun, 'dependencies.out', /get_lun

Nnames = 0
Nsubr = 0
Nfunc = 0
repeat begin

   oldNnames = Nnames

   ;; Compile all subprograms found so far
   if Nsubr ne 0 then resolve_routine,subrnames,is_function=0
   if Nfunc ne 0 then resolve_routine,funcnames,is_function=1
   
   ;; Find
   subrnames = strlowcase(routine_info(/unresolved,functions=0))
   funcnames = strlowcase(routine_info(/unresolved,functions=1))
   Nsubr = (size(subrnames, /dim))[0]
   Nfunc = (size(funcnames, /dim))[0]
 
   ;; Concatenate lists
   if Nsubr eq 0 and Nfunc eq 0 then begin
      Nnames = 0
      names = ''
   endif else begin
      if Nsubr gt 0 and Nfunc eq 0 then begin
         names = [subrnames, funcnames]
      endif else begin
         if Nsubr gt 0 then names = subrnames else names = funcnames
      endelse
      Nnames = (size(names, /dim))[0] 
   endelse
   
   if Nnames gt 0 then begin

      prnt, '=============== Found '+strtrim(string(Nnames), 2)+' subprograms!', lun
      prnt, '', lun

      if Nsubr gt 0 then subrinpath = bytarr(Nsubr)+1 else subrinpath = 0
      if Nfunc gt 0 then funcinpath = bytarr(Nfunc)+1 else funcinpath = 0
      
      for i = 0, Nnames-1 do begin
         findpro, names[i], /NoPrint, dirlist = dirlist ;, prolist = prolist
         ; Could also use 
         ;   pth=routine_info(names[i],/source)
         ; but then we don't get info about name collisions.

         if dirlist[0] eq '' then begin
            if i lt Nsubr then begin
               prnt, 'Subroutine not in path: '+names[i], lun
               subrinpath[i] = 0
            endif else begin
               prnt, 'Function not in path: '+names[i], lun
               funcinpath[i-Nsubr] = 0
               printblankline = 1
            endelse
         endif else begin
            Ndirs = (size(dirlist, /dim))[0]
            printblankline = 0
            for j = 0, Ndirs-1 do begin
               isOK = 0 ; OK if one of the okstrings matches the directory
               for ii = 0, (size(okstrings, /dim))[0]-1 do begin
;                  print, strmatch(dirlist[j], okstrings[ii], /fold)
                  isOK = isOK or strmatch(dirlist[j], okstrings[ii], /fold)
               endfor
;               islocal = strmatch(dirlist[j], pwdstring, /fold) 
;               isidl = strmatch(dirlist[j], idlstring, /fold)
               if ~isOK then begin
                  prnt, dirlist[j]+names[i], lun
                  spawn, 'echo "----- "'+names[i]+'" -----" >> dependencies.grep'
                  spawn, 'grep -i '+names[i]+' *.pro >> dependencies.grep'
                  printblankline = 1
               endif
            endfor
         endelse
         if printblankline then prnt, '', lun
      endfor
      

      if total(subrinpath) gt 0 then begin
         subrnames = subrnames[where(subrinpath)]
         Nsubr = (size(subrnames, /dim))[0] 
      endif else begin
         subrnames = ''
         Nsubr = 0
      endelse

      if total(funcinpath) gt 0 then begin
         funcnames = funcnames[where(funcinpath)]
         Nfunc = (size(funcnames, /dim))[0] 
      endif else begin
         funcnames = ''
         Nfunc = 0
      endelse

      Nnames = Nsubr+Nfunc

   endif

endrep until Nnames eq 0

free_lun, lun

end

