;+
; Open a mySQL database for operations via a pipe.
; 
; Assumes that your .my.cnf file points to the correct mysql server
; and that it includes the login information. Also, the command
; 'mysql' must appear in your default path.
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
;    Based on openmysql.pro by Marc W. Buie, Lowell Observatory
; 
; :Returns:
; 
; 
; :Params:
; 
;     dbname : in, type=string
;
;        Name of database to open at start.
;
;     error : out, type=integer
;
;        Return value indicating if the open call succeeded. If error
;        is zero, the open was good and the lun is ready for use. If
;        error is not zero, the lun will not point to an open file and
;        the !error_state system variable will contain more
;        information if you want it.
; 
;     lun : out, integer
;
;        The logical unit of the pipe (use free_lun to close).
; 
; :Keywords:
; 
;     host : out, optional, type=string
;
;           Host name of the server to connect to for queries.  The default
;           is to use the host name specified in your .my.conf configuration
;           file.
;
;     user  : out, optional, type=string
;
;           User name to be used for the connection.  This is the user name
;           as understood by mysql and has no relationship to the user name
;           of the calling process.  The default is to use the user name
;           specified in your .my.conf configuration file.
;
;    socket : out, optional, type=string
;
;           String, identifies socket for client to read and write from.
;           default is to not use a socket.
;   
;     PW : in, optional, type=string
;
;           The mysql password.
; 
; :History:
; 
;  2002-01-09 : MWB. First version.
;
;  2006-12-08 : MWB. Added error output argument, now can indicate if
;               open failed
;
;  2007-09-06 : MWB. Added HOST/USER keywords
;
;  2010-03-05 : MWB. Added the EXTRACONFIG keyword
;
;  2011-11-23 : Jarle Brinchmann, Leiden Univ., Add SOCKET keyword to
;               support reading from sockets
; 
;  2018-06-28 : MGL. Renamed from original "openmysql" to
;               "red_mysql_open" for inclusion in SST data pipelines.
;               Adaptions to SST style. New keyword pw.
; 
;-
pro red_mysql_open, lun, dbname, error $
                    , host=host $
                    , user=user $
                    , pw = pw $
                    , extraconfig=extraconfig $
                    , socket=socket
  
  inam = red_subprogram(/low, calling = inam1)

  if red_badpar(lun,[0,2,3],0,caller=inam+' : (lun) ') then return
  if red_badpar(dbname,7,0,caller=inam+' : (dbname) ') then return
  if red_badpar(host,[0,7],0,caller=inam+' : (host) ',default='') then return
  if red_badpar(user,[0,7],0,caller=inam+' : (user) ',default='') then return
  if red_badpar(pw,[0,7],0,caller=inam+' : (pw) ',default='') then return
  if red_badpar(extraconfig,[0,7],0,caller=inam+' : (extraconfig) ', $
                default='') then return
  if (red_badpar(socket,[0,7],0,caller=inam+' : (socket)', $
                 default='')) then return

  cmd = 'mysql '
  if host ne '' then cmd += '-h '+host+' '
  if user ne '' then cmd += '-u '+user+' '
  if pw ne '' then cmd += '-p'+pw+' '
  if extraconfig ne '' then cmd += '--defaults-extra-file='+extraconfig+' '
  if socket ne '' then cmd += '--socket='+socket+' '
  cmd += '-B -q -n '
  cmd += dbname
  spawn,cmd,unit=lun

  catch,error

  if error ne 0 then begin
    print,inam+' : ', !error_state.msg
    print,inam+' : Unable to open database '+dbname+', closing pipe and quitting.'
    free_lun,lun
    catch,/cancel
    return
  endif
  
  ;; This command is needed to trigger the catch and facilitate
  ;; closing of the lun. We don't need these results, just the
  ;; activity on the pipe.
  red_mysql_cmd, lun, 'show tables;', result, nlines

end
