; docformat = 'rst'

;+
; Figure out the current site and default data directories and
; similar.
;
; Old 'hostname -I' doesn't work on La Palma under an old OpenSuSE
; installation therefore the only supported keyword there is -i that gives
; the right IP-address on La Palma but not in Stockholm.  There are more
; robust but less trivial ways to get the IPv4-address of the machine.  As
; for 31 July 2017, the following works:
;   $ ip addr show | grep -Po 'inet \K[\d.]+'
; or
;   $ /sbin/ifconfig -a | grep -Po 't addr:\K[\d.]+'
; Both produce a list of IPv4-addresses for each network interface on La
; Palma as well as in Stockholm.  You can narrow the list to only one
; address by specifying 'ip addr show eth0' or '/sbin/ifconfig eth0' but it
; is better to use full listing as some interfaces might be temporarily off
; the network.
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
;   dnsdomainname :  out, optional, type=string
;
;      The output of shell commend dnsdomainname.
;
;   ipv4addresses : out, optional, type=strarr
;
;      The found ipv4 addresses.
;
;   search_dirs : out, optional, type=strarr
;
;      Where to search for raw data.  Could be an array or a single
;      string. 
;
;   site :  out, optional, type=string
;
;      A string identifying the current site, like "AlbaNova", "La
;      Palma", etc.
;
; :History:
;
;    2017-08-15 : MGL. First version based on logic by Andrii Sukhorukov.
;
;    2017-08-29 : MGL. Use dnsdomainname to match also nodes that do
;                 not have world-wide ipv4 addresses. New keywords
;                 ipv4addresses and dnsdomainname.
;
;-
pro red_currentsite, site = site $
                     , search_dirs = search_dirs $
                     , ipv4addresses = ipv4addresses $
                     , dnsdomainname = dnsdomainname

  ;; Name of this subroutine
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  spawn, "ip addr show | grep -Po 'inet \K[\d.]+'", ipv4addresses
  spawn, 'dnsdomainname', dnsdomainname
  
  ;; Check if any of the addresses matches the local range of
  ;; IPv4-addresses.
  case 1 of
    strmatch(dnsdomainname, '*royac.iac.es' ) $
       || max( strmatch( ipv4addresses, '*161.72.15.*' ) ) : begin
      ;; SST network address range on La Palma.
      site = "La Palma"
      ;; Observed data is in any folder mounted at /data/, either disk? or
      ;; camera?, or some other drive.  As Pit said on 28 July 2017, we should
      ;; search /data/camera? folders only if it is really necessary.  Observed
      ;; data must be first transfered to /data/disk? folders, by default to
      ;; disk1.  The inner directory level is either CHROMIS folder for CHROMIS
      ;; data or something different, usually the name of the insitute like
      ;; ISP, UK, UIO, for CRISP data.  For example,
      ;;   /data/disk1/CHROMIS/2017.04.05
      ;;   /data/disk1/ISP/2017.04.05
      ;;search_dirs = "/data/*/*/" ; Full contents of /data folder.
      search_dirs = "/data/disk?/*/" ; Only /data/disk? storages.
    end
    max( strmatch( ipv4addresses, '*130.237.166.*' ) ) : begin
      ;; ISP network range at AlbaNova in Stockholm.
      site =  "AlbaNova"
      ;; Transfered data is stored in all the sandboxes at /storage/sand*.
      ;; The inner directory level is either empty, or Incoming, or
      ;; Incoming/Checked.  For example:
      ;;   /storage/sand04n/2017.04.05
      ;;   /storage/sand05/Incoming/2017.04.05
      search_dirs = '/storage/sand*/' + ['', 'Incoming/', 'Incoming/Checked/']
    end
    else : begin
      message, 'No matching IPv4-address in ' + strjoin( ipv4addresses, ', ' )
      retall
    end
  endcase

  print, inam + ' : We are in ', site, '.'

end
