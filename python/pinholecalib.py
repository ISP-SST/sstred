#!/usr/bin/python

# Mats Lofdahl 2009-10-26

# Changes by Pit:
#    Catch 'no config files supplied' error
#    Replace call to imgadd with pyana
#    Replace some os.system calls with their os.method
#    implement plane fit in ana
#    create usable median filter

from numpy import *
from scipy.signal.signaltools import medfilt2d, medfilt
import pyana

import time
import sys
import os
import commands
import subprocess
import string
#import fnmatch    # not used?
import re
from optparse import OptionParser


# Kill subprocesses without the version 2.6 .kill() method.
def killsubprocess(sp):
    #sp.kill() Need python 2.6 for this
    #pid = sp.pid
    #print 'Kill process ',pid
    #os.system('kill -9 '+str(pid))
    os.kill(sp.pid, 9)

## Get tag value from config file after optionally skipping past strings in <skip>:
def getcfgtag(cfglines,tag,options,skip=None):
    if options.verbose: print 'getcfgtag'
    i=-1
    # First the skipping
    if skip:
        for skipstring in skip:
            if options.verbose: print "-----Skip past ",skipstring
            for i in range(i+1,len(cfglines)):
                #if options.verbose: print cfglines[i].rstrip("\n")
                if re.search(skipstring,cfglines[i]):
                    #if options.verbose: print "-----Found ",skipstring
                    break
    # Then the location of the tag
    for i in range(i+1,len(cfglines)):
        #if options.verbose: print cfglines[i].rstrip("\n")
        if re.search(tag,cfglines[i]):
            return cfglines[i].split('=')[1].rstrip("\n")
    else:
        print "No such tag "+tag+" in "
        print cfglines
        return None

## Write cfglines to cfg_file. Optionally replace the amount of
## diversity given for the non-anchor channel. Optionally replace the
## modes. 
def writeconfig(cfglines,cfg_file,options,diversity=None,modes=None):
    if options.verbose: print "Write config file: ",cfg_file
    #print "with contents: ",cfglines
    if modes: 
        if options.verbose: print 'Modes: ',modes
        for i in range(len(cfglines)):
            if re.search('MODES=',cfglines[i]):
                if options.verbose: print 'Change: ',cfglines[i].rstrip("\n")
                cfglines[i] = 'MODES='+modes+'\n'
                if options.verbose: print '  to: ',cfglines[i].rstrip("\n")
                break
        else:
            cfglines.append('MODES='+modes+'\n')
            if options.verbose: print 'Add line: ', cfglines[-1].rstrip("\n")
    if diversity:
        if options.verbose: print 'Diversity: ',diversity
        i=0
        # Skip past beginning of second channel def:
        while not re.search('channel',cfglines[i]): i += 1
        i += 1
        while not re.search('channel',cfglines[i]): i += 1
        i += 1
        # Now search for the DIVERSITY tag:
        for j in range(i,len(cfglines)):
            if re.search('DIVERSITY=',cfglines[j]):
                if options.verbose: print 'Change: ',cfglines[j].rstrip("\n")
                cfglines[j] = '    DIVERSITY='+diversity+'\n'
                if options.verbose: print '  to: ',cfglines[j].rstrip("\n")
                break
        else:
            sys.exit('No DIVERSITY tag in second channel def of '+cfg_file)

    f = open(cfg_file,'w')
    f.writelines(cfglines)
    f.close()

    return cfglines

def jstat(jORs,options):
    cmdarr = ["jstat","-p",str(options.momfbd_port),"-"+jORs]
    p = subprocess.Popen(cmdarr, stdout=subprocess.PIPE)
    if options.verbose: print " ".join(cmdarr)
    output = p.stdout.readlines()
    if len(output) > 1: print output[1].rstrip("\n")
    return output

## Read alpha file:
def readalpha(aname):
    f = open(aname,'r')
    alines = f.readlines()
    f.close()
    
    xd,yd,ch,sx,sy = [long(a) for a in alines[0].split()]
    xy = zeros((xd,yd,2),int32)
    alpha = zeros((xd,yd,ch))
    
    n=0
    for i in range(xd):
        for j in range(yd):
            n+=1
            linearray = alines[n].split()
            for k in range(2): xy[i,j,k] = long(linearray[k])
            for k in range(ch): alpha[i,j,k] = float(linearray[2+k])
    
    return xy,alpha,sx,sy

def median_2d(data, kernel_size=3):
    
    if data.ndim != 2:
        sys.exit("Only for 2D data!")
    
    sx,sy = shape(data)
    k2 = int(kernel_size)/2
    if 2*k2 == kernel_size:
        sys.exit("Kernel size must be odd!")
    
    result = copy(data)
    for k in range(1,k2+1):
        result[k:sx-k, k:sy-k] = medfilt2d(data, kernel_size=2*k+1)[k:sx-k, k:sy-k]
    return result


def fitplane(afile1, afile2, xfile1, xfile2, options):

    xy,alpha1,sx,sy = readalpha(afile1)
    xy,alpha2,sx,sy = readalpha(afile2)
    alpha2 -= alpha1

    alpha2 *= options.gain

    if options.removespikes:
        for m in range(ch):
            # that's what you would like, but medfilt2d returns bad data
            # for the edge (i.e., zero or random).  Use a wrapper
            # alpha2[:,:,m]= medfilt2d(alpha2[:,:,m], kernel_size=ks)
            alpha2[:,:,m] = median_2D(alpha2[:,:,m], kernel_size=options.median)

    clip = options.clip

    xd,yd,ch = shape(alpha2)
    nr = (xd-2*clip)*(yd-2*clip)
    x = xy[clip:xd-clip,clip:yd-clip,0].reshape(nr)
    y = xy[clip:xd-clip,clip:yd-clip,1].reshape(nr)
    shx = alpha2[clip:xd-clip,clip:yd-clip,0].reshape(nr)
    shy = alpha2[clip:xd-clip,clip:yd-clip,1].reshape(nr)
    if ch == 3:
        defoc = alpha2[clip:xd-clip, clip:yd-clip, 2].mean()
    else:
        defoc = 0.0

    ut = zeros((4, nr), float64)
    k = 0
    #for i in range(2):
    #    for j in range(2):
    #        ut[k,:] = x**i * y**j
    #        k += 1
    ut[0,:] = ones(nr)
    ut[1,:] = y
    ut[2,:] = x
    ut[3,:] = x*y

    kk = inner(linalg.inv(inner(ut,ut)),ut.transpose())
    kxx = inner(kk,shx)
    kxy = inner(kk,shy)

    X,Y = mgrid[0:sy, 0:sx]
    fitx = int16((100 * (kxx[1]*X + kxx[2]*Y + kxx[3]*Y*X + kxx[0])).round())
    fity = int16((100 * (kxy[1]*X + kxy[2]*Y + kxy[3]*Y*X + kxy[0])).round())
    mm = max(abs(array([fitx.min(),fitx.max(),fity.min(),fity.max()])))

    pyana.writeto(xfile1, pyana.getdata(xfile1)+fitx)
    pyana.writeto(xfile2, pyana.getdata(xfile2)+fity)
#    pyana.writeto(xfile1, pyana.getdata(xfile1)+fitx, comments=' ')
#    pyana.writeto(xfile2, pyana.getdata(xfile2)+fity, comments=' ')

    return mm,defoc


def iterate(cfglines,cfg_file,options,tiltprogress,focusprogress,modes='2-3',maxiter=None):

    ## Get some info from cfg file:  
    xoffsetfile = getcfgtag(cfglines,'XOFFSET',options)
    if not xoffsetfile : sys.exit('No XOFFSET specified')

    yoffsetfile = getcfgtag(cfglines,'YOFFSET',options)
    if not yoffsetfile: sys.exit('No YOFFSET specified')

    channel_number = getcfgtag(cfglines,'IMAGE_NUM',options)
    if not channel_number: sys.exit('No IMAGE_NUM specified')


    # Make base part of output file names:
    filename_template = getcfgtag(cfglines,'FILENAME_TEMPLATE',options) 
    # Replace place holder for image number "channel_number..channel_number".
    p=re.compile('%0[1-9]d')
    fname_base = p.sub("".join([channel_number,"..",channel_number]),filename_template)

    cfglines = writeconfig(cfglines,cfg_file,options,modes=modes)

    ii=0
    done=False
    while not done:
        ii += 1


        # Output file names:
        outfile=fname_base+".f0"
        alphafile1=fname_base+".alpha.1.1"
        alphafile2=fname_base+".alpha.2.1"
        logfile = fname_base+'.log'

        # Make sure the output files do not exist already
        os.system('rm -f '+outfile)
        os.system('rm -f '+alphafile1)
        os.system('rm -f '+alphafile2)

        # Submit the MOMFDB job:
        if modes == '2-3':
            name = 'Tilts:'
        else:
            name = 'Diversity:'
        name += cfg_file+":it"+str(ii)

        # jsub command from calibration.sh:
        # jsub -p 7171 -pri -2 -v -cfg momfbd.calib.98.cfg -f -name 98_it2 -lg camXX_im17Jul2008.98..98.pinh.log
        JSUB=["jsub","-p",str(options.momfbd_port)]
        cmd = " ".join(JSUB+["-pri","-2","-v","-cfg",cfg_file,"-f","-name",name,"-lg",logfile])
        os.system(cmd)

        if options.verbose: print cmd        

        while not ( os.path.isfile(outfile) and os.path.isfile(alphafile1) and os.path.isfile(alphafile2) ):
            time.sleep(float(options.tsleep))
            a=jstat("j",options)


        maxoffset, focuserror = fitplane(alphafile1, alphafile2,
                                         xoffsetfile, yoffsetfile,options)
        tiltprogress.append(float(maxoffset))
        focusprogress.append(float(focuserror))

        if modes == '2-4':
            print ''
            # Update the amount of diversity in cfglines
            # new diversity = old diversity minus the defocus error from this iteration.
            olddiversity = float(getcfgtag(cfglines,'DIVERSITY',options,skip=['channel','channel']).split()[0])
            newdiversity = olddiversity - focuserror
            if options.verbose: print 'Focus (old,delta,new): ',olddiversity,focuserror,newdiversity
            
            cfglines = writeconfig(cfglines,cfg_file,options,diversity=str(newdiversity)+' mm')

        
        print 'Tilt convergence: ',tiltprogress
        print 'Focus convergence: ',focusprogress

        done = ((maxoffset <= options.shiftcriterion) and (focuserror <= options.focuscriterion))
        done = done or (ii >= maxiter) 

    return cfglines,tiltprogress,focusprogress


## Find a free port to use for MOMFBD, need this before the option parsing
def freeport():
    # Get list of used ports from netstat:
    #for line in commands.getoutput('netstat -atn').splitlines()[2:]: print line.split()[3].split(':')[-1]
    usedports = [long(line.split()[3].split(':')[-1]) for line in commands.getoutput('netstat -atn').splitlines()[2:]]
    #print usedports
    # Step down from 32768 until we find an unused port:
    for port in xrange(32768,10000,-1):
        if port not in usedports:
            return port
        else:
            port += -1
    else:
        sys.exit('No free ports')

def main():

    momfbd_port_default = freeport()

    ## Parse command line options
    parser = OptionParser(usage="usage: %prog [options] cfg-file [cfg_file ...]")

    parser.add_option("-d", "--diversity", dest="diversity",
                      default=False, action="store_true", 
                      help="Calibrate diversity (not implemented yet)")

    parser.add_option("-s", "--slaves", type="int" ,dest="Nslaves",
                      default=1, metavar="#", 
                      help="Number of slaves to start [default=%default]")

    parser.add_option("-p", "--port", dest="momfbd_port",
                      default=momfbd_port_default,metavar="#",
                      help="MOMFBD port, defaults to a free port")

    parser.add_option("-T", "--Tsleep", type="int", dest="tsleep", 
                      default=3,metavar="#",
                      help="""Time to sleep between jstats in seconds [default=%default]""")

    parser.add_option("-i", "--iterations", type="int", dest="maxiterations", 
                      default=10, metavar="#",
                      help="Max number of iterations [default=%default]")

    parser.add_option("-S", "--shiftcrit", type="int", dest="shiftcriterion", 
                      default=10, metavar="#",
                      help="Shift criterion for stopping iterations. In hundedths of a pixel [default: #=%default]")

    parser.add_option("-F", "--focuscrit", type="float", dest="focuscriterion", 
                      default=0.1, metavar="#",
                      help="Focus criterion for stopping iterations. In mm. [default: #=%default]")

    parser.add_option("-a", "--align", dest="align", default=False,
                      action="store_true", 
                      help="""Align by cross correlation to get starting offsets (not implemented yet)""") 

    parser.add_option("-r", "--removespikes", dest="removespikes",
                      default=False, action="store_true", 
                      help="Remove spikes from estimated quantities by median filtering")

    parser.add_option("-m", "--median", type="int", dest="median",
                      default=3, metavar="#",
                      help="Kernel size for remove spikes [default: #=%default]")

    parser.add_option("-c", "--clip", type="int", dest="clip", 
                      default=1, metavar="#",
                      help="Borderclip for plane fit to shifts. [default: #=%default]")

    parser.add_option("-v", "--verbose", action="store_true",
                      dest="verbose", default=False, 
                      help="Print status messages to stdout") 

    parser.add_option("-t", "--test", action="store_true",
                      dest="test", default=False, 
                      help="Run test() indstead of main()") 

    parser.add_option("-g", "--gain", type="float", dest="gain", 
                      default=1.0, 
                      help="Gain applied to measured Z coeffs (silly debug thingie) [default=%default]")

    (options, cfg_files) = parser.parse_args()

    if len(cfg_files) == 0:
        parser.error('No config files specified')

    if options.verbose: 
        print cfg_files
        if options.diversity:
            print 'Do diversity processing'
        print "MOMFBD port: ",options.momfbd_port

    ## Make data subdir
    if not os.path.isdir('data'):
        os.mkdir('data')

    # Start MOMFBD manager
    MANAGER = ["manager","-v","-p",str(options.momfbd_port)]
    if options.verbose: print " ".join(MANAGER)
    manager = subprocess.Popen(MANAGER)

    # Start MOMFBD slaves
    SLAVE=["momfbd_slave","-p",str(options.momfbd_port)]
    slaves = []
    for i in range(options.Nslaves):
        if options.verbose: print " ".join(SLAVE)
        slaves.append(subprocess.Popen(SLAVE))

    #if options.verbose: print " ".join(JSTATJ)
    #jstatj = subprocess.Popen(JSTATJ, stdout=subprocess.PIPE)
    #for line in jstatj.stdout.readlines():
    #    if options.verbose: print line.rstrip("\n")
    #
    #if options.verbose: print " ".join(JSTATS)
    #jstats = subprocess.Popen(JSTATS, stdout=subprocess.PIPE)
    #for line in jstats.stdout.readlines():
    #    if options.verbose: print line.rstrip("\n")



    for cfg_file in cfg_files:
        print "--------------------------------------------------Config file: ",cfg_file

        ## Read the config file. We will keep the current version in
        ## memory during the process.
        f = open(cfg_file,'r')
        cfglines = f.readlines()
        f.close()

        ## Write backup of original config file
        writeconfig(cfglines,cfg_file+'.bak',options)

        # Make sure the config file has DATA_TYPE=FLOAT
        # Without it, you can run into overflow problems.
        for cfgline in cfglines:
            if cfgline.rstrip("\n") == 'DATA_TYPE=FLOAT':
                break
        else:
            cfglines.append('DATA_TYPE=FLOAT\n')
            writeconfig(cfglines,cfg_file,options)

        # Initialize progress logs
        tiltprogress = []
        focusprogress = []


        if options.diversity:

            print 'Iterate tilts a couple of times'
            cfglines,tiltprogress,focusprogress = iterate(cfglines, cfg_file, options, tiltprogress,
                                                  focusprogress,modes='2-3',maxiter=min((2,options.maxiterations)))

            print 'Iterate until diversity is OK'
            cfglines,tiltprogress,focusprogress = iterate(cfglines, cfg_file, options, tiltprogress, 
                                                  focusprogress,modes='2-4', maxiter=options.maxiterations)

        print 'Iterate until tilts are OK'
        cfglines,tiltprogress,focusprogress = iterate(cfglines, cfg_file, options, tiltprogress, 
                                                      focusprogress,modes='2-3', maxiter=options.maxiterations)

        # Write progress logs to disk
        print cfg_file.replace('.cfg','.progress')
        f = open(cfg_file+'.progress','w')
        f.write("%10s %10s\n" % ('tilt','focus'))
        for k in range(len(tiltprogress)):
            f.write("%10d %10d\n" % (tiltprogress[k],focusprogress[k]))
        f.close()


    ## Kill the slaves and the manager
    for i in range(options.Nslaves):
        killsubprocess(slaves[i])
        killsubprocess(manager)

def test():
    momfbd_port_default = freeport()

    ## Parse command line options
    parser = OptionParser(usage="usage: %prog [options] cfg-file [cfg_file ...]")

    parser.add_option("-d", "--diversity", dest="diversity",
                      default=False, action="store_true", 
                      help="Calibrate diversity (not implemented yet)")

    parser.add_option("-s", "--slaves", type="int" ,dest="Nslaves",
                      default=1, metavar="#", 
                      help="Number of slaves to start [default=%default]")

    parser.add_option("-p", "--port", dest="momfbd_port",
                      default=momfbd_port_default,metavar="#",
                      help="MOMFBD port, defaults to a free port")

    parser.add_option("-T", "--Tsleep", type="int", dest="tsleep", 
                      default=3,metavar="#",
                      help="""Time to sleep between jstats in seconds [default=%default]""")

    parser.add_option("-i", "--iterations", type="int", dest="maxiterations", 
                      default=10, metavar="#",
                      help="Max number of iterations [default=%default]")

    parser.add_option("-S", "--shiftcrit", type="int", dest="shiftcriterion", 
                      default=10, metavar="#",
                      help="Shift criterion for stopping iterations. In hundedths of a pixel [default: #=%default]")

    parser.add_option("-F", "--focuscrit", type="float", dest="focuscriterion", 
                      default=0.1, metavar="#",
                      help="Focus criterion for stopping iterations. In mm. [default: #=%default]")

    parser.add_option("-a", "--align", dest="align", default=False,
                      action="store_true", 
                      help="""Align by cross correlation to get starting offsets (not implemented yet)""") 

    parser.add_option("-r", "--removespikes", dest="removespikes",
                      default=False, action="store_true", 
                      help="Remove spikes from estimated quantities by median filtering")

    parser.add_option("-v", "--verbose", action="store_true",
                      dest="verbose", default=False, 
                      help="Print status messages to stdout") 

    parser.add_option("-t", "--test", action="store_true",
                      dest="test", default=False, 
                      help="Run test() indstead of main()") 


    (options, cfg_files) = parser.parse_args()

    cfg_file = cfg_files[0]
    print cfg_file
    
    f = open(cfg_file,'r')
    cfglines = f.readlines()
    f.close()

    

    print getcfgtag(cfglines,'DIVERSITY',options,skip=['channel','channel']).split()[0]

    
  

if(__name__ == "__main__"):
#    print sys.argv
    if "-t" in sys.argv or "--test" in sys.argv:
        test()
    else:
        main()
