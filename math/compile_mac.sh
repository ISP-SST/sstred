#!/bin/bash
CC=gcc-mp-4.7
CXX=g++-mp-4.7
OPTS='-O3 -fPIC -fpermissive -flto -ffast-math'
LDOPTS=''
rm *.o *.so
$CXX $OPTS -c mymath.cc
$CXX $OPTS -c math_wrap.c
$CXX -bundle $LDOPTS -o libmath.so mymath.o math_wrap.o -lm -lc -flat_namespace -undefined suppress #-macosx_version_min 10.7
echo copying result
cp libmath.so libmath.dlm ~/idl/dlm/x86_64/
