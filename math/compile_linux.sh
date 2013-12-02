#!/bin/bash
CC=gcc
CXX=g++
OPTS='-O3 -fPIC -fpermissive -ffast-math'
LDOPTS='-flto -march=native'
rm *.o *.so
$CXX $OPTS -c mymath.cc
$CXX $OPTS -c math_wrap.c
$CXX -shared $LDOPTS -o libmath.so mymath.o math_wrap.o -lm -lc  #-macosx_version_min 10.7
echo copying result
cp libmath.so libmath.dlm ~/idl/dlm/x86_64/
