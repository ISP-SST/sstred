ARCH     = $(shell uname -m)
ARCHF    = $(shell uname -a)
CC = gcc-mp-4.6
CXX = g++-mp-4.6
MDIR = /opt/local/
CPPFLAGS =-O3 -fPIC -march=core2 -I$(MDIR)/include #-no-cpp-precomp -dynamic -m64
OFILES = mymath.o math_wrap.o
INSTFILES = math.so math.dlm
DLM_PATH = $(HOME)/idl/dlm/$(ARCH)

ifneq (,$(findstring Darwin,$(ARCHF))) #
  #CPPFLAGS += -arch x86_64
  LIBTOOL = /usr/bin/glibtool
  INSTALL = install
  MACOSX_DEPLOYMENT_TARGET=10.7
  export MACOSX_DEPLOYMENT_TARGET
endif
export bla

all: clean dlm install

clean:
	rm -rf *.o *.lo *.la .libs
        
dlm:
	$(LIBTOOL) --mode=compile --tag=CC $(CC) $(CFLAGS) -o math_wrap.lo -c math_wrap.c
	$(LIBTOOL) --mode=link    --tag=CC $(CC) -bundle -o libmath.la math_wrap.lo -rpath $(DLM_PATH)

install:
	$(LIBTOOL) --mode=install --tag=CC $(INSTALL) -c libmath.la $(DLM_PATH)
	$(LIBTOOL) --mode=finish  --tag=CC --finish $(DLM_PATH)
	rm -f $(DLM_PATH)/libmath.dlm
	cp libmath.dlm $(DLM_PATH)

