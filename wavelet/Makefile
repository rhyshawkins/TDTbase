
PHDLIBBASE=..
INCLUDES = -I$(PHDLIBBASE)/log \
	$(shell gsl-config --cflags)

CC ?= gcc
CFLAGS = -c -g -Wall $(INCLUDES) -fPIC

CFLAGS += -O2

AR = ar
ARFLAGS = -r

INSTALL = install
INSTALLFLAGS = -D

LIBS = 

TARGETS = libwavelet.a

OBJS = boundary.o \
	cdf97_lift.o \
	cdf97_lift_periodic.o \
	cdf97_lift_impulse.o \
	cdf97_matrix.o \
	daub4_lift.o \
	daub4_matrix.o \
	daub4_dwt.o \
	daub6_dwt.o \
	daub8_dwt.o \
	haar_lift.o \
	haar_matrix.o \
	generic_lift.o \
	generic_matrix.o

SRCS = Makefile \
	boundary.c \
	boundary.h \
	cdf97_healpix.c \
	cdf97_healpix.h \
	cdf97_lift.c \
	cdf97_lift.h \
	cdf97_lift_periodic.c \
	cdf97_lift_periodic.h \
	cdf97_lift_impulse.c \
	cdf97_lift_impulse.h \
	cdf97_matrix.c \
	cdf97_matrix.h \
	daub4_dwt.c \
	daub4_dwt.h \
	daub6_dwt.c \
	daub6_dwt.h \
	daub8_dwt.c \
	daub8_dwt.h \
	daub4_lift.c \
	daub4_lift.h \
	daub4_matrix.c \
	daub4_matrix.h \
	daubechies.c \
	daubechies.h \
	generic_lift.c \
	generic_lift.h \
	generic_matrix.c \
	generic_matrix.h \
	haar_lift.c \
	haar_lift.h \
	haar_matrix.c \
	haar_matrix.h \
	wavelet2d_tile.c \
	wavelet2d_tile.h \
	doc/matrix/matrix_templates.tex	\
	Makefile \
	tests/cdf97_healpix_tests.c \
	tests/cdf97_impulse_tests.c \
	tests/cdf97_lift_tests.c \
	tests/cdf97_matrix_tests.c \
	tests/daub4_dwt_tests.c \
	tests/daub4_lift_tests.c \
	tests/daub4_matrix_tests.c \
	tests/daub4_tests.c \
	tests/generic_matrix_tests.c \
	tests/haar_lift_tests.c \
	tests/haar_matrix_tests.c 

all : $(TARGETS)

libwavelet.a : $(OBJS)
	$(AR) $(ARFLAGS) libwavelet.a $(OBJS) 

%.o : %.c
	$(CC) $(CFLAGS) -o $*.o $*.c

DATE = $(shell date +"%Y%m%d%H%M")
DIR = wavelet
TGZ = $(DIR).tar.gz

dist :
	mkdir -p $(DIR)
	echo $(DATE) > $(DIR)/Version
	for f in $(SRCS) $(EXTRADIST); do \
	    $(INSTALL) $(INSTALLFLAGS) $$f $(DIR)/$$f ; \
	done
	tar -czf $(TGZ) $(DIR)/*
	rm -rf $(DIR)

clean : 
	rm -f $(TARGETS) *.o
