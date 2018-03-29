INCLUDES = -I../log

CC ?= gcc
CFLAGS = -c -g -Wall  -fPIC $(INCLUDES)

CXX ?= g++
CXXFLAGS = -c -g -Wall  -fPIC $(INCLUDES)

AR = ar
ARFLAGS = -r

LIBS = -lm $(shell gsl-config --libs) -lgmp \
	-L../log -llog

OBJS = hnk.o \
	hnk_butterfly.o \
	hnk_cartesian.o \
	hnk_cartesian_nonsquare.o \
	hnk_standard.o \
	hnk_facewaveletsphere.o \
	hnk_healpix.o \
	hnk_triangle.o \
	hnk_s2.o

TARGETS = libhnk.a \
	hnk_cartesian_timing

INSTALL = install
INSTALLFLAGS = -D

SRCS = Makefile \
	hnk.h \
	hnk_butterfly.h \
	hnk_cartesian.h \
	hnk_cartesian_nonsquare.h \
	hnk_standard.h \
	hnk_healpix.h \
	hnk_s2.h \
	hnk.c \
	hnk_butterfly.c \
	hnk_cartesian.c \
	hnk_cartesian_nonsquare.c \
	hnk_cartesian_timing.c \
	hnk_standard.c \
	hnk_facewaveletsphere.c \
	hnk_facewaveletsphere.h \
	hnk_healpix.c \
	hnk_s2.c \
	hnk_triangle.h \
	hnk_triangle.c \
	tests/Makefile \
	tests/hnk_aggregate_tests.c \
	tests/hnk_btree_tests.c \
	tests/hnk_butterfly_tests.c \
	tests/hnk_cart12_tests.c \
	tests/hnk_cart34_tests.c \
	tests/hnk_cart78_tests.c \
	tests/hnk_cartesian_nonsquare_tests.c \
	tests/hnk_facewaveletsphere_tests.c \
	tests/hnk_octree_tests.c \
	tests/hnk_qtree_tests.c \
	tests/hnk_unary_tests.c

EXTRADIST =

all : $(TARGETS)

libhnk.a : $(OBJS)
	$(AR) $(ARFLAGS) libhnk.a $(OBJS)

hnk_cartesian_timing : hnk_cartesian_timing.o $(OBJS)
	$(CC) -o hnk_cartesian_timing hnk_cartesian_timing.o $(OBJS) $(LIBS)

%.o : %.c
	$(CC) $(CFLAGS) -o $*.o $*.c

DATE = $(shell date +"%Y%m%d%H%M")
DIR = hnk
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

