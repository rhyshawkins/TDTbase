
PHDBASE=..

INCLUDES = -I$(PHDBASE)/log \
	-I$(PHDBASE)/wavelet \
	$(shell gsl-config --cflags)

CC ?= gcc
CFLAGS = -c -g -Wall $(INCLUDES)

OBJS = sphericalwavelet.o \
	edge.o \
	face_subdivision.o \
	face_wavelet.o \
	manifold.o \
	manifold_validate.o \
	icosahedron.o \
	octahedron.o \
	radialbasisfunction.o \
	spherical.o \
	tetrahedron.o \
	transform.o \
	triangle.o \
	vertex3.o \
	vertex_wavelet.o

AR = ar
ARFLAGS = -r

INSTALL = install
INSTALLFLAGS = -D

LIBS = 

TARGETS = libsphericalwavelet.a

SRCS = Makefile \
	edge.c \
	edge.h \
	face_subdivision.c \
	face_subdivision.h \
	face_wavelet.c \
	face_wavelet.h \
	icosahedron.c \
	icosahedron.h \
	manifold.c \
	manifold.h \
	manifold_validate.c \
	octahedron.c \
	octahedron.h \
	radialbasisfunction.c \
	radialbasisfunction.h \
	spherical.c \
	spherical.h \
	sphericalwavelet.c \
	sphericalwavelet.h \
	tetrahedron.c \
	tetrahedron.h \
	transform.c \
	transform.h \
	triangle.c \
	triangle.h \
	vertex3.c \
	vertex3.h \
	vertex_wavelet.c \
	vertex_wavelet.h \
	tests/Makefile \
	tests/butterfly_lifted_tests.c \
	tests/butterfly_tests.c \
	tests/icosahedron_tests.c \
	tests/octahedron_tests.c \
	tests/tetrahedron_tests.c \
	tests/triangle_tests.c 

all : $(TARGETS)

libsphericalwavelet.a : $(OBJS)
	$(AR) $(ARFLAGS) libsphericalwavelet.a $(OBJS) 

%.o : %.c
	$(CC) $(CFLAGS) -o $*.o $*.c

DATE = $(shell date +"%Y%m%d%H%M")
DIR = sphericalwavelet
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
