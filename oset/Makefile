
INCLUDES = $(shell gsl-config --cflags) \
	-I../log

CC ?= gcc
CFLAGS = -c -g -Wall -fPIC $(INCLUDES)

AR = ar
ARFLAGS = -r

LIBS = -L../log -llog

TARGETS = liboset.a

OBJS = oset_gmpz.o \
	oset_int.o \
	oset_int_double.o \
	multiset_int.o \
	multiset_int_double.o \
	ohist_64.o \
	ohist_64_set.o \
	ttree_int.o

INSTALL = install
INSTALLFLAGS = -D

SRCS = Makefile \
	multiset_int.c \
	multiset_int.h \
	multiset_int_double.c \
	multiset_int_double.h \
	ohist_64.c \
	ohist_64.h \
	ohist_64_set.c \
	ohist_64_set.h \
	oset_gmpz.c \
	oset_gmpz.h \
	oset_int.c \
	oset_int.h \
	oset_int_double.c \
	oset_int_double.h \
	ttree_int.c \
	ttree_int.h \
	tests/Makefile \
	tests/multiset_int_double_tests.c \
	tests/multiset_int_tests.c \
	tests/ohist_64_tests.c \
	tests/oset_int_tests.c \
	tests/ttree_int_tests.c

all : $(TARGETS)

liboset.a : $(OBJS)
	$(AR) $(ARFLAGS) liboset.a $(OBJS) 

%.o : %.c
	$(CC) $(CFLAGS) -o $*.o $*.c

DATE = $(shell date +"%Y%m%d%H%M")
DIR = oset
TGZ = $(DIR).tar.gz

dist :
	mkdir -p $(DIR)
	echo $(DATA) > $(DIR)/Version
	for f in $(SRCS) $(EXTRADIST); do \
	    $(INSTALL) $(INSTALLFLAGS) $$f $(DIR)/$$f ; \
	done
	tar -czf $(TGZ) $(DIR)/*
	rm -rf $(DIR)

clean : 
	rm -f $(TARGETS) *.o
