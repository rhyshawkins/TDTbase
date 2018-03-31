
CC ?= gcc
CFLAGS = -c -g -Wall

AR = ar
ARFLAGS = -r

INSTALL = install
INSTALLFLAGS = -D

LIBS = 

TARGETS = libtracking.a

OBJS = tracking.o

SRCS = Makefile \
	tracking.c \
	tracking.h 
all : $(TARGETS)

libtracking.a : $(OBJS)
	$(AR) $(ARFLAGS) libtracking.a $(OBJS) 

%.o : %.c
	$(CC) $(CFLAGS) -o $*.o $*.c

DATE = $(shell date +"%Y%m%d%H%M")
DIR = tracking
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
