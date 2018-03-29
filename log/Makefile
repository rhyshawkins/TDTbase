
CC ?= gcc
CFLAGS = -c -g -Wall -fPIC

AR = ar
ARFLAGS = -r

INSTALL = install
INSTALLFLAGS = -D

LIBS = 

TARGETS = liblog.a

OBJS = slog.o

SRCS = Makefile \
	slog.c \
	slog.h 
all : $(TARGETS)

liblog.a : $(OBJS)
	$(AR) $(ARFLAGS) liblog.a $(OBJS) 

%.o : %.c
	$(CC) $(CFLAGS) -o $*.o $*.c

DATE = $(shell date +"%Y%m%d%H%M")
DIR = log
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
