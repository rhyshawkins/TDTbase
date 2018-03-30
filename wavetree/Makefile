
INCLUDES = -I../log \
	-I../oset \
	-I../sphericalwavelet

CC ?= gcc
CFLAGS = -c -g -Wall -fPIC $(INCLUDES)

CFLAGS += -O2

AR = ar
ARFLAGS = -r

INSTALL = install
INSTALLFLAGS = -D

LIBS = 

TARGETS = libwavetree.a

OBJS = \
	chain_history.o \
	coefficient_histogram.o \
	subdivisiontree2d.o \
	wavetree2d.o \
	wavetree2d_sub.o \
	wavetree3d.o \
	wavetree3d_sub.o \
	wavetreepp.o \
	wavetreesphere3d.o \
	wavetreesphereface2d.o \
	wavetreesphereface3d.o \
	wavetree_prior.o \
	wavetree_prior_globally_uniform.o \
	wavetree_prior_globally_laplacian.o \
	wavetree_prior_depth_uniform.o \
	wavetree_prior_depth_generalised_gaussian.o \
	wavetree_birth_proposal.o \
	wavetree_value_proposal.o \
	wavetree_value_proposal_cauchy_am.o \
	wavetree_value_proposal_gaussian_am.o \
	wavetree_value_proposal_global_gaussian.o \
	wavetree_value_proposal_depth_gaussian.o \
	wavetree_value_proposal_depth_cauchy.o \
	wavetree_value_proposal_depth_gaussian_scam.o

SRCS = Makefile \
	chain_history.c \
	chain_history.h \
	coefficient_histogram.c \
	coefficient_histogram.h \
	subdivisiontree2d.c \
	subdivisiontree2d.h \
	wavetree.h \
	wavetree2d.c \
	wavetree2d.h \
	wavetree2d_sub.c \
	wavetree2d_sub.h \
	wavetree3d.c \
	wavetree3d.h \
	wavetree3d_sub.c \
	wavetree3d_sub.h \
	wavetree_birth_proposal.c \
	wavetree_birth_proposal.h \
	wavetree_prior.c \
	wavetree_prior.h \
	wavetree_prior_depth_generalised_gaussian.c \
	wavetree_prior_depth_uniform.c \
	wavetree_prior_globally_uniform.c \
	wavetree_prior_globally_laplacian.c \
	wavetree_value_proposal.c \
	wavetree_value_proposal.h \
	wavetree_value_proposal_cauchy_am.c \
	wavetree_value_proposal_depth_cauchy.c \
	wavetree_value_proposal_depth_gaussian.c \
	wavetree_value_proposal_depth_gaussian_am.c \
	wavetree_value_proposal_depth_gaussian_scam.c \
	wavetree_value_proposal_gaussian_am.c \
	wavetree_value_proposal_global_gaussian.c \
	wavetreepp.c \
	wavetreepp.h \
	wavetreesphere2d.h \
	wavetreesphere3d.c \
	wavetreesphere3d.h \
	wavetreesphereface2d.c \
	wavetreesphereface2d.h \
	wavetreesphereface3d.c \
	wavetreesphereface3d.h \
	Makefile \
	tests/lanczos_images.c \
	tests/pyramid_images.c \
	tests/subdivisiontree2d_basis_tests.c \
	tests/subdivisiontree2d_lanczos_tests.c \
	tests/subdivisiontree2d_pyramid_tests.c \
	tests/subdivisiontree2d_tests.c \
	tests/wavetree2d_tests.c \
	tests/wavetree3d_tests.c \
	tests/wavetree_prior_tests.c \
	tests/wavetree_value_proposal_tests.c \
	tests/wavetreesphere3d_tests.c

all : $(TARGETS)

libwavetree.a : $(OBJS)
	$(AR) $(ARFLAGS) libwavetree.a $(OBJS) 

%.o : %.c
	$(CC) $(CFLAGS) -o $*.o $*.c

DATE = $(shell date +"%Y%m%d%H%M")
DIR = wavetree
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
