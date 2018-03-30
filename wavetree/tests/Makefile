
INCLUDES = -I../

CC = gcc
CFLAGS = -c -g -Wall $(INCLUDES) \
	-I../../sphericalwavelet \
	-I../../oset \
	$(shell gsl-config --libs) \
	$(shell pkg-config --cflags check)

LIBS = -L../ -lwavetree \
	-L../../oset -loset \
	-L../../sphericalwavelet -lsphericalwavelet \
	-L../../log -llog \
	-lm -lgmp \
	$(shell gsl-config --libs) \
	$(shell pkg-config --libs check)

TARGETS = wavetree2d_tests \
	wavetree2d_sub_tests \
	wavetree3d_tests \
	wavetree3d_sub_tests \
	wavetreesphere3d_tests \
	wavetreesphereface2d_tests \
	wavetreesphereface3d_tests \
	subdivisiontree2d_tests \
	subdivisiontree2d_basis_tests \
	subdivisiontree2d_pyramid_tests \
	subdivisiontree2d_lanczos_tests \
	wavetree_prior_tests \
	wavetree_value_proposal_tests \
	pyramid_images \
	lanczos_images

all : $(TARGETS)

wavetree2d_tests: wavetree2d_tests.o
	$(CC) -o wavetree2d_tests wavetree2d_tests.o $(LIBS)

wavetree2d_sub_tests: wavetree2d_sub_tests.o
	$(CC) -o wavetree2d_sub_tests wavetree2d_sub_tests.o $(LIBS)

wavetree3d_tests: wavetree3d_tests.o
	$(CC) -o wavetree3d_tests wavetree3d_tests.o $(LIBS)

wavetree3d_sub_tests: wavetree3d_sub_tests.o
	$(CC) -o wavetree3d_sub_tests wavetree3d_sub_tests.o $(LIBS)

wavetreesphere3d_tests: wavetreesphere3d_tests.o
	$(CC) -o wavetreesphere3d_tests wavetreesphere3d_tests.o $(LIBS)

wavetreesphereface2d_tests: wavetreesphereface2d_tests.o
	$(CC) -o wavetreesphereface2d_tests wavetreesphereface2d_tests.o $(LIBS)

wavetreesphereface3d_tests: wavetreesphereface3d_tests.o
	$(CC) -o wavetreesphereface3d_tests wavetreesphereface3d_tests.o $(LIBS)

subdivisiontree2d_tests: subdivisiontree2d_tests.o
	$(CC) -o subdivisiontree2d_tests subdivisiontree2d_tests.o $(LIBS)

subdivisiontree2d_basis_tests: subdivisiontree2d_basis_tests.o
	$(CC) -o subdivisiontree2d_basis_tests subdivisiontree2d_basis_tests.o $(LIBS)

subdivisiontree2d_pyramid_tests: subdivisiontree2d_pyramid_tests.o
	$(CC) -o subdivisiontree2d_pyramid_tests subdivisiontree2d_pyramid_tests.o $(LIBS)

subdivisiontree2d_lanczos_tests: subdivisiontree2d_lanczos_tests.o
	$(CC) -o subdivisiontree2d_lanczos_tests subdivisiontree2d_lanczos_tests.o $(LIBS)

wavetree_prior_tests : wavetree_prior_tests.o
	$(CC) -o wavetree_prior_tests wavetree_prior_tests.o $(LIBS)

wavetree_value_proposal_tests : wavetree_value_proposal_tests.o
	$(CC) -o wavetree_value_proposal_tests wavetree_value_proposal_tests.o $(LIBS)

pyramid_images : pyramid_images.o
	$(CC) -o pyramid_images pyramid_images.o $(LIBS)

lanczos_images : lanczos_images.o
	$(CC) -o lanczos_images lanczos_images.o $(LIBS)

%.o : %.c
	$(CC) $(CFLAGS) -o $*.o $*.c

clean :
	rm -f $(TARGETS) *.o
