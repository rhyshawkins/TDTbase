
INCLUDES = -I../ \
	-I../../wavelet

CC = gcc
CFLAGS = -c -g -Wall $(INCLUDES) \
	$(shell pkg-config --cflags check) \
	$(shell gsl-config --cflags)

LIBS = -L../ -lsphericalwavelet \
	-L../../wavelet/ -l wavelet \
	-L../../log -llog \
	-lm \
	$(shell pkg-config --libs check) \
	$(shell gsl-config --libs)

TARGETS = icosahedron_tests \
	octahedron_tests \
	tetrahedron_tests \
	butterfly_tests \
	butterfly_lifted_tests \
	triangle_tests \
	face_subdivision_tests \
	face_wavelet_tests \
	icosahedron_resolution

all : $(TARGETS)

icosahedron_tests: icosahedron_tests.o
	$(CC) -o icosahedron_tests icosahedron_tests.o $(LIBS)

tetrahedron_tests: tetrahedron_tests.o
	$(CC) -o tetrahedron_tests tetrahedron_tests.o $(LIBS)

octahedron_tests: octahedron_tests.o
	$(CC) -o octahedron_tests octahedron_tests.o $(LIBS)

butterfly_tests: butterfly_tests.o
	$(CC) -o butterfly_tests butterfly_tests.o $(LIBS)

butterfly_lifted_tests: butterfly_lifted_tests.o
	$(CC) -o butterfly_lifted_tests butterfly_lifted_tests.o $(LIBS)

triangle_tests: triangle_tests.o
	$(CC) -o triangle_tests triangle_tests.o $(LIBS)

face_subdivision_tests: face_subdivision_tests.o
	$(CC) -o face_subdivision_tests face_subdivision_tests.o $(LIBS)

face_wavelet_tests: face_wavelet_tests.o
	$(CC) -o face_wavelet_tests face_wavelet_tests.o $(LIBS)

icosahedron_resolution: icosahedron_resolution.o
	$(CC) -o icosahedron_resolution icosahedron_resolution.o $(LIBS)

%.o : %.c
	$(CC) $(CFLAGS) -o $*.o $*.c

clean :
	rm -f $(TARGETS) *.o
