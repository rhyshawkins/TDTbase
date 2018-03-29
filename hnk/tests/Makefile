
INCLUDES = -I../

CC = gcc
CFLAGS = -c -g -Wall $(INCLUDES) $(shell pkg-config --cflags check)

LIBS = -L../ -lhnk -lm -lgmp $(shell pkg-config --libs check) -L../../log -llog

TARGETS = hnk_aggregate_tests \
	hnk_btree_tests \
	hnk_qtree_tests \
	hnk_octree_tests \
	hnk_unary_tests \
	hnk_cart12_tests \
	hnk_cart34_tests \
	hnk_cartesian_nonsquare_tests \
	hnk_cart78_tests \
	hnk_butterfly_tests \
	hnk_facewaveletsphere_tests \
	hnk_triangle_tests \
	hnk_s2_tests

all : $(TARGETS)

hnk_aggregate_tests: hnk_aggregate_tests.o
	$(CC) -o hnk_aggregate_tests hnk_aggregate_tests.o $(LIBS)

hnk_btree_tests: hnk_btree_tests.o
	$(CC) -o hnk_btree_tests hnk_btree_tests.o $(LIBS)

hnk_qtree_tests: hnk_qtree_tests.o
	$(CC) -o hnk_qtree_tests hnk_qtree_tests.o $(LIBS)

hnk_octree_tests: hnk_octree_tests.o
	$(CC) -o hnk_octree_tests hnk_octree_tests.o $(LIBS)

hnk_unary_tests: hnk_unary_tests.o
	$(CC) -o hnk_unary_tests hnk_unary_tests.o $(LIBS)

hnk_cart12_tests: hnk_cart12_tests.o
	$(CC) -o hnk_cart12_tests hnk_cart12_tests.o $(LIBS)

hnk_cart34_tests: hnk_cart34_tests.o
	$(CC) -o hnk_cart34_tests hnk_cart34_tests.o $(LIBS)

hnk_s2_tests: hnk_s2_tests.o
	$(CC) -o hnk_s2_tests hnk_s2_tests.o $(LIBS)

hnk_cartesian_nonsquare_tests: hnk_cartesian_nonsquare_tests.o
	$(CC) -o hnk_cartesian_nonsquare_tests hnk_cartesian_nonsquare_tests.o $(LIBS)

hnk_cart78_tests: hnk_cart78_tests.o
	$(CC) -o hnk_cart78_tests hnk_cart78_tests.o $(LIBS)

hnk_butterfly_tests: hnk_butterfly_tests.o
	$(CC) -o hnk_butterfly_tests hnk_butterfly_tests.o $(LIBS)

hnk_facewaveletsphere_tests: hnk_facewaveletsphere_tests.o
	$(CC) -o hnk_facewaveletsphere_tests hnk_facewaveletsphere_tests.o $(LIBS)

hnk_triangle_tests: hnk_triangle_tests.o
	$(CC) -o hnk_triangle_tests hnk_triangle_tests.o $(LIBS)

%.o : %.c
	$(CC) $(CFLAGS) -o $*.o $*.c

clean :
	rm -f $(TARGETS) *.o
