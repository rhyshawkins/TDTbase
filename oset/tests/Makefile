
INCLUDES = -I../

CC = gcc
CFLAGS = -c -g -Wall $(INCLUDES) $(shell pkg-config --cflags check)

LIBS = -L../ -loset -lm -lgmp $(shell pkg-config --libs check)

TARGETS = oset_int_tests \
	ohist_64_tests \
	multiset_int_tests \
	multiset_int_double_tests \
	ttree_int_tests

all : $(TARGETS)

oset_int_tests: oset_int_tests.o
	$(CC) -o oset_int_tests oset_int_tests.o $(LIBS)

ohist_64_tests: ohist_64_tests.o
	$(CC) -o ohist_64_tests ohist_64_tests.o $(LIBS)

multiset_int_tests: multiset_int_tests.o
	$(CC) -o multiset_int_tests multiset_int_tests.o $(LIBS)

multiset_int_double_tests: multiset_int_double_tests.o
	$(CC) -o multiset_int_double_tests multiset_int_double_tests.o $(LIBS)

ttree_int_tests: ttree_int_tests.o
	$(CC) -o ttree_int_tests ttree_int_tests.o $(LIBS)

%.o : %.c
	$(CC) $(CFLAGS) -o $*.o $*.c

clean :
	rm -f $(TARGETS) *.o
