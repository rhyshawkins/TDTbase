
INCLUDES = -I../

CC = gcc
CFLAGS = -c -g -Wall $(INCLUDES) \
	$(shell gsl-config --cflags) \
	$(shell pkg-config --cflags check)

LIBS = -L../ -lwavelet \
	-L../../log -llog \
	-lm \
	$(shell gsl-config --libs) \
	$(shell pkg-config --libs check)

TARGETS = cdf97_lift_tests \
	cdf97_lift_periodic_tests \
	cdf97_matrix_tests \
	cdf97_impulse_tests \
	daub4_dwt_tests \
	daub6_dwt_tests \
	daub8_dwt_tests \
	daub4_lift_tests \
	daub4_matrix_tests \
	haar_lift_tests \
	haar_matrix_tests \
	generic_matrix_tests \
	generic_lift_tests 

all : $(TARGETS)

cdf97_lift_tests: cdf97_lift_tests.o
	$(CC) -o cdf97_lift_tests cdf97_lift_tests.o $(LIBS)

cdf97_lift_periodic_tests: cdf97_lift_periodic_tests.o
	$(CC) -o cdf97_lift_periodic_tests cdf97_lift_periodic_tests.o $(LIBS)

cdf97_matrix_tests: cdf97_matrix_tests.o
	$(CC) -o cdf97_matrix_tests cdf97_matrix_tests.o $(LIBS)

cdf97_impulse_tests: cdf97_impulse_tests.o
	$(CC) -o cdf97_impulse_tests cdf97_impulse_tests.o $(LIBS)

daub4_dwt_tests: daub4_dwt_tests.o
	$(CC) -o daub4_dwt_tests daub4_dwt_tests.o $(LIBS)

daub6_dwt_tests: daub6_dwt_tests.o
	$(CC) -o daub6_dwt_tests daub6_dwt_tests.o $(LIBS)

daub8_dwt_tests: daub8_dwt_tests.o
	$(CC) -o daub8_dwt_tests daub8_dwt_tests.o $(LIBS)

haar_lift_tests: haar_lift_tests.o
	$(CC) -o haar_lift_tests haar_lift_tests.o $(LIBS)

haar_matrix_tests: haar_matrix_tests.o
	$(CC) -o haar_matrix_tests haar_matrix_tests.o $(LIBS)

daub4_lift_tests: daub4_lift_tests.o
	$(CC) -o daub4_lift_tests daub4_lift_tests.o $(LIBS)

daub4_matrix_tests: daub4_matrix_tests.o
	$(CC) -o daub4_matrix_tests daub4_matrix_tests.o $(LIBS)

generic_matrix_tests: generic_matrix_tests.o
	$(CC) -o generic_matrix_tests generic_matrix_tests.o $(LIBS)

generic_lift_tests: generic_lift_tests.o
	$(CC) -o generic_lift_tests generic_lift_tests.o $(LIBS)

%.o : %.c
	$(CC) $(CFLAGS) -o $*.o $*.c

clean :
	rm -f $(TARGETS) *.o
