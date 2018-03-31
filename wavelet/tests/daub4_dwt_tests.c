
#include <stdio.h>
#include <stdlib.h>
#include <check.h>
#include <math.h>

#include "daub4_dwt.h"
#include "boundary.h"

#define WIDTH 32

#define NS1_WIDTH 32
#define NS1_HEIGHT 16

#define NS2_WIDTH 16
#define NS2_HEIGHT 64

#define NS3_WIDTH 16
#define NS3_HEIGHT 16
#define NS3_DEPTH 32

#define NS4_WIDTH 16
#define NS4_HEIGHT 64
#define NS4_DEPTH 32

static inline double within1pc(double a, double b)
{
  double delta;

  delta = fabs(a - b);
  if (delta < 1.0e-6) {
    return -1;
  }

  if (fabs(a) > fabs(b)) {
    if ((delta/fabs(a)) < 0.01) {
      return -1;
    } else {
      fprintf(stderr, "within1p: %f %f\n", a, b);
      return 0;
    }
  } else {
    if ((delta/fabs(b)) < 0.01) {
      return -1;
    } else {
      fprintf(stderr, "within1p: %f %f\n", a, b);
      return 0;
    }
  }
}

START_TEST (test_daub4_1d_simple_step)
{
#define OVERLAP 8
  double data[WIDTH + 8];
  double expected[WIDTH];
  double work[WIDTH];
  int i;


#define H0 0.34150635094610965
#define H1 0.5915063509461096
#define H2 0.15849364905389032
#define H3 -9.150635094610965e-2

#define G0 H3
#define G1 -H2
#define G2 H1
#define G3 -H0
  
  for (i = 0; i < WIDTH; i ++) {
    data[i] = 0.1 * (double)i;
  }
  for (i = WIDTH; i < (WIDTH + OVERLAP); i ++) {
    data[i] = data[wavelet_boundary_periodic(i, WIDTH)];
  }

  for (i = 0; i < WIDTH/2; i ++) {
    expected[i] =
      H0*data[2*i] + 
      H1*data[2*i + 1] + 
      H2*data[2*i + 2] + 
      H3*data[2*i + 3];

    expected[WIDTH/2 + i] =
      G0*data[2*i] + 
      G1*data[2*i + 1] + 
      G2*data[2*i + 2] + 
      G3*data[2*i + 3];
  }

  ck_assert(daub4_dwt_forward1d_daub4_step(data, 
					   WIDTH,
					   1, 
					   work) >= 0);

  /* for (i = 0; i < WIDTH; i ++) { */
  /*   printf("%2d: %15.9f %15.9f %15.9f\n", i, data[i], expected[i], fabs(data[i] - expected[i])); */
  /* } */

  for (i = 0; i < WIDTH; i ++) {
    ck_assert(fabs(data[i] - expected[i]) < 1.0e-6);
  }
}
END_TEST

START_TEST (test_daub4_1d_constant_step)
{
  double data[WIDTH];
  double work[WIDTH];
  int i;

  for (i = 0; i < (WIDTH); i ++) {
    data[i] = 3.1415;
  }

  ck_assert(daub4_dwt_forward1d_daub4_step(data, 
					    WIDTH,
					    1, 
					    work) >= 0);

  for (i = 1; i < (WIDTH); i ++) {
    ck_assert(isnormal(data[i]) || data[i] == 0.0);
  }

  for (i = 0; i < (WIDTH/2); i ++) {
    ck_assert(within1pc(data[i], 3.1415));
  }

  for (i = WIDTH/2; i < WIDTH; i ++) {
    ck_assert(fabs(data[i]) < 1.0e-6);
  }
  
}
END_TEST

START_TEST (test_daub4_1d_sinusoid_step)
{
  double data[WIDTH];
  double recon[WIDTH];
  double work[WIDTH];
  int i;

  for (i = 0; i < (WIDTH); i ++) {
    data[i] = sin((double)i/64.0 * M_PI * 2.0);
    recon[i] = data[i];
  }

  ck_assert(daub4_dwt_forward1d_daub4_step(recon, 
					   WIDTH,
					   1, 
					   work) >= 0);

  for (i = 1; i < (WIDTH); i ++) {
    ck_assert(isnormal(recon[i]) || recon[i] == 0.0);
  }

  /* printf("1D Forward Step Result: sinusoid\n"); */
  /* for (i = 0; i < (WIDTH); i ++) { */
  /*   printf("%3d %10.6f %10.6f\n", i, recon[i], data[i]); */
  /* } */

  ck_assert(daub4_dwt_inverse1d_daub4_step(recon, 
					   WIDTH,
					   1, 
					   work) >= 0);

  /* printf("1D Inverse Step Result: sinusoid\n"); */
  /* for (i = 0; i < (WIDTH); i ++) { */
  /*   printf("%3d %10.6f %10.6f\n", i, recon[i], data[i]); */
  /* } */

  for (i = 0; i < WIDTH; i ++) {
    ck_assert(fabs(recon[i] - data[i]) < 1.0e-6);
  }
  
}
END_TEST

START_TEST (test_daub4_1d_constant)
{
  double data[WIDTH];
  double work[WIDTH];
  int i;

  for (i = 0; i < (WIDTH); i ++) {
    data[i] = 3.1415;
  }

  ck_assert(daub4_dwt_forward1d_daub4(data, 
				       WIDTH,
				       1, 
				       work) >= 0);

  for (i = 1; i < (WIDTH); i ++) {
    ck_assert(isnormal(data[i]) || data[i] == 0.0);
  }

  ck_assert(within1pc(data[0], 3.1415));

  for (i = 1; i < (WIDTH); i ++) {
    ck_assert(within1pc(data[i], 0.0));
  }
  
  ck_assert(daub4_dwt_inverse1d_daub4(data, 
				       WIDTH,
				       1, 
				       work) >= 0);
  
  for (i = 0; i < (WIDTH); i ++) {
    ck_assert(within1pc(data[i], 3.14));
  }
}
END_TEST

START_TEST (test_daub4_1d_sinusoid)
{
  double data[WIDTH];
  double recon[WIDTH];
  double work[WIDTH];
  int i;

  for (i = 0; i < WIDTH; i ++) {
    data[i] = sin((double)i/64.0 * M_PI * 2.0);
    recon[i] = data[i];
  }

  for (i = 0; i < (WIDTH); i ++) {
    ck_assert(isnormal(recon[i]) || isfinite(recon[i]));
  }

  ck_assert(daub4_dwt_forward1d_daub4(recon, 
				      WIDTH,
				      1, 
				      work) >= 0);

  for (i = 0; i < (WIDTH); i ++) {
    ck_assert(isnormal(recon[i]) || isfinite(recon[i]));
  }
  
  ck_assert(daub4_dwt_inverse1d_daub4(recon, 
				      WIDTH,
				      1,
				      work) >= 0);

  /* printf("1D Inverse Result: sinusoid\n"); */
  /* for (i = 0; i < (WIDTH); i ++) { */
  /*   printf("%3d %10.6f %10.6f\n", i, recon[i], data[i]); */
  /* } */

  for (i = 0; i < (WIDTH); i ++) {
    ck_assert(isnormal(recon[i]) || isfinite(recon[i]));
  }

  for (i = 0; i < (WIDTH); i ++) {
    ck_assert(within1pc(data[i], recon[i]));
  }
}
END_TEST

START_TEST (test_daub4_2d_constant)
{
  double data[WIDTH * WIDTH];
  double work[WIDTH];
  int i;

  for (i = 0; i < (WIDTH * WIDTH); i ++) {
    data[i] = 3.14;
  }

  ck_assert(daub4_dwt_forward2d_daub4(data, 
				      WIDTH,
				      WIDTH,
				      WIDTH, 
				      work,
				      0) >= 0);

  for (i = 1; i < (WIDTH * WIDTH); i ++) {
    ck_assert(isnormal(data[i]) || data[i] == 0.0);
  }

  ck_assert(within1pc(data[0], 3.14));

  for (i = 1; i < (WIDTH * WIDTH); i ++) {
    ck_assert(within1pc(data[i], 0.0));
  }
  
  ck_assert(daub4_dwt_inverse2d_daub4(data, 
				      WIDTH,
				      WIDTH,
				      WIDTH, 
				      work,
				      0) >= 0);
  
  for (i = 0; i < (WIDTH * WIDTH); i ++) {
    ck_assert(within1pc(data[i], 3.14));
  }
}
END_TEST

START_TEST (test_daub4_2d_sinusoid)
{
  double data[WIDTH * WIDTH];
  double recon[WIDTH * WIDTH];
  double work[WIDTH];
  int i;
  int j;

  for (j = 0; j < WIDTH; j ++) {
    
    for (i = 0; i < WIDTH; i ++) {
      data[j * WIDTH + i] = sin((double)i/64.0 * M_PI * 2.0) * cos((double)j/8.0 * M_PI * 2.0);
      recon[j * WIDTH + i] = data[j * WIDTH + i];
    }

  }

  for (i = 0; i < (WIDTH * WIDTH); i ++) {
    ck_assert(isnormal(recon[i]) || isfinite(recon[i]));
  }

  ck_assert(daub4_dwt_forward2d_daub4(recon, 
				      WIDTH,
				      WIDTH,
				      WIDTH, 
				      work,
				      0) >= 0);

  for (i = 0; i < (WIDTH * WIDTH); i ++) {
    ck_assert(isnormal(recon[i]) || isfinite(recon[i]));
  }
  
  ck_assert(daub4_dwt_inverse2d_daub4(recon, 
				      WIDTH,
				      WIDTH,
				      WIDTH, 
				      work,
				      0) >= 0);

  for (i = 0; i < (WIDTH * WIDTH); i ++) {
    ck_assert(isnormal(recon[i]) || isfinite(recon[i]));
  }

  for (i = 0; i < (WIDTH * WIDTH); i ++) {
    ck_assert(within1pc(data[i], recon[i]));
  }
}
END_TEST

START_TEST (test_daub4_2d_constant_nonsquare1)
{
  double data[NS1_WIDTH * NS1_HEIGHT];
  double work[NS1_WIDTH];
  int i;
  int size;

  size = (NS1_WIDTH * NS1_HEIGHT);

  for (i = 0; i < size; i ++) {
    data[i] = 3.14;
  }

  ck_assert(daub4_dwt_forward2d_daub4(data, 
				      NS1_WIDTH,
				      NS1_HEIGHT,
				      NS1_WIDTH, 
				      work,
				      0) >= 0);

  for (i = 1; i < size; i ++) {
    ck_assert(isnormal(data[i]) || data[i] == 0.0);
  }

  ck_assert(within1pc(data[0], 3.14));

  for (i = 1; i < size; i ++) {
    ck_assert(within1pc(data[i], 0.0));
  }
  
  ck_assert(daub4_dwt_inverse2d_daub4(data, 
				      NS1_WIDTH,
				      NS1_HEIGHT,
				      NS1_WIDTH, 
				      work,
				      0) >= 0);
  
  for (i = 0; i < size; i ++) {
    ck_assert(within1pc(data[i], 3.14));
  }
}
END_TEST

START_TEST (test_daub4_2d_constant_nonsquare2)
{
  double data[NS2_WIDTH * NS2_HEIGHT];
  double work[NS2_HEIGHT];
  int i;
  int size;

  size = (NS2_WIDTH * NS2_HEIGHT);

  for (i = 0; i < size; i ++) {
    data[i] = 3.14;
  }

  ck_assert(daub4_dwt_forward2d_daub4(data, 
				      NS2_WIDTH,
				      NS2_HEIGHT,
				      NS2_WIDTH, 
				      work,
				      0) >= 0);

  for (i = 1; i < size; i ++) {
    ck_assert(isnormal(data[i]) || data[i] == 0.0);
  }

  ck_assert(within1pc(data[0], 3.14));

  for (i = 1; i < size; i ++) {
    ck_assert(within1pc(data[i], 0.0));
  }
  
  ck_assert(daub4_dwt_inverse2d_daub4(data, 
				      NS2_WIDTH,
				      NS2_HEIGHT,
				      NS2_WIDTH, 
				      work,
				      0) >= 0);
  
  for (i = 0; i < size; i ++) {
    ck_assert(within1pc(data[i], 3.14));
  }
}
END_TEST

START_TEST (test_daub4_2d_sinusoid_nonsquare1)
{
  double data[NS1_WIDTH * NS1_HEIGHT];
  double recon[NS1_WIDTH * NS1_HEIGHT];
  double work[NS1_WIDTH];
  int i;
  int j;

  int size;

  size = NS1_WIDTH * NS1_HEIGHT;

  for (j = 0; j < NS1_HEIGHT; j ++) {
    
    for (i = 0; i < NS1_WIDTH; i ++) {
      data[j * NS1_WIDTH + i] = sin((double)i/64.0 * M_PI * 2.0) * cos((double)j/8.0 * M_PI * 2.0);
      recon[j * NS1_WIDTH + i] = data[j * WIDTH + i];
    }

  }

  for (i = 0; i < size; i ++) {
    ck_assert(isnormal(recon[i]) || isfinite(recon[i]));
  }

  ck_assert(daub4_dwt_forward2d_daub4(recon, 
				      NS1_WIDTH,
				      NS1_HEIGHT,
				      NS1_WIDTH, 
				      work,
				      0) >= 0);

  for (i = 0; i < size; i ++) {
    ck_assert(isnormal(recon[i]) || isfinite(recon[i]));
  }
  
  ck_assert(daub4_dwt_inverse2d_daub4(recon, 
				      NS1_WIDTH,
				      NS1_HEIGHT,
				      NS1_WIDTH, 
				      work,
				      0) >= 0);

  for (i = 0; i < size; i ++) {
    ck_assert(isnormal(recon[i]) || isfinite(recon[i]));
  }

  for (i = 0; i < size; i ++) {
    ck_assert(within1pc(data[i], recon[i]));
  }
}
END_TEST

START_TEST (test_daub4_2d_sinusoid_nonsquare2)
{
  double data[NS2_WIDTH * NS2_HEIGHT];
  double recon[NS2_WIDTH * NS2_HEIGHT];
  double work[NS2_HEIGHT];
  int i;
  int j;

  int size;

  size = NS2_WIDTH * NS2_HEIGHT;

  for (j = 0; j < NS2_HEIGHT; j ++) {
    
    for (i = 0; i < NS2_WIDTH; i ++) {
      data[j * NS2_WIDTH + i] = sin((double)i/64.0 * M_PI * 2.0) * cos((double)j/8.0 * M_PI * 2.0);
      recon[j * NS2_WIDTH + i] = data[j * NS2_WIDTH + i];
    }

  }

  for (i = 0; i < size; i ++) {
    ck_assert(isnormal(recon[i]) || isfinite(recon[i]));
  }

  ck_assert(daub4_dwt_forward2d_daub4(recon, 
				      NS2_WIDTH,
				      NS2_HEIGHT,
				      NS2_WIDTH, 
				      work,
				      0) >= 0);

  for (i = 0; i < size; i ++) {
    ck_assert(isnormal(recon[i]) || isfinite(recon[i]));
  }
  
  ck_assert(daub4_dwt_inverse2d_daub4(recon, 
				      NS2_WIDTH,
				      NS2_HEIGHT,
				      NS2_WIDTH, 
				      work,
				      0) >= 0);

  for (i = 0; i < size; i ++) {
    ck_assert(isnormal(recon[i]) || isfinite(recon[i]));
  }

  for (i = 0; i < size; i ++) {
    ck_assert(within1pc(data[i], recon[i]));
  }
}
END_TEST

START_TEST (test_daub4_3d_sinusoid)
{
  double data[WIDTH * WIDTH * WIDTH];
  double recon[WIDTH * WIDTH * WIDTH];
  double work[WIDTH];
  int i;
  int j;
  int k;

  for (k = 0; k < WIDTH; k ++) {
    for (j = 0; j < WIDTH; j ++) {
      for (i = 0; i < WIDTH; i ++) {
	data[(k * WIDTH + j) * WIDTH + i] = 
	  sin((double)k/32.0 * M_PI * 2.0) * sin((double)i/64.0 * M_PI * 2.0) * cos((double)j/8.0 * M_PI * 2.0);
	recon[(k * WIDTH + j) * WIDTH + i] = data[(k * WIDTH + j) * WIDTH + i];
      }
    }
  }


  for (i = 0; i < (WIDTH * WIDTH * WIDTH); i ++) {
    ck_assert(isnormal(recon[i]) || isfinite(recon[i]));
  }

  /*
   * First check a single step
   */
  ck_assert(daub4_dwt_forward3d_daub4_step(recon,
					   WIDTH,
					   WIDTH,
					   WIDTH,
					   WIDTH,
					   WIDTH * WIDTH,
					   work) >= 0);
  
  for (i = 0; i < (WIDTH * WIDTH * WIDTH); i ++) {
    ck_assert(isnormal(recon[i]) || isfinite(recon[i]));
  }

  ck_assert(daub4_dwt_inverse3d_daub4_step(recon, 
					   WIDTH,
					   WIDTH,
					   WIDTH, 
					   WIDTH,
					   WIDTH * WIDTH,
					   work) >= 0);
  
  for (i = 0; i < (WIDTH * WIDTH * WIDTH); i ++) {
    ck_assert(isnormal(recon[i]) || isfinite(recon[i]));
  }

  for (i = 0; i < (WIDTH * WIDTH * WIDTH); i ++) {
    ck_assert(within1pc(data[i], recon[i]));
  }

  /*
   * Check a full multiresolution decomposition.
   */
  for (i = 0; i < (WIDTH * WIDTH * WIDTH); i ++) {
    recon[i] = data[i];
  }

  ck_assert(daub4_dwt_forward3d_daub4(recon, 
				      WIDTH,
				      WIDTH,
				      WIDTH, 
				      WIDTH,
				      WIDTH * WIDTH,
				      work,
				      0) >= 0);

  for (i = 0; i < (WIDTH * WIDTH * WIDTH); i ++) {
    ck_assert(isnormal(recon[i]) || isfinite(recon[i]));
  }
  
  ck_assert(daub4_dwt_inverse3d_daub4(recon, 
				      WIDTH,
				      WIDTH,
				      WIDTH, 
				      WIDTH,
				      WIDTH * WIDTH,
				      work,
				      0) >= 0);

  for (i = 0; i < (WIDTH * WIDTH * WIDTH); i ++) {
    ck_assert(isnormal(recon[i]) || isfinite(recon[i]));
  }

  for (i = 0; i < (WIDTH * WIDTH * WIDTH); i ++) {
    ck_assert(within1pc(data[i], recon[i]));
  }
}
END_TEST
 
START_TEST (test_daub4_3d_sinusoid_nonsquare1)
{
  double data[NS3_WIDTH * NS3_HEIGHT * NS3_DEPTH];
  double recon[NS3_WIDTH * NS3_HEIGHT * NS3_DEPTH];
  double work[NS3_DEPTH];
  int i;
  int j;
  int k;
  int size;

  for (k = 0; k < NS3_DEPTH; k ++) {
    for (j = 0; j < NS3_HEIGHT; j ++) {
      for (i = 0; i < NS3_WIDTH; i ++) {
	data[(k * NS3_HEIGHT + j) * NS3_WIDTH + i] = 
	  sin((double)k/32.0 * M_PI * 2.0) * sin((double)i/64.0 * M_PI * 2.0) * cos((double)j/8.0 * M_PI * 2.0);
      }
    }
  }

  size = NS3_WIDTH * NS3_HEIGHT * NS3_DEPTH;

  for (i = 0; i < size; i ++) {
    recon[i] = data[i];
    ck_assert(isnormal(recon[i]) || isfinite(recon[i]));
  }

  /*
   * First check a single step
   */
  ck_assert(daub4_dwt_forward3d_daub4_step(recon,
					    NS3_WIDTH,
					    NS3_HEIGHT,
					    NS3_DEPTH, 
					    NS3_WIDTH,
					    NS3_WIDTH * NS3_HEIGHT,
					    work) >= 0);

  for (i = 0; i < size; i ++) {
    ck_assert(isnormal(recon[i]) || isfinite(recon[i]));
  }

  ck_assert(daub4_dwt_inverse3d_daub4_step(recon, 
					    NS3_WIDTH,
					    NS3_HEIGHT,
					    NS3_DEPTH, 
					    NS3_WIDTH,
					    NS3_WIDTH * NS3_HEIGHT,
					    work) >= 0);
  
  for (i = 0; i < size; i ++) {
    ck_assert(isnormal(recon[i]) || isfinite(recon[i]));
  }

  for (i = 0; i < size; i ++) {
    ck_assert(within1pc(data[i], recon[i]));
  }

  /*
   * Check a full multiresolution decomposition.
   */
  for (i = 0; i < size; i ++) {
    recon[i] = data[i];
  }

  ck_assert(daub4_dwt_forward3d_daub4(recon, 
				      NS3_WIDTH,
				      NS3_HEIGHT,
				      NS3_DEPTH, 
				      NS3_WIDTH,
				      NS3_WIDTH * NS3_HEIGHT,
				      work,
				      0) >= 0);
  
  for (i = 0; i < size; i ++) {
    ck_assert(isnormal(recon[i]) || isfinite(recon[i]));
  }
  
  ck_assert(daub4_dwt_inverse3d_daub4(recon, 
				      NS3_WIDTH,
				      NS3_HEIGHT,
				      NS3_DEPTH, 
				      NS3_WIDTH,
				      NS3_WIDTH * NS3_HEIGHT,
				      work,
				      0) >= 0);

  for (i = 0; i < size; i ++) {
    ck_assert(isnormal(recon[i]) || isfinite(recon[i]));
  }

  for (i = 0; i < size; i ++) {
    ck_assert(within1pc(data[i], recon[i]));
  }
}
END_TEST

START_TEST (test_daub4_3d_sinusoid_nonsquare2)
{
  double data[NS4_WIDTH * NS4_HEIGHT * NS4_DEPTH];
  double recon[NS4_WIDTH * NS4_HEIGHT * NS4_DEPTH];
  double work[NS4_HEIGHT];
  int i;
  int j;
  int k;
  int size;

  for (k = 0; k < NS4_DEPTH; k ++) {
    for (j = 0; j < NS4_HEIGHT; j ++) {
      for (i = 0; i < NS4_WIDTH; i ++) {
	data[(k * NS4_HEIGHT + j) * NS4_WIDTH + i] = 
	  sin((double)k/32.0 * M_PI * 2.0) * sin((double)i/64.0 * M_PI * 2.0) * cos((double)j/8.0 * M_PI * 2.0);
      }
    }
  }

  size = NS4_WIDTH * NS4_HEIGHT * NS4_DEPTH;

  for (i = 0; i < size; i ++) {
    recon[i] = data[i];
    ck_assert(isnormal(recon[i]) || isfinite(recon[i]));
  }

  /*
   * First check a single step
   */
  ck_assert(daub4_dwt_forward3d_daub4_step(recon,
					  NS4_WIDTH,
					  NS4_HEIGHT,
					  NS4_DEPTH, 
					  NS4_WIDTH,
					  NS4_WIDTH * NS4_HEIGHT,
					  work) >= 0);

  for (i = 0; i < size; i ++) {
    ck_assert(isnormal(recon[i]) || isfinite(recon[i]));
  }

  ck_assert(daub4_dwt_inverse3d_daub4_step(recon, 
					    NS4_WIDTH,
					    NS4_HEIGHT,
					    NS4_DEPTH, 
					    NS4_WIDTH,
					    NS4_WIDTH * NS4_HEIGHT,
					    work) >= 0);
  
  for (i = 0; i < size; i ++) {
    ck_assert(isnormal(recon[i]) || isfinite(recon[i]));
  }

  for (i = 0; i < size; i ++) {
    ck_assert(within1pc(data[i], recon[i]));
  }

  /*
   * Check a full multiresolution decomposition.
   */
  for (i = 0; i < size; i ++) {
    recon[i] = data[i];
  }

  ck_assert(daub4_dwt_forward3d_daub4(recon, 
				      NS4_WIDTH,
				      NS4_HEIGHT,
				      NS4_DEPTH, 
				      NS4_WIDTH,
				      NS4_WIDTH * NS4_HEIGHT,
				      work,
				      0) >= 0);

  for (i = 0; i < size; i ++) {
    ck_assert(isnormal(recon[i]) || isfinite(recon[i]));
  }
  
  ck_assert(daub4_dwt_inverse3d_daub4(recon, 
				      NS4_WIDTH,
				      NS4_HEIGHT,
				      NS4_DEPTH, 
				      NS4_WIDTH,
				      NS4_WIDTH * NS4_HEIGHT,
				      work,
				      0) >= 0);

  for (i = 0; i < size; i ++) {
    ck_assert(isnormal(recon[i]) || isfinite(recon[i]));
  }

  for (i = 0; i < size; i ++) {
    ck_assert(within1pc(data[i], recon[i]));
  }
}
END_TEST

Suite *
daub4_suite (void)
{
  Suite *s = suite_create ("DAUB4");

  /* Core test case */
  TCase *tc_core = tcase_create ("Core");

  tcase_add_test (tc_core, test_daub4_1d_simple_step);
  
  tcase_add_test (tc_core, test_daub4_1d_constant_step);
  tcase_add_test (tc_core, test_daub4_1d_sinusoid_step);

  tcase_add_test (tc_core, test_daub4_1d_constant);
  tcase_add_test (tc_core, test_daub4_1d_sinusoid);
  
  tcase_add_test (tc_core, test_daub4_2d_constant);
  tcase_add_test (tc_core, test_daub4_2d_sinusoid);
  tcase_add_test (tc_core, test_daub4_2d_constant_nonsquare1);
  tcase_add_test (tc_core, test_daub4_2d_constant_nonsquare2);
  tcase_add_test (tc_core, test_daub4_2d_sinusoid_nonsquare1);
  tcase_add_test (tc_core, test_daub4_2d_sinusoid_nonsquare2);
  tcase_add_test (tc_core, test_daub4_3d_sinusoid);
  tcase_add_test (tc_core, test_daub4_3d_sinusoid_nonsquare1);
  tcase_add_test (tc_core, test_daub4_3d_sinusoid_nonsquare2);

  suite_add_tcase (s, tc_core);

  return s;
}

int main (void) 
{
  int number_failed;
  Suite *s = daub4_suite ();
  SRunner *sr = srunner_create (s);

  srunner_set_fork_status (sr, CK_NOFORK);

  srunner_run_all (sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
