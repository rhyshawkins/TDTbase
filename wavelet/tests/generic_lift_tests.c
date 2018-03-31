
#include <stdio.h>
#include <stdlib.h>
#include <check.h>
#include <math.h>

#include "cdf97_lift.h"
#include "generic_lift.h"
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


START_TEST (test_generic_2d_sinusoid_nonsquare1)
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

  ck_assert(generic_lift_forward2d(recon, 
				   NS1_WIDTH,
				   NS1_HEIGHT,
				   NS1_WIDTH, 
				   work,
				   cdf97_lift_forward1d_cdf97_step,
				   cdf97_lift_forward1d_cdf97_step,
				   1) >= 0);

  for (i = 0; i < size; i ++) {
    ck_assert(isnormal(recon[i]) || isfinite(recon[i]));
  }
  
  ck_assert(generic_lift_inverse2d(recon, 
				   NS1_WIDTH,
				   NS1_HEIGHT,
				   NS1_WIDTH, 
				   work,
				   cdf97_lift_inverse1d_cdf97_step,
				   cdf97_lift_inverse1d_cdf97_step,
				   1) >= 0);

  for (i = 0; i < size; i ++) {
    ck_assert(isnormal(recon[i]) || isfinite(recon[i]));
  }

  for (i = 0; i < size; i ++) {
    ck_assert(within1pc(data[i], recon[i]));
  }
}
END_TEST

START_TEST (test_cdf97_superresolution)
{
  double data[NS1_WIDTH * NS1_HEIGHT];
  double work[NS1_WIDTH];

  int i;
  int j;
  int size;

  double dy;
  double dx;
  
  /*
   * Constant image
   */
  size = NS1_WIDTH * NS1_HEIGHT;
  memset(data, 0, sizeof(double) * NS1_WIDTH * NS1_HEIGHT);
  for (j = 0; j < NS1_HEIGHT/2; j ++) {
    for (i = 0; i < NS1_WIDTH/2; i ++) {
      data[j * NS1_WIDTH + i] = 1.0;
    }
  }

  ck_assert(generic_lift_inverse2d_step(data, NS1_WIDTH, NS1_HEIGHT, NS1_WIDTH, work,
					cdf97_lift_inverse1d_cdf97_step,
					cdf97_lift_inverse1d_cdf97_step) >= 0);

  for (i = 0; i < size; i ++) {
    ck_assert(isnormal(data[i]) || isfinite(data[i]));
  }

  for (j = 0; j < NS1_HEIGHT; j ++) {
    for (i = 0; i < NS1_WIDTH; i ++) {
      printf("%3.1f ", data[j * NS1_WIDTH + i]);
    }
    printf("\n");
  }
  
  for (i = 0; i < size; i ++) {
    ck_assert(within1pc(data[i], 1.0));
  }

  /*
   * Gaussian blob image
   */
  size = NS1_WIDTH * NS1_HEIGHT;
  memset(data, 0, sizeof(double) * NS1_WIDTH * NS1_HEIGHT);
  for (j = 0; j < NS1_HEIGHT/2; j ++) {
    dy = (double)j - (double)NS1_HEIGHT/4;
    
    for (i = 0; i < NS1_WIDTH/2; i ++) {
      dx = (double)i - (double)NS1_WIDTH/4;

      
      data[j * NS1_WIDTH + i] = exp(-(dx*dx + dy*dy));

      printf("%3.1f ", data[j * NS1_WIDTH + i]);
    }
    printf("\n");
  }

  ck_assert(generic_lift_inverse2d_step(data, NS1_WIDTH, NS1_HEIGHT, NS1_WIDTH, work,
					cdf97_lift_inverse1d_cdf97_step,
					cdf97_lift_inverse1d_cdf97_step) >= 0);

  for (j = 0; j < NS1_HEIGHT; j ++) {
    for (i = 0; i < NS1_WIDTH; i ++) {
      printf("%3.1f ", data[j * NS1_WIDTH + i]);
    }
    printf("\n");
  }
  
  for (i = 0; i < size; i ++) {
    ck_assert(isnormal(data[i]) || isfinite(data[i]));
  }

  
}
END_TEST

#if 0
START_TEST (test_cdf97_2d_sinusoid_nonsquare2)
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

  ck_assert(cdf97_lift_forward2d_cdf97(recon, 
				       NS2_WIDTH,
				       NS2_HEIGHT,
				       NS2_WIDTH, 
				       work,
				       0) >= 0);

  for (i = 0; i < size; i ++) {
    ck_assert(isnormal(recon[i]) || isfinite(recon[i]));
  }
  
  ck_assert(cdf97_lift_inverse2d_cdf97(recon, 
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

START_TEST (test_cdf97_3d_sinusoid)
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
  ck_assert(cdf97_lift_forward3d_cdf97_step(recon,
					    WIDTH,
					    WIDTH,
					    WIDTH,
					    WIDTH,
					    WIDTH * WIDTH,
					    work) >= 0);

  for (i = 0; i < (WIDTH * WIDTH * WIDTH); i ++) {
    ck_assert(isnormal(recon[i]) || isfinite(recon[i]));
  }

  ck_assert(cdf97_lift_inverse3d_cdf97_step(recon, 
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

  ck_assert(cdf97_lift_forward3d_cdf97(recon, 
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
  
  ck_assert(cdf97_lift_inverse3d_cdf97(recon, 
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
 
START_TEST (test_cdf97_3d_sinusoid_nonsquare1)
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
  ck_assert(cdf97_lift_forward3d_cdf97_step(recon,
					    NS3_WIDTH,
					    NS3_HEIGHT,
					    NS3_DEPTH, 
					    NS3_WIDTH,
					    NS3_WIDTH * NS3_HEIGHT,
					    work) >= 0);

  for (i = 0; i < size; i ++) {
    ck_assert(isnormal(recon[i]) || isfinite(recon[i]));
  }

  ck_assert(cdf97_lift_inverse3d_cdf97_step(recon, 
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

  ck_assert(cdf97_lift_forward3d_cdf97(recon, 
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
  
  ck_assert(cdf97_lift_inverse3d_cdf97(recon, 
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

START_TEST (test_cdf97_3d_sinusoid_nonsquare2)
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
  ck_assert(cdf97_lift_forward3d_cdf97_step(recon,
					    NS4_WIDTH,
					    NS4_HEIGHT,
					    NS4_DEPTH, 
					    NS4_WIDTH,
					    NS4_WIDTH * NS4_HEIGHT,
					    work) >= 0);

  for (i = 0; i < size; i ++) {
    ck_assert(isnormal(recon[i]) || isfinite(recon[i]));
  }

  ck_assert(cdf97_lift_inverse3d_cdf97_step(recon, 
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

  ck_assert(cdf97_lift_forward3d_cdf97(recon, 
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
  
  ck_assert(cdf97_lift_inverse3d_cdf97(recon, 
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

START_TEST (test_cdf97_3d_sinusoid_nonsquare2_subtile)
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
   * Check a full multiresolution decomposition.
   */
  ck_assert(cdf97_lift_forward3d_cdf97(recon, 
				       NS4_WIDTH,
				       NS4_HEIGHT,
				       NS4_DEPTH, 
				       NS4_WIDTH,
				       NS4_WIDTH * NS4_HEIGHT,
				       work,
				       1) >= 0);

  for (i = 0; i < size; i ++) {
    ck_assert(isnormal(recon[i]) || isfinite(recon[i]));
  }
  
  ck_assert(cdf97_lift_inverse3d_cdf97(recon, 
				       NS4_WIDTH,
				       NS4_HEIGHT,
				       NS4_DEPTH, 
				       NS4_WIDTH,
				       NS4_WIDTH * NS4_HEIGHT,
				       work,
				       1) >= 0);
  
  for (i = 0; i < size; i ++) {
    ck_assert(isnormal(recon[i]) || isfinite(recon[i]));
  }

  for (i = 0; i < size; i ++) {
    ck_assert(within1pc(data[i], recon[i]));
  }
}
END_TEST

START_TEST (test_generic_cdf97_3d_sinusoid)
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
  ck_assert(generic_lift_forward3d(recon,
				   WIDTH,
				   WIDTH,
				   WIDTH,
				   WIDTH,
				   WIDTH * WIDTH,
				   work,
				   cdf97_lift_forward1d_cdf97_step,
				   cdf97_lift_forward1d_cdf97_step,
				   cdf97_lift_forward1d_cdf97_step,
				   1) >= 0);

  for (i = 0; i < (WIDTH * WIDTH * WIDTH); i ++) {
    ck_assert(isnormal(recon[i]) || isfinite(recon[i]));
  }

  ck_assert(generic_lift_inverse3d(recon, 
				   WIDTH,
				   WIDTH,
				   WIDTH, 
				   WIDTH,
				   WIDTH * WIDTH,
				   work,
				   cdf97_lift_inverse1d_cdf97_step,
				   cdf97_lift_inverse1d_cdf97_step,
				   cdf97_lift_inverse1d_cdf97_step,
				   1) >= 0);
  
  for (i = 0; i < (WIDTH * WIDTH * WIDTH); i ++) {
    ck_assert(isnormal(recon[i]) || isfinite(recon[i]));
  }

  for (i = 0; i < (WIDTH * WIDTH * WIDTH); i ++) {
    ck_assert(within1pc(data[i], recon[i]));
  }
}
END_TEST

START_TEST (test_generic_cdf97_3d_sinusoid_nonsquare2_subtile)
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
   * Check a full multiresolution decomposition.
   */
  ck_assert(generic_lift_forward3d(recon, 
				   NS4_WIDTH,
				   NS4_HEIGHT,
				   NS4_DEPTH, 
				   NS4_WIDTH,
				   NS4_WIDTH * NS4_HEIGHT,
				   work,
				   cdf97_lift_forward1d_cdf97_step,
				   cdf97_lift_forward1d_cdf97_step,
				   cdf97_lift_forward1d_cdf97_step,
				   1) >= 0);

  for (i = 0; i < size; i ++) {
    ck_assert(isnormal(recon[i]) || isfinite(recon[i]));
  }
  
  ck_assert(generic_lift_inverse3d(recon, 
				   NS4_WIDTH,
				   NS4_HEIGHT,
				   NS4_DEPTH, 
				   NS4_WIDTH,
				   NS4_WIDTH * NS4_HEIGHT,
				   work,
				   cdf97_lift_inverse1d_cdf97_step,
				   cdf97_lift_inverse1d_cdf97_step,
				   cdf97_lift_inverse1d_cdf97_step,
				   1) >= 0);
  
  for (i = 0; i < size; i ++) {
    ck_assert(isnormal(recon[i]) || isfinite(recon[i]));
  }

  for (i = 0; i < size; i ++) {
    /* if (!within1pc(data[i], recon[i])) { */
    /*   printf("%d %f %f\n", i, data[i], recon[i]); */
    /* } */
    ck_assert(within1pc(data[i], recon[i]));
  }
}
END_TEST
#endif

Suite *
cdf97_suite (void)
{
  Suite *s = suite_create ("CDF97");

  /* Core test case */
  TCase *tc_core = tcase_create ("Core");

  tcase_add_test (tc_core, test_generic_2d_sinusoid_nonsquare1);

  tcase_add_test (tc_core, test_cdf97_superresolution);

  /* tcase_add_test (tc_core, test_cdf97_2d_sinusoid_nonsquare2); */
  /* tcase_add_test (tc_core, test_cdf97_3d_sinusoid); */
  /* tcase_add_test (tc_core, test_cdf97_3d_sinusoid_nonsquare1); */
  /* tcase_add_test (tc_core, test_cdf97_3d_sinusoid_nonsquare2); */
  /* tcase_add_test (tc_core, test_cdf97_3d_sinusoid_nonsquare2_subtile); */
  /* tcase_add_test (tc_core, test_generic_cdf97_3d_sinusoid); */
  /* tcase_add_test (tc_core, test_generic_cdf97_3d_sinusoid_nonsquare2_subtile); */

  suite_add_tcase (s, tc_core);

  return s;
}

int main (void) 
{
  int number_failed;
  Suite *s = cdf97_suite ();
  SRunner *sr = srunner_create (s);

  srunner_set_fork_status (sr, CK_NOFORK);

  srunner_run_all (sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
