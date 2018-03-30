
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "coefficient_histogram.h"

#include "slog.h"

static int bin_index(double v, double vmin, double vmax, int nbins);

static double bin_center(const coefficient_histogram_t *c,
			 double vmin,
			 double vmax,
			 int i);

static double percent(int proposed, int accepted);


coefficient_histogram_t *
coefficient_histogram_create(int ncoeff,
			     int nbins,
			     double vmin,
			     double vmax,
			     ch_coord_to_index_t coordtoindex,
			     ch_index_to_coord_t indextocoord,
			     void *ch_user)
{
  coefficient_histogram_t *c;
  int i;

  c = malloc(sizeof(coefficient_histogram_t));
  if (c == NULL) {
    ERROR("failed to allocate struct");
    return NULL;
  }

  c->coordtoindex = coordtoindex;
  c->indextocoord = indextocoord;
  c->ch_user = ch_user;

  c->ncoeff = ncoeff;
  c->nbins = nbins;
  c->gvmin = vmin;
  c->gvmax = vmax;

  c->counts = malloc(sizeof(int*) * ncoeff);
  if (c->counts == NULL) {
    ERROR("failed to allocate hist array");
    return NULL;
  }

  for (i = 0; i < ncoeff; i ++) {
    c->counts[i] = malloc(sizeof(int) * nbins);
    if (c->counts[i] == NULL) {
      ERROR("failed to allocate hist array %d", i);
      return NULL;
    }
  }

#define ALLOC1D(var, type, n, name) \
  var = malloc(sizeof(type) * n); \
  if (var == NULL) { \
  ERROR("failed to allocate %s", name); \
  return NULL; \
  }

  ALLOC1D(c->vmin, double, ncoeff, "vmin");
  ALLOC1D(c->vmax, double, ncoeff, "vmax");

  ALLOC1D(c->under, int, ncoeff, "underflow");
  ALLOC1D(c->over, int, ncoeff, "overflow");
  ALLOC1D(c->rmin, double, ncoeff, "min");
  ALLOC1D(c->rmax, double, ncoeff, "max");
  ALLOC1D(c->rmean, double, ncoeff, "mean");
  ALLOC1D(c->rstd, double, ncoeff, "std. dev.");
  ALLOC1D(c->n, int, ncoeff, "count");

  ALLOC1D(c->valpha, double, ncoeff, "valpha");
  ALLOC1D(c->valpha_mean, double, ncoeff, "valpha mean");
  ALLOC1D(c->valpha_n, int, ncoeff, "valpha n");
  
  ALLOC1D(c->pb, int, ncoeff, "pb");
  ALLOC1D(c->ab, int, ncoeff, "ab");
  ALLOC1D(c->pd, int, ncoeff, "pd");
  ALLOC1D(c->ad, int, ncoeff, "ad");
  ALLOC1D(c->pv, int, ncoeff, "pv");
  ALLOC1D(c->av, int, ncoeff, "av");

  if (coefficient_histogram_reset(c) < 0) {
    ERROR("failed to reset values");
    return NULL;
  }

  return c;
}

void
coefficient_histogram_destroy(coefficient_histogram_t *c)
{
  int i;

  if (c != NULL) {

    free(c->av);
    free(c->pv);
    free(c->ad);
    free(c->pd);
    free(c->ab);
    free(c->pb);

    free(c->valpha_n);
    free(c->valpha_mean);
    free(c->valpha);

    free(c->n);
    free(c->rstd);
    free(c->rmean);
    free(c->rmax);
    free(c->rmin);
    free(c->over);
    free(c->under);

    free(c->vmin);
    free(c->vmax);

    for (i = 0; i < c->ncoeff; i ++) {
      free(c->counts[i]);
    }
    free(c->counts);
    free(c);
  }
}

int
coefficient_histogram_save(coefficient_histogram_t *c,
			   const char *filename)
{
  FILE *fp;
  int i;

  fp = fopen(filename, "w");
  if (fp == NULL) {
    ERROR("failed to create file");
    return -1;
  }

  fwrite(&c->ncoeff, sizeof(int), 1, fp);
  fwrite(&c->nbins, sizeof(int), 1, fp);

  fwrite(&c->gvmin, sizeof(double), 1, fp);
  fwrite(&c->gvmax, sizeof(double), 1, fp);

  fwrite(c->vmin, sizeof(double), c->ncoeff, fp);
  fwrite(c->vmax, sizeof(double), c->ncoeff, fp);

  for (i = 0; i < c->ncoeff; i ++) {
    fwrite(c->counts[i], sizeof(int), c->nbins, fp);
  }
  fwrite(c->under, sizeof(int), c->ncoeff, fp);
  fwrite(c->over, sizeof(int), c->ncoeff, fp);
  
  fwrite(c->rmin, sizeof(double), c->ncoeff, fp);
  fwrite(c->rmax, sizeof(double), c->ncoeff, fp);
  fwrite(c->rmean, sizeof(double), c->ncoeff, fp);
  fwrite(c->rstd, sizeof(double), c->ncoeff, fp);
  fwrite(c->n, sizeof(int), c->ncoeff, fp);

  fwrite(c->valpha, sizeof(double), c->ncoeff, fp);
  fwrite(c->valpha_mean, sizeof(double), c->ncoeff, fp);
  fwrite(c->valpha_n, sizeof(int), c->ncoeff, fp);

  fwrite(c->pb, sizeof(int), c->ncoeff, fp);
  fwrite(c->ab, sizeof(int), c->ncoeff, fp);
  fwrite(c->pd, sizeof(int), c->ncoeff, fp);
  fwrite(c->ad, sizeof(int), c->ncoeff, fp);
  fwrite(c->pv, sizeof(int), c->ncoeff, fp);
  fwrite(c->av, sizeof(int), c->ncoeff, fp);

  fclose(fp);

  return 0;
}

int
coefficient_histogram_load(coefficient_histogram_t *c,
			   const char *filename)
{
  FILE *fp;
  int i;
  int ncoeff;
  int nbins;

  fp = fopen(filename, "r");
  if (fp == NULL) {
    ERROR("failed to open file");
    return -1;
  }

  if (fread(&ncoeff, sizeof(int), 1, fp) != 1) {
    ERROR("failed to read ncoeff");
    return -1;
  }

  if (fread(&nbins, sizeof(int), 1, fp) != 1) {
    ERROR("failed to read nbins");
    return -1;
  }

  if (ncoeff != c->ncoeff ||
      nbins != c->nbins) {
    ERROR("size mismatch ncoeff %d != %d, nbins %d != %d",
	    ncoeff, c->ncoeff,
	    nbins, c->nbins);
    return -1;
  }

  if (fread(&c->gvmin, sizeof(double), 1, fp) != 1) {
    ERROR("failed to read gvmin");
    return -1;
  }
  
  if (fread(&c->gvmax, sizeof(double), 1, fp) != 1) {
    ERROR("failed to read gvmax");
    return -1;
  }

  if (fread(c->vmin, sizeof(double), c->ncoeff, fp) != c->ncoeff) {
    ERROR("failed to read vmin");
    return -1;
  }
  if (fread(c->vmax, sizeof(double), c->ncoeff, fp) != c->ncoeff) {
    ERROR("failed to read vmax");
    return -1;
  }

  for (i = 0; i < c->ncoeff; i ++) {
    if (fread(c->counts[i], sizeof(int), nbins, fp) != nbins) {
      ERROR("failed to read row of counts");
      return -1;
    }
  }
  
  if (fread(c->under, sizeof(int), ncoeff, fp) != c->ncoeff) {
    ERROR("failed to read under");
    return -1;
  }
  if (fread(c->over, sizeof(int), ncoeff, fp) != c->ncoeff) {
    ERROR("failed to read over");
    return -1;
  }
  
  if (fread(c->rmin, sizeof(double), c->ncoeff, fp) != c->ncoeff) {
    ERROR("failed to read rmin");
    return -1;
  }
  if (fread(c->rmax, sizeof(double), c->ncoeff, fp) != c->ncoeff) {
    ERROR("failed to read rmax");
    return -1;
  }
  if (fread(c->rmean, sizeof(double), c->ncoeff, fp) != c->ncoeff) {
    ERROR("failed to read rmean");
    return -1;
  }
  if (fread(c->rstd, sizeof(double), c->ncoeff, fp) != c->ncoeff) {
    ERROR("failed to read rstd");
    return -1;
  }
  if (fread(c->n, sizeof(int), c->ncoeff, fp) != c->ncoeff) {
    ERROR("failed to read n");
    return -1;
  }

  if (fread(c->valpha, sizeof(double), c->ncoeff, fp) != c->ncoeff) {
    ERROR("failed to read valpha");
    return -1;
  }
  if (fread(c->valpha_mean, sizeof(double), c->ncoeff, fp) != c->ncoeff) {
    ERROR("failed to read valpha_mean");
    return -1;
  }
  if (fread(c->valpha_n, sizeof(int), c->ncoeff, fp) != c->ncoeff) {
    ERROR("failed to read valpha_n");
    return -1;
  }

  if (fread(c->pb, sizeof(int), c->ncoeff, fp) != c->ncoeff) {
    ERROR("failed to read pb");
    return -1;
  }
  if (fread(c->ab, sizeof(int), c->ncoeff, fp) != c->ncoeff) {
    ERROR("failed to read ab");
    return -1;
  }
  if (fread(c->pd, sizeof(int), c->ncoeff, fp) != c->ncoeff) {
    ERROR("failed to read pd");
    return -1;
  }
  if (fread(c->ad, sizeof(int), c->ncoeff, fp) != c->ncoeff) {
    ERROR("failed to read ad");
    return -1;
  }
  if (fread(c->pv, sizeof(int), c->ncoeff, fp) != c->ncoeff) {
    ERROR("failed to read pv");
    return -1;
  }
  if (fread(c->av, sizeof(int), c->ncoeff, fp) != c->ncoeff) {
    ERROR("failed to read av");
    return -1;
  }

  fclose(fp);

  return 0;
}


int
coefficient_histogram_coord_to_index(const coefficient_histogram_t *c, int i, int j, int k, int depth)
{
  return c->coordtoindex(c->ch_user, i, j, k, depth);
}

int
coefficient_histogram_index_to_coord(const coefficient_histogram_t *c,
				     int index, int *i, int *j, int *k, int *depth)
{
  return c->indextocoord(c->ch_user, index, i, j, k, depth);
}

int
coefficient_histogram_reset(coefficient_histogram_t *c)
{
  int i;

  if (c == NULL) {
    return -1;
  }

  for (i = 0; i < c->ncoeff; i ++) {
    memset(c->counts[i], 0, sizeof(int) * c->nbins);

    c->vmin[i] = c->gvmin;
    c->vmax[i] = c->gvmax;
  }


  memset(c->under, 0, sizeof(int) * c->ncoeff);
  memset(c->over, 0, sizeof(int) * c->ncoeff);
  
  memset(c->rmin, 0, sizeof(double) * c->ncoeff);
  memset(c->rmax, 0, sizeof(double) * c->ncoeff);
  memset(c->rmean, 0, sizeof(double) * c->ncoeff);
  memset(c->rstd, 0, sizeof(double) * c->ncoeff);

  memset(c->n, 0, sizeof(int) * c->ncoeff);

  memset(c->valpha, 0, sizeof(double) * c->ncoeff);
  memset(c->valpha_mean, 0, sizeof(double) * c->ncoeff);
  memset(c->valpha_n, 0, sizeof(int) * c->ncoeff);
  
  memset(c->pb, 0, sizeof(int) * c->ncoeff);
  memset(c->ab, 0, sizeof(int) * c->ncoeff);
  memset(c->pd, 0, sizeof(int) * c->ncoeff);
  memset(c->ad, 0, sizeof(int) * c->ncoeff);
  memset(c->pv, 0, sizeof(int) * c->ncoeff);
  memset(c->av, 0, sizeof(int) * c->ncoeff);

  return 0;
}

int 
coefficient_histogram_set_range(coefficient_histogram_t *c, 
				int index,
				double vmin,
				double vmax)
{
  if (c == NULL ||
      index < 0 || index >= c->ncoeff) {
    ERROR("invalid parameters");
    return -1;
  }

  c->vmin[index] = vmin;
  c->vmax[index] = vmax;

  return 0;
}

int
coefficient_histogram_sample(coefficient_histogram_t *c, int index, double value)
{
  double delta;

  if (c == NULL ||
      index < 0 || index >= c->ncoeff) {
    ERROR("invalid parameters");
    return -1;
  }

  if (value < c->vmin[index]) {
    c->under[index] ++;
  } else if (value > c->vmax[index]) {
    c->over[index] ++;
  } else {
    c->counts[index][bin_index(value, c->vmin[index], c->vmax[index], c->nbins)] ++;
  }

  if (c->n[index] == 0) {
    /* Initialisation of range */
    c->rmin[index] = value;
    c->rmax[index] = value;
  } else {
    if (value < c->rmin[index]) {
      c->rmin[index] = value;
    }

    if (value > c->rmax[index]) {
      c->rmax[index] = value;
    }
  }

  c->n[index] ++;
  delta = value - c->rmean[index];
  c->rmean[index] += delta/(double)(c->n[index]);

  /* Note that this is computing variance and we need to sqrt at the 
   * end to get std. dev.
   */
  c->rstd[index] += delta*(value - c->rmean[index]);

  return 0;
}

int 
coefficient_histogram_propose_birth(coefficient_histogram_t *c, int index)
{
  if (c == NULL ||
      index < 0 || index >= c->ncoeff) {
    ERROR("invalid parameters (%p, %d, %d)", c, index, c == NULL ? -1 : c->ncoeff);
    return -1;
  }
  
  c->pb[index] ++;
  return 0;
}

int 
coefficient_histogram_accept_birth(coefficient_histogram_t *c, int index, double value)
{
  if (c == NULL ||
      index < 0 || index >= c->ncoeff) {
    ERROR("invalid parameters");
    return -1;
  }
  
  c->ab[index] ++;
  return 0;
}
			     
int 
coefficient_histogram_reject_birth(coefficient_histogram_t *c, int index, double value)
{
  if (c == NULL ||
      index < 0 || index >= c->ncoeff) {
    ERROR("invalid parameters %p %d",
	  c,
	  index);
    return -1;
  }
  
  /* TODO: record rejected birth values */
  return 0;
}

/*
 * For recording per coefficient death information
 */

int 
coefficient_histogram_propose_death(coefficient_histogram_t *c, int index)
{
  if (c == NULL ||
      index < 0 || index >= c->ncoeff) {
    ERROR("invalid parameters");
    return -1;
  }
  
  c->pd[index] ++;
  return 0;
}

int 
coefficient_histogram_accept_death(coefficient_histogram_t *c, int index)
{
  if (c == NULL ||
      index < 0 || index >= c->ncoeff) {
    ERROR("invalid parameters");
    return -1;
  }
  
  c->ad[index] ++;
  return 0;
}
			     
/*
 * For recording per coefficient value information
 */
int 
coefficient_histogram_propose_value(coefficient_histogram_t *c, int index)
{
  if (c == NULL ||
      index < 0 || index >= c->ncoeff) {
    ERROR("invalid parameters");
    return -1;
  }
  
  c->pv[index] ++;
  return 0;
}

int 
coefficient_histogram_accept_value(coefficient_histogram_t *c, int index, double value)
{
  if (c == NULL ||
      index < 0 || index >= c->ncoeff) {
    ERROR("coefficient_histogram_accept_value: invalid parameters");
    return -1;
  }
  
  c->av[index] ++;
  return 0;
}
			     
int 
coefficient_histogram_reject_value(coefficient_histogram_t *c, int index, double value)
{
  if (c == NULL ||
      index < 0 || index >= c->ncoeff) {
    ERROR("invalid parameters");
    return -1;
  }

  /* TODO: record rejected values */
  return 0;
}

int
coefficient_histogram_sample_value_alpha(coefficient_histogram_t *c, int index, double alpha)
{
  double delta;
  
  if (c == NULL ||
      index < 0 || index >= c->ncoeff) {
    ERROR("invalid parameters");
    return -1;
  }

  /*
   * Force alpha = log(max(1, ...)) here
   */
  if (alpha > 0.0) {
    alpha = 0.0;
  }
  
  c->valpha[index] = alpha;

  c->valpha_n[index] ++;

  delta = alpha - c->valpha_mean[index];
  c->valpha_mean[index] += delta/(double)(c->valpha_n[index]);
  
  return 0;
}

int 
coefficient_histogram_get_coefficient_mean_std(coefficient_histogram_t *c,
					       int index,
					       double *mean,
					       double *std)
{
  if (c == NULL || 
      index < 0 || index >= c->ncoeff) {
    ERROR("invalid parameters");
    return -1;
  }
  
  if (c->n[index] > 2) {

    *mean = c->rmean[index];
    *std = sqrt(c->rstd[index]/(double)(c->n[index] - 1));
    return c->n[index];

  } else {

    *mean = 0.0;
    *std = 0.0;
    return 0;

  }
}

int
coefficient_histogram_get_accept_reject(coefficient_histogram_t *c,
					int index,
					int *propose,
					int *accept)
{
  if (c == NULL || 
      index < 0 || index >= c->ncoeff) {
    ERROR("invalid parameters");
    return -1;
  }

  *propose = c->pv[index];
  *accept = c->av[index];

  return 0;
}


int
coefficient_histogram_finalise(coefficient_histogram_t *c)
{
  int i;

  if (c == NULL) {
    return -1;
  }

  for (i = 0; i < c->ncoeff; i ++) {
    if (c->n[i] > 2) {
      c->rstd[i] = sqrt(c->rstd[i]/(double)(c->n[i] - 1));
    } else {
      c->rstd[i] = 0.0;
    }
  }
  
  return 0;
}

int
coefficient_histogram_save_2D(const coefficient_histogram_t *c,
			      const char *filename)
{
  FILE *fp;
  int i;
  int j;

  fp = fopen(filename, "w");
  if (fp == NULL) {
    ERROR("failed to create file");
    return -1;
  }

  for (j = 0; j < c->ncoeff; j ++) {
    for (i = 0; i < c->nbins; i ++) {
      fprintf(fp, "%d ", c->counts[j][i]);
    }
    fprintf(fp, "\n");
  }

  fclose(fp);
  return 0;
}

static int save_summary_stat(const char *filename,
			     const double *stat,
			     int n)
{
  FILE *fp;
  int i;

  fp = fopen(filename, "w");
  if (fp == NULL) {
    return -1;
  }

  for (i = 0; i < n; i ++) {
    fprintf(fp, "%d %g\n", i, stat[i]);
  }

  fclose(fp);
  return 0;
}

int
coefficient_histogram_save_min(const coefficient_histogram_t *c,
			       const char *filename)
{
  return save_summary_stat(filename, c->rmin, c->ncoeff);
}


int
coefficient_histogram_save_max(const coefficient_histogram_t *c,
			       const char *filename)
{
  return save_summary_stat(filename, c->rmax, c->ncoeff);
}



int
coefficient_histogram_save_mean(const coefficient_histogram_t *c,
				const char *filename)
{
  return save_summary_stat(filename, c->rmean, c->ncoeff);
}


int
coefficient_histogram_save_std(const coefficient_histogram_t *c,
			       const char *filename)
{
  return save_summary_stat(filename, c->rstd, c->ncoeff);
}

int
coefficient_histogram_save_aggregated_histogram(const coefficient_histogram_t *c,
						const char *filename,
						coefficient_histogram_index_filter_t filter,
						void *user,
						int subset)
{
  int *hist;
  FILE *fp;
  int i;
  int j;
  double avmin;
  double avmax;

  fp = fopen(filename, "w");
  if (fp == NULL) {
    ERROR("failed to create file");
    return -1;
  }

  hist = malloc(sizeof(int) * c->nbins);
  if (hist == NULL) {
    ERROR("failed to allocate temporary hist");
    return -1;
  }
  memset(hist, 0, sizeof(int) * c->nbins);

  avmin = 0.0;
  avmax = 0.0;
  for (i = 0; i < c->ncoeff; i ++) {
    if (filter(user, subset, i)) {

      if (avmin == avmax) {
	avmin = c->vmin[i];
	avmax = c->vmax[i];
	printf("  av: %f %f\n", avmin, avmax);
      } else {
	if (avmin != c->vmin[i] ||
	    avmax != c->vmax[i]) {
	  ERROR("different histogram ranges");
	  return -1;
	}
      }

      for (j = 0; j < c->nbins; j ++) {
	hist[j] += c->counts[i][j];
      }
    }
  }

  for (i = 0; i < c->nbins; i ++) {
    fprintf(fp, "%f %d\n", bin_center(c, avmin, avmax, i), hist[i]);
  }

  free(hist);
  fclose(fp);

  return 0;
}

int
coefficient_histogram_save_aggregated_histogram_image(const coefficient_histogram_t *c,
						      const char *filename,
						      coefficient_histogram_index_filter_t filter,
						      void *user,
						      int subset)
{
  FILE *fp;
  int i;
  int j;
  double avmin;
  double avmax;

  fp = fopen(filename, "w");
  if (fp == NULL) {
    ERROR("failed to create file");
    return -1;
  }

  avmin = 0.0;
  avmax = 0.0;

  for (i = 0; i < c->ncoeff; i ++) {
    if (filter(user, subset, i)) {

      if (avmin == avmax) {
	avmin = c->vmin[i];
	avmax = c->vmax[i];

	/*
	 * Make the first line the bin centres
	 */
	for (j = 0; j < c->nbins; j ++) {
	  fprintf(fp, "%f ", bin_center(c, avmin, avmax, j));
	}
	fprintf(fp, "\n");

      } else {
	if (avmin != c->vmin[i] ||
	    avmax != c->vmax[i]) {
	  ERROR("coefficient_histogram_save_aggregated_histogram_image: different histogram ranges");
	  return -1;
	}
      }

      for (j = 0; j < c->nbins; j ++) {
	fprintf(fp, "%d ", c->counts[i][j]);
      }
      fprintf(fp, "\n");
    }
  }

  fclose(fp);

  return 0;
}

int
coefficient_histogram_save_acceptance(const coefficient_histogram_t *c,
				      const char *filename)
{
  FILE *fp;
  int i;
  

  fp = fopen(filename, "w");
  if (fp == NULL) {
    ERROR("failed to create file");
    return -1;
  }

  for (i = 0; i < c->ncoeff; i ++) {
    fprintf(fp, "%d %d %g %d %d %g %d %d %g\n",
	    c->pb[i], c->ab[i], percent(c->pb[i], c->ab[i]),
	    c->pd[i], c->ad[i], percent(c->pd[i], c->ad[i]),
	    c->pv[i], c->av[i], percent(c->pv[i], c->av[i]));
  }

  fclose(fp);
  return 0;
}

static int bin_index(double v, double vmin, double vmax, int nbins)
{
  return (int)((v - vmin)/(vmax - vmin) * (double)nbins);
}

static double bin_center(const coefficient_histogram_t *c,
			 double vmin,
			 double vmax,
			 int i)
{
  double dx = (vmax - vmin)/(double)(c->nbins);

  return vmin + dx/2.0 + (double)i*dx;
}

static double percent(int proposed, int accepted)
{
  if (proposed > 0) {
    return (double)accepted/(double)proposed;
  }
  return 0.0;
}
   
