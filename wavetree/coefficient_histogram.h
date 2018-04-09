//
//    Wavetree Library : A library for performed trans-dimensional tree inversion,
//    See
//
//      R Hawkins and M Sambridge, "Geophysical imaging using trans-dimensional trees",
//      Geophysical Journal International, 2015, 203:2, 972 - 1000,
//      https://doi.org/10.1093/gji/ggv326
//    
//    Copyright (C) 2014 - 2018 Rhys Hawkins
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

#ifndef coefficient_histogram_h
#define coefficient_histogram_h

typedef int (*ch_coord_to_index_t)(void *user, int i, int j, int k, int depth);
typedef int (*ch_index_to_coord_t)(void *user, int index, int *i, int *j, int *k, int *depth);

struct coefficient_histogram {

  ch_coord_to_index_t coordtoindex;
  ch_index_to_coord_t indextocoord;
  void *ch_user;
  
  int ncoeff;
  int nbins;

  double gvmin;
  double gvmax;

  double *vmin; /* [ncoeff] */
  double *vmax; /* [ncoeff] */

  int **counts; /* [ncoeff][nbins] */
  int *under;   /* [ncoeff] */
  int *over;    /* [ncoeff] */

  double *rmin; /* [ncoeff] */
  double *rmax; /* [ncoeff] */
  double *rmean;/* [ncoeff] */
  double *rstd; /* [ncoeff] */
  int *n;       /* [ncoeff] */

  double *valpha;      /* [ncoeff] */
  double *valpha_mean; /* [ncoeff] */
  int *valpha_n;       /* [ncoeff] */
  
  int *pb;      /* [ncoeff] */
  int *ab;      /* [ncoeff] */
  int *pd;      /* [ncoeff] */
  int *ad;      /* [ncoeff] */
  int *pv;      /* [ncoeff] */
  int *av;      /* [ncoeff] */

};
typedef struct coefficient_histogram coefficient_histogram_t;

coefficient_histogram_t *
coefficient_histogram_create(int ncoeff,
			     int nbins,
			     double vmin,
			     double vmax,
			     ch_coord_to_index_t coordtoindex,
			     ch_index_to_coord_t indextocoord,
			     void *ch_user);

void
coefficient_histogram_destroy(coefficient_histogram_t *c);

int
coefficient_histogram_save(coefficient_histogram_t *c,
			   const char *filename);

int
coefficient_histogram_load(coefficient_histogram_t *c,
			   const char *filename);


int
coefficient_histogram_coord_to_index(const coefficient_histogram_t *c, int i, int j, int k, int depth);

int
coefficient_histogram_index_to_coord(const coefficient_histogram_t *c,
				     int index, int *i, int *j, int *k, int *depth);

int
coefficient_histogram_reset(coefficient_histogram_t *c);

int 
coefficient_histogram_set_range(coefficient_histogram_t *c, 
				int index,
				double vmin,
				double vmax);

int
coefficient_histogram_sample(coefficient_histogram_t *c, int index, double value);

/*
 * For recording per coefficient birth information
 */
int 
coefficient_histogram_propose_birth(coefficient_histogram_t *c, int index);

int 
coefficient_histogram_accept_birth(coefficient_histogram_t *c, int index, double value);
			     
int 
coefficient_histogram_reject_birth(coefficient_histogram_t *c, int index, double value);

/*
 * For recording per coefficient death information
 */
int 
coefficient_histogram_propose_death(coefficient_histogram_t *c, int index);

int 
coefficient_histogram_accept_death(coefficient_histogram_t *c, int index);
			     
/*
 * For recording per coefficient value information
 */
int 
coefficient_histogram_propose_value(coefficient_histogram_t *c, int index);

int 
coefficient_histogram_accept_value(coefficient_histogram_t *c, int index, double value);
			     
int 
coefficient_histogram_reject_value(coefficient_histogram_t *c, int index, double value);

int
coefficient_histogram_sample_value_alpha(coefficient_histogram_t *c, int index, double alpha);

int 
coefficient_histogram_get_coefficient_mean_std(coefficient_histogram_t *c,
					       int index,
					       double *mean,
					       double *std);

int
coefficient_histogram_get_accept_reject(coefficient_histogram_t *c,
					int index,
					int *propose,
					int *accept);

int
coefficient_histogram_finalise(coefficient_histogram_t *c);

int
coefficient_histogram_save_2D(const coefficient_histogram_t *c,
			      const char *filename);

int
coefficient_histogram_save_min(const coefficient_histogram_t *c,
			       const char *filename);

int
coefficient_histogram_save_max(const coefficient_histogram_t *c,
			       const char *filename);

int
coefficient_histogram_save_mean(const coefficient_histogram_t *c,
				const char *filename);

int
coefficient_histogram_save_std(const coefficient_histogram_t *c,
			       const char *filename);


typedef int (*coefficient_histogram_index_filter_t)(void *user, int subset, int index);
int
coefficient_histogram_save_aggregated_histogram(const coefficient_histogram_t *c,
						const char *filename,
						coefficient_histogram_index_filter_t filter,
						void *user,
						int subset);

int
coefficient_histogram_save_aggregated_histogram_image(const coefficient_histogram_t *c,
						      const char *filename,
						      coefficient_histogram_index_filter_t filter,
						      void *user,
						      int subset);

int
coefficient_histogram_save_acceptance(const coefficient_histogram_t *c,
				      const char *filename);


#endif /* coefficient_histogram */
