#ifndef chain_history_h
#define chain_history_h

#include "wavetree.h"
#include "multiset_int_double.h"

typedef enum {
  CH_NOCHANGE = WT_PERTURB_NONE,
  CH_BIRTH = WT_PERTURB_BIRTH,
  CH_DEATH = WT_PERTURB_DEATH,
  CH_VALUE = WT_PERTURB_VALUE,
  CH_MOVE = WT_PERTURB_MOVE, 
  CH_HIERARCHICAL = WT_PERTURB_HIERARCHICAL,
  CH_PTEXCHANGE = WT_PERTURB_PTEXCHANGE,
  CH_PTMODELEXCHANGE = WT_PERTURB_PTMODELEXCHANGE,
  CH_HYPER,
  CH_INITIALISE
} chain_history_step_t;

/*
 * This structure records most steps (except for initialisation where we need the
 * whole structure. Example values for each type of step.
 *
 */
typedef struct {
  struct {
    chain_history_step_t type;
    int accepted;
    double likelihood;
    double temperature;
    double hierarchical;
  } header;

  union {
    struct {
      int node_depth;
      int node_id;
      double new_value;
    } birth;
    struct {
      int node_depth;
      int node_id;
      double old_value;
    } death;
    struct {
      int node_depth;
      int node_id;
      double new_value;
      double old_value;
    } value;
    struct {
      int node_depth;
      int node_id;
      int new_node_id;
      double new_value;
      double old_value;
    } move;
    struct {
      double old_value;
      double new_value;
    } hierarchical;
    struct {
      double old_temperature;
    } ptexchange;
    struct {
      int index;
      double old_value;
      double new_value;
    } hyper;
  } perturbation;

} chain_history_change_t;

typedef struct _chain_history chain_history_t;

chain_history_t *
chain_history_create(int maxsteps);

void
chain_history_destroy(chain_history_t *ch);

/*
 * Set the starting model from a set (usually obtained from wavetomo2d_t or wavetomo3d_t).
 */
int
chain_history_initialise(chain_history_t *ch,
			 const multiset_int_double_t *S_v,
			 double likelihood,
			 double temperature,
			 double hierarhical);
/*
 * Set the current model to the initial. Clear all deltas
 */
int
chain_history_reset(chain_history_t *ch);

/*
 * Add a step to the current history
 */
int
chain_history_add_step(chain_history_t *ch,
		       const chain_history_change_t *step);


typedef size_t (*ch_read_t)(void *ptr, size_t size, size_t nmemb, void *fp);
typedef size_t (*ch_write_t)(const void *ptr, size_t size, size_t nmemb, void *fp);

int
chain_history_write(chain_history_t *ch,
		    ch_write_t write_function,
		    void *fp);

int
chain_history_read(chain_history_t *ch,
		   ch_read_t read_function,
		   void *fp);

typedef int (*chain_history_replay_function_t)(int i,
					       void *user,
					       const chain_history_change_t *step,
					       const multiset_int_double_t *S_v);

int
chain_history_replay(chain_history_t *ch,
		     multiset_int_double_t *S_v,
		     chain_history_replay_function_t cb,
		     void *user);

int
chain_history_nsteps(chain_history_t *ch);

int
chain_history_full(chain_history_t *ch);

#endif /* chain_history_h */
