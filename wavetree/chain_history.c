
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "chain_history.h"

#include "slog.h"

struct _chain_history {

  multiset_int_double_t *S_v_initial;
  multiset_int_double_t *S_v_current;

  int maxsteps;
  int nsteps;
  chain_history_change_t *steps;
};

static int do_step(multiset_int_double_t *S_v,
		   chain_history_change_t *step);

chain_history_t *
chain_history_create(int maxsteps)
{
  chain_history_t *ch;

  ch = malloc(sizeof(chain_history_t));
  if (ch == NULL) {
    ERROR("failed to allocate structure");
    return NULL;
  }

  ch->S_v_initial = multiset_int_double_create();
  ch->S_v_current = multiset_int_double_create();
  if (ch->S_v_initial == NULL ||
      ch->S_v_current == NULL) {
    ERROR("failed to allocate multisets");
    return NULL;
  }

  ch->maxsteps = maxsteps;
  ch->nsteps = 0;
  ch->steps = malloc(sizeof(chain_history_change_t) * maxsteps);
  if (ch->steps == NULL) {
    ERROR("failed to allocate steps");
    return NULL;
  }

  return ch;
}

void
chain_history_destroy(chain_history_t *ch)
{
  if (ch != NULL) {
    multiset_int_double_destroy(ch->S_v_initial);
    multiset_int_double_destroy(ch->S_v_current);
    free(ch->steps);
    free(ch);
  }
}

int
chain_history_initialise(chain_history_t *ch,
			 const multiset_int_double_t *S_v,
			 double likelihood,
			 double temperature,
			 double hierarchical)
{
  if (multiset_int_double_clone(ch->S_v_initial, S_v) < 0) {
    ERROR("failed to clone S_v");
    return -1;
  }
  
  if (multiset_int_double_clone(ch->S_v_current, S_v) < 0) {
    ERROR("failed to clone S_v");
    return -1;
  }

  memset(&(ch->steps[0]), 0, sizeof(chain_history_change_t));
  ch->steps[0].header.type = CH_INITIALISE;
  ch->steps[0].header.likelihood = likelihood;
  ch->steps[0].header.temperature = temperature;
  ch->steps[0].header.hierarchical = hierarchical;

  ch->nsteps = 1;
  
  return 0;
}

int
chain_history_reset(chain_history_t *ch)
{
  if (ch->nsteps <= 0) {
    ERROR("uninitialised");
    return -1;
  }
  
  if (ch->nsteps == 1) {
    return 0;
  }

  if (multiset_int_double_clone(ch->S_v_initial, ch->S_v_current) < 0) {
    ERROR("failed to clone");
    return -1;
  }

  memset(&(ch->steps[0]), 0, sizeof(chain_history_change_t));
  ch->steps[0].header.type = CH_INITIALISE;
  ch->steps[0].header.likelihood = ch->steps[ch->nsteps - 1].header.likelihood;
  ch->steps[0].header.temperature = ch->steps[ch->nsteps - 1].header.temperature;
  ch->steps[0].header.hierarchical = ch->steps[ch->nsteps - 1].header.hierarchical;
  
  ch->nsteps = 1;

  return 0;
}

int
chain_history_add_step(chain_history_t *ch,
		       const chain_history_change_t *step)
{
  if (ch->nsteps == ch->maxsteps) {
    ERROR("full");
    return -1;
  }

  memcpy(&(ch->steps[ch->nsteps]), step, sizeof(chain_history_change_t));

  if (do_step(ch->S_v_current, &(ch->steps[ch->nsteps])) < 0) {
    ERROR("failed to update model");
    return -1;
  }
  
  ch->nsteps ++;

  return 0;
}

int
chain_history_write(chain_history_t *ch,
		    ch_write_t write_function,
		    void *fp)
{
  int i;
  
  if (ch->steps[0].header.type != CH_INITIALISE) {
    ERROR("first step is not initialisation");
    return -1;
  }

  if (ch->nsteps <= 0) {
    ERROR("no steps");
    return -1;
  }
  
  if (write_function(&(ch->nsteps), sizeof(int), 1, fp) != 1) {
    ERROR("failed to write nsteps");
    return -1;
  }

  /*
   * Write initial S_v
   */
  if (multiset_int_double_write_binary(ch->S_v_initial, write_function, fp) < 0) {
    ERROR("failed to write initial S_v");
    return -1;
  }

  /*
   * Write initial likelihood
   */
  if (write_function(&(ch->steps[0].header.likelihood), sizeof(double), 1, fp) != 1) {
    ERROR("failed to write likelihood");
    return -1;
  }
  
  /*
   * Write initial temperature
   */
  if (write_function(&(ch->steps[0].header.temperature), sizeof(double), 1, fp) != 1) {
    ERROR("failed to write temperature");
    return -1;
  }

  /*
   * Write initial hierarchical
   */
  if (write_function(&(ch->steps[0].header.hierarchical), sizeof(double), 1, fp) != 1) {
    ERROR("failed to write hierarchical");
    return -1;
  }

  for (i = 1; i < ch->nsteps; i ++) {
    if (write_function(&(ch->steps[i]), sizeof(chain_history_change_t), 1, fp) != 1) {
      ERROR("failed to write step %d", i);
      return -1;
    }
  }

  return 0;
}

int 
chain_history_read(chain_history_t *ch,
		   ch_read_t read_function,
		   void *fp)
{
  int i;
  
  if (read_function(&(ch->nsteps), sizeof(int), 1, fp) != 1) {
    /* ERROR("failed to read nsteps"); */
    return -1;
  }

  if (ch->nsteps > ch->maxsteps) {
    ERROR("invalid number of steps (%d > %d)", ch->nsteps, ch->maxsteps);
    return -1;
  }

  /*
   * Read initial S_v
   */
  if (multiset_int_double_read_binary(ch->S_v_initial, read_function, fp) < 0) {
    ERROR("failed to read initial model");
    return -1;
  }

  /*
   * Set initial step
   */
  memset(&(ch->steps[0]), 0, sizeof(chain_history_change_t));
  ch->steps[0].header.type = CH_INITIALISE;

  /*
   * Read initial likelihood
   */
  if (read_function(&(ch->steps[0].header.likelihood), sizeof(double), 1, fp) != 1) {
    ERROR("failed to read likelihood");
    return -1;
  }
  
  /*
   * Read initial temperature
   */
  if (read_function(&(ch->steps[0].header.temperature), sizeof(double), 1, fp) != 1) {
    ERROR("failed to read hierarchical");
    return -1;
  }

  /*
   * Read initial hierarchical
   */
  if (read_function(&(ch->steps[0].header.hierarchical), sizeof(double), 1, fp) != 1) {
    ERROR("failed to read hierarchical");
    return -1;
  }

  for (i = 1; i < ch->nsteps; i ++) {
    
    if (read_function(ch->steps + i, sizeof(chain_history_change_t), 1, fp) != 1) {
      if (feof(fp)) {
	/*
	 * Truncated file
	 */
	ch->nsteps = i;
	return 0;
      } else {
	ERROR("failed to read step %d", i);
	return -1;
      }
    }

  }

  return 0;
    
}

int
chain_history_replay(chain_history_t *ch,
		     multiset_int_double_t *S_v,
		     chain_history_replay_function_t cb,
		     void *user)
{
  int i;

  if (ch->nsteps < 0) {
    return -1;
  }
  
  if (ch->nsteps == 0) {
    return 0;
  }

  /*
   * Initialise to the first
   */
  if (multiset_int_double_clone(S_v, ch->S_v_initial) < 0) {
    ERROR("failed to clone initial model");
    return -1;
  }

  for (i = 1; i < ch->nsteps; i ++) {

    /* printf("  %10.6f %10.6f\n", ch->steps[i].new_value, ch->steps[i].old_value); */
    
    /*
     * Alter model
     */
    if (do_step(S_v, &(ch->steps[i])) < 0) {
      ERROR("failed to do step");
      return -1;
    }

    if (cb(i, user, &(ch->steps[i]), S_v) < 0) {
      return -1;
    }
  }

  return 0;
}

int
chain_history_nsteps(chain_history_t *ch)
{
  return ch->nsteps;
}

int
chain_history_full(chain_history_t *ch)
{
  return ch->nsteps == ch->maxsteps;
}

static int do_step(multiset_int_double_t *S_v,
		   chain_history_change_t *step)
{
  switch (step->header.type) {
    
  case CH_NOCHANGE:
    /* printf("ch: none\n"); */
    break;
    
  case CH_BIRTH:
    /* printf("ch: birth %d %d %f\n", */
    /* 	   step->perturbation.birth.node_id, */
    /* 	   step->perturbation.birth.node_depth, */
    /* 	   step->perturbation.birth.new_value); */

    if (step->header.accepted) {
      if (multiset_int_double_insert(S_v,
				     step->perturbation.birth.node_id,
				     step->perturbation.birth.node_depth,
				     step->perturbation.birth.new_value) < 0) {
	ERROR("failed to recreate birth");
	return -1;
      }
    }
    break;
    
  case CH_DEATH:
    /* printf("ch: death %d %d\n", */
    /* 	   step->perturbation.death.node_id, */
    /* 	   step->perturbation.death.node_depth); */
    if (step->header.accepted) {
      if (multiset_int_double_remove(S_v,
				     step->perturbation.death.node_id,
				     step->perturbation.death.node_depth) < 0) {
	ERROR("failed to recreate death");
	return -1;
      }
    }
    break;
    
  case CH_VALUE:
    /* printf("ch: value %d %d %f (%d)\n", */
    /* 	   step->perturbation.value.node_id, */
    /* 	   step->perturbation.value.node_depth, */
    /* 	   step->perturbation.value.new_value, */
    /* 	   multiset_int_double_total_count(S_v)); */
    if (step->header.accepted) {
      if (multiset_int_double_set(S_v,
				  step->perturbation.value.node_id,
				  step->perturbation.value.node_depth,
				  step->perturbation.value.new_value) < 0) {
	ERROR("failed to recreate change value (%d %d)",
		step->perturbation.value.node_id, step->perturbation.value.node_depth);
	multiset_int_double_dump(S_v);
	return -1;
      }
    }
    break;
    
  case CH_MOVE:

    if (step->header.accepted) {
      if (multiset_int_double_remove(S_v,
    				     step->perturbation.move.node_id,
    				     step->perturbation.move.node_depth) < 0) {
    	ERROR("failed to recreate move (delete)");
    	return -1;
      }
      
      if (multiset_int_double_insert(S_v,
    				     step->perturbation.move.new_node_id,
    				     step->perturbation.move.node_depth,
    				     step->perturbation.move.new_value) < 0) {
    	ERROR("failed to recreate move (insert)");
    	return -1;
      }
    }
    break;

  case CH_HIERARCHICAL:
    /*
     * Does nothing to model S_v!
     */
    break;
    

  case CH_HYPER:
    /*
     * Do nothing - user controlled
     */
    break;
    
  case CH_INITIALISE:
    ERROR("initialisation unimplemented");
    return -1;
    
  case CH_PTEXCHANGE:
    ERROR("exchange unimplemented");
    return -1;

    
  default:
    ERROR("invalid step type %d", step->header.type);
    return -1;

  }

  return 0;
}

