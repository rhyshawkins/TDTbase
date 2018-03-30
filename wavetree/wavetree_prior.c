
#include <stdio.h>
#include <stdlib.h>

#include "wavetree_prior.h"

void
wavetree_prior_destroy(wavetree_prior_t *p)
{
  p->destroy(p->user);

  free(p);
}



