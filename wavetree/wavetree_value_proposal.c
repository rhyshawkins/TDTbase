
#include <stdio.h>
#include <stdlib.h>

#include "wavetree_value_proposal.h"

void
wavetree_value_destroy(wavetree_value_t *p)
{
  p->destroy(p->user);
  free(p);
}

