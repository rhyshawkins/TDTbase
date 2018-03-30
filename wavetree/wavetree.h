#ifndef wavetree_h
#define wavetree_h

typedef enum {

  WT_PERTURB_INVALID = -1,
  WT_PERTURB_NONE = 0,
  WT_PERTURB_BIRTH = 1,
  WT_PERTURB_DEATH = 2,
  WT_PERTURB_VALUE = 3,
  WT_PERTURB_MOVE  = 4,
  WT_PERTURB_HIERARCHICAL = 5,
  WT_PERTURB_PTEXCHANGE = 6,
  WT_PERTURB_PTMODELEXCHANGE = 7

} wavetree_perturb_t;

#endif /* wavetree_h */
