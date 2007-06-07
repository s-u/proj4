/* Copyright (c) 2007 Simon Urbanek

   Part of proj4 R package, license: GPL v2 */

#include <R.h>
#include <proj_api.h>

#define F_INVERSE 1
#define F_DEG     2

/* .C is more efficient here, because we can re-use the vectors */

void project(char **proj, int *n, double *x, double *y, int *f){
  projUV c;
  projPJ pj;
  int i = 0, j = *n, deg = (*f & F_DEG), inv = (*f & F_INVERSE);
  
  if (!(pj = pj_init_plus(*proj))) 
    error(pj_strerrno(*pj_get_errno_ref()));

  while (i < j){
    if (ISNA(x[i]) || ISNA(y[i])){ /* if either is NA, set both to NA */
      x[i] = R_NaReal;
      y[i] = R_NaReal;
    } else {
      if (!inv && deg) {
	c.u = x[i] * DEG_TO_RAD;
	c.v = y[i] * DEG_TO_RAD;
      } else {
	c.u = x[i];
	c.v = y[i];
      }
      c = inv?pj_inv(c, pj):pj_fwd(c, pj);
      if (c.u == HUGE_VAL) {
	pj_free(pj);
	error(pj_strerrno(*pj_get_errno_ref()));
      }
      if (inv && deg) {
	x[i] = c.u * RAD_TO_DEG;
	y[i] = c.v * RAD_TO_DEG;
      } else {
	x[i] = c.u;
	y[i] = c.v;
      }
    }
    i++;
  }

  pj_free(pj);
}

void transform(char **psrc, char **pdst, int *n, double *x, double *y, double *z) {
  projPJ ps, pd;
  int r;

  if (!(ps = pj_init_plus(*psrc)))
    error(pj_strerrno(*pj_get_errno_ref()));
  if (!(pd = pj_init_plus(*pdst))) {
    pj_free(ps);
    error(pj_strerrno(*pj_get_errno_ref()));
  }
  
  r = pj_transform(ps, pd, *n, 0, x, y, z);

  pj_free(ps);
  pj_free(pd);

  if (r)
    error(pj_strerrno(*pj_get_errno_ref()));    
}
