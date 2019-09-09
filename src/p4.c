/* Copyright (c) 2007,2019 Simon Urbanek

   Part of proj4 R package, license: GPL v2 */

#include <R.h>
#ifdef USE_PROJ6_API
#include <proj.h>
#else
#include <proj_api.h>
#endif

#define F_INVERSE 1
#define F_DEG     2

/* .C is more efficient here, because we can re-use the vectors */

void project_(char **proj, int *n, double *x, double *y, int *f){
  int i = 0, j = *n, deg = (*f & F_DEG), inv = (*f & F_INVERSE);
#ifdef USE_PROJ6_API
  PJ *pj;
  PJ_COORD c;

  if (!(pj = proj_create(0, *proj)))
    error(proj_errno_string(proj_errno(0)));

  while (i < j){
    if (ISNA(x[i]) || ISNA(y[i])){ /* if either is NA, set both to NA */
      x[i] = R_NaReal;
      y[i] = R_NaReal;
    } else {
      if (!inv && deg) {
	c.uv.u = proj_torad(x[i]);
	c.uv.v = proj_torad(y[i]);
      } else {
	c.uv.u = x[i];
	c.uv.v = y[i];
      }

      c = proj_trans(pj, inv ? PJ_INV : PJ_FWD, c);

      if (c.uv.u == HUGE_VAL) {
	int r = proj_errno(pj);
	proj_destroy(pj);
	error(proj_errno_string(r));
      }
      if (inv && deg) {
	x[i] = proj_todeg(c.uv.u);
	y[i] = proj_todeg(c.uv.v);
      } else {
	x[i] = c.uv.u;
	y[i] = c.uv.v;
      }
    }
    i++;
  }

  proj_destroy(pj);
#else
  projUV c;
  projPJ pj;
  
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

      c = inv ? pj_inv(c, pj) : pj_fwd(c, pj);

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
#endif
}

void transform_(char **psrc, char **pdst, int *n, double *x, double *y, double *z) {
#ifdef USE_PROJ6_API
  PJ *pj, *pj2;
  int i, N = *n; /* if we want to support large vectors we'd have to change type of n first */
  int r;
  int in_rad2deg = 0, out_rad2deg = 0;

  if (!(pj = proj_create_crs_to_crs(0, *psrc, *pdst, 0)))
    error(proj_errno_string(proj_errno(0)));

  if (!(pj2 = proj_normalize_for_visualization(0, pj))) {
    int r = proj_errno(pj);
    proj_destroy(pj);
    error(proj_errno_string(r));
  }

  proj_destroy(pj);
  pj = pj2;

  /* this is the nasty bit - the new API auto-coverts deg to rad, so
     in order to be compatible, we have to identify cases when rad
     was previously used and convert it to deg so that the auto-convert
     can convert it back. Yes, utterly stupid, and requires nasty heuristic,
     but there is not much we can do .. */
  if (proj_angular_input(pj, PJ_FWD) == 0 && (pj2 = proj_get_source_crs(0, pj))) {
    const char *pstr = proj_as_proj_string(0, pj2, PJ_PROJ_4, 0);
    if (strstr(pstr, "=longlat")) in_rad2deg = 1;
    proj_destroy(pj2);
  }

  if (proj_angular_output(pj, PJ_FWD) == 0 && (pj2 = proj_get_target_crs(0, pj))) {
    const char *pstr = proj_as_proj_string(0, pj2, PJ_PROJ_4, 0);
    if (strstr(pstr, "=longlat")) out_rad2deg = 1;
    proj_destroy(pj2);
  }

  if (in_rad2deg) for (i = 0; i < N; i++) {
      x[i] = proj_todeg(x[i]);
      y[i] = proj_todeg(y[i]);
    }

  proj_trans_generic(pj, PJ_FWD,
		     x, sizeof(*x), N,
		     y, sizeof(*y), N,
		     z, sizeof(*z), N,
		     0, 0, 0);

  if (out_rad2deg) for (i = 0; i < N; i++) {
      x[i] = proj_torad(x[i]);
      y[i] = proj_torad(y[i]);
    }

  r = proj_errno(pj);

  proj_destroy(pj);

  if (r)
    error(proj_errno_string(r));

#else
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
#endif
}
