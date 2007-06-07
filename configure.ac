AC_INIT(proj4, 1.0, Simon.Urbanek@r-project.org)
AC_CONFIG_SRCDIR(src/p4.c)

# find R home and set correct compiler + flags
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  AC_MSG_ERROR([cannot determine R_HOME. Make sure you use R CMD INSTALL!])
fi
RBIN="${R_HOME}/bin/R"

: ${CC=`"${RBIN}" CMD config CC`}
: ${CPP=`"${RBIN}" CMD config CPP`}
: ${CFLAGS=`"${RBIN}" CMD config CFLAGS`}

# honor PKG_xx overrides
LIBS="${LIBS} ${PKG_LIBS}"
CPPFLAGS="${CPPFLAGS} ${PKG_CPPFLAGS}"

proj4ok=yes
AC_CHECK_HEADERS(proj_api.h,,proj4ok=no)
AC_CHECK_LIB(proj,pj_init_plus,,proj4ok=no)
if test "${proj4ok}" = no; then
   AC_MSG_ERROR([libproj and/or proj_api.h not found in standard search locations.

*** Install PROJ.4 and if necessary set PKG_CPPFLAGS/PKG_LIBS accordingly.])
fi

AC_ARG_VAR([PKG_CPPFLAGS],[custom C preprocessor flags for packages compilation])
AC_ARG_VAR([PKG_LIBS],[custom libraries for package compilation])

AC_CONFIG_FILES(src/Makevars)
AC_OUTPUT


