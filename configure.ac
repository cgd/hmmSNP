#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.


AC_PREREQ(2.59)
AC_INIT(hmmSNP, 1.3.2, glen.beane@jax.org)
AC_CONFIG_SRCDIR([src/CLI.c])
AC_CONFIG_AUX_DIR([buildutils])
AM_CONFIG_HEADER([src/include/config.h])
AM_INIT_AUTOMAKE(foreign)

# Checks for programs.
AC_PROG_CC

# Checks for libraries.
AC_CHECK_LIB([m], [sqrt])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([limits.h stdlib.h string.h sys/param.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STAT
AC_FUNC_STRTOD
AC_CHECK_FUNCS([memmove sqrt strcspn strdup strerror strtol])

AC_CONFIG_FILES([Makefile
                 src/Makefile
		 src/include/Makefile])
AC_OUTPUT
