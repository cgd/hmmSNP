AM_CC = gcc
AM_CFLAGS =  -std=gnu99  -Wall -pedantic -Werror
AM_LIBS = -lm 
AM_CPPFLAGS = -I$(top_srcdir)/src/include

SUBDIRS = include

bin_PROGRAMS = hmmSNP
hmmSNP_SOURCES = main.c dataIO.c CLI.c hmm_em.c hmmMatrix.c dataFilter.c \
	filter.c path.c fill.c graph.c prune.c hmmUtil.c rand.c sort_max_trace.c \
	csv.c temp_file_cleanup.c version.c state_sorting.c


#force version.c to be compiled every time, this way __DATE__ gets updated
version.o: FORCE
	$(COMPILE) -c $(top_srcdir)/src/version.c

FORCE:
