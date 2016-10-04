
AM_CPPFLAGS = \
	-I$(blitzdir)/include \
	-I$(libconfigdir)/include \
	-I$(top_srcdir)/src

AM_LDFLAGS = -static

noinst_LTLIBRARIES = libupd.la
libupd_la_SOURCES = update.cpp update.h \
	update_full.cpp update_full.h \
	update_no_dispersion.cpp update_no_dispersion.h
