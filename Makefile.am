AUTOMAKE_OPTIONS = gnu

SUBDIRS = src/apps/TestFastaPlus \
          src/apps/FilterFasta \
          src/apps/GetRandFasta \
          src/apps/SplitFasta
EXTRA_DIST = bin  src/include doc/doxy/html doc/doxy/latex/refman.pdf doc/usage/FastaPlus.pdf demo Doxyfile
 
LDFLAGS_BOOST = \
    $(BOOST_SYSTEM_LDFLAGS) \
    $(BOOST_PROGRAM_OPTIONS_LDFLAGS)

LIBS_BOOST = \
    $(BOOST_SYSTEM_LIBS) \
    $(BOOST_PROGRAM_OPTIONS_LIBS) \
    $(BOOST_FILESYSTEM_LIBS)

AM_CPPFLAGS = -I$(top_srcdir)/src/include $(BOOST_CPPFLAGS)
AM_LDFLAGS = $(LDFLAGS_BOOST)
LIBS += $(LIBS_BOOST)

