###################################################################################################

.PHONY:	r d p sh cr cd cp csh lr ld lp lsh config all install install-headers install-lib\
        install-bin clean distclean
all:	lr lsh

## Load Previous Configuration ####################################################################

-include config.mk

## Configurable options ###########################################################################

# Directory to store object files, libraries, executables, and dependencies:
BUILD_DIR      ?= build

# Include debug-symbols in release builds
MAPLESAT_RELSYM ?= -g

# Sets of compile flags for different build types
MAPLESAT_REL    ?= -O3 -D NDEBUG
MAPLESAT_DEB    ?= -O0 -D DEBUG 
MAPLESAT_PRF    ?= -O3 -D NDEBUG
MAPLESAT_FPIC   ?= -fpic

# GNU Standard Install Prefix
prefix         ?= /usr/local

## Write Configuration  ###########################################################################

config:
	@( echo 'BUILD_DIR?=$(BUILD_DIR)'           ; \
	   echo 'MAPLESAT_RELSYM?=$(MAPLESAT_RELSYM)' ; \
	   echo 'MAPLESAT_REL?=$(MAPLESAT_REL)'       ; \
	   echo 'MAPLESAT_DEB?=$(MAPLESAT_DEB)'       ; \
	   echo 'MAPLESAT_PRF?=$(MAPLESAT_PRF)'       ; \
	   echo 'MAPLESAT_FPIC?=$(MAPLESAT_FPIC)'     ; \
	   echo 'prefix?=$(prefix)'                 ) > config.mk

## Configurable options end #######################################################################

INSTALL ?= install

# GNU Standard Install Variables
exec_prefix ?= $(prefix)
includedir  ?= $(prefix)/include
bindir      ?= $(exec_prefix)/bin
libdir      ?= $(exec_prefix)/lib
datarootdir ?= $(prefix)/share
mandir      ?= $(datarootdir)/man

# Target file names
MAPLESAT      = maplesat#       Name of MiniSat main executable.
MAPLESAT_CORE = maplesat_core#  Name of simplified MiniSat executable (only core solver support).
MAPLESAT_SLIB = lib$(MAPLESAT).a#  Name of MiniSat static library.
MAPLESAT_DLIB = lib$(MAPLESAT).so# Name of MiniSat shared library.

# Shared Library Version
SOMAJOR=2
SOMINOR=1
SORELEASE?=.0#   Declare empty to leave out from library file name.

MAPLESAT_CXXFLAGS = -I. -D __STDC_LIMIT_MACROS -D __STDC_FORMAT_MACROS -Wall -Wno-parentheses -Wextra
MAPLESAT_LDFLAGS  = -Wall -lz

ECHO=@echo
ifeq ($(VERB),)
VERB=@
else
VERB=
endif

SRCS = $(wildcard maplesat/core/*.cc) $(wildcard maplesat/simp/*.cc) $(wildcard maplesat/utils/*.cc)
HDRS = $(wildcard maplesat/mtl/*.h) $(wildcard maplesat/core/*.h) $(wildcard maplesat/simp/*.h) $(wildcard maplesat/utils/*.h)
OBJS = $(filter-out %Main.o, $(SRCS:.cc=.o))

r:	$(BUILD_DIR)/release/bin/$(MAPLESAT)
d:	$(BUILD_DIR)/debug/bin/$(MAPLESAT)
p:	$(BUILD_DIR)/profile/bin/$(MAPLESAT)
sh:	$(BUILD_DIR)/dynamic/bin/$(MAPLESAT)

cr:	$(BUILD_DIR)/release/bin/$(MAPLESAT_CORE)
cd:	$(BUILD_DIR)/debug/bin/$(MAPLESAT_CORE)
cp:	$(BUILD_DIR)/profile/bin/$(MAPLESAT_CORE)
csh:	$(BUILD_DIR)/dynamic/bin/$(MAPLESAT_CORE)

lr:	$(BUILD_DIR)/release/lib/$(MAPLESAT_SLIB)
ld:	$(BUILD_DIR)/debug/lib/$(MAPLESAT_SLIB)
lp:	$(BUILD_DIR)/profile/lib/$(MAPLESAT_SLIB)
lsh:	$(BUILD_DIR)/dynamic/lib/$(MAPLESAT_DLIB).$(SOMAJOR).$(SOMINOR)$(SORELEASE)

## Build-type Compile-flags:
$(BUILD_DIR)/release/%.o:			MAPLESAT_CXXFLAGS +=$(MAPLESAT_REL) $(MAPLESAT_RELSYM)
$(BUILD_DIR)/debug/%.o:				MAPLESAT_CXXFLAGS +=$(MAPLESAT_DEB) -g
$(BUILD_DIR)/profile/%.o:			MAPLESAT_CXXFLAGS +=$(MAPLESAT_PRF) -pg
$(BUILD_DIR)/dynamic/%.o:			MAPLESAT_CXXFLAGS +=$(MAPLESAT_REL) $(MAPLESAT_FPIC)

## Build-type Link-flags:
$(BUILD_DIR)/profile/bin/$(MAPLESAT):		MAPLESAT_LDFLAGS += -pg
$(BUILD_DIR)/release/bin/$(MAPLESAT):		MAPLESAT_LDFLAGS += --static $(MAPLESAT_RELSYM)
$(BUILD_DIR)/profile/bin/$(MAPLESAT_CORE):	MAPLESAT_LDFLAGS += -pg
$(BUILD_DIR)/release/bin/$(MAPLESAT_CORE):	MAPLESAT_LDFLAGS += --static $(MAPLESAT_RELSYM)

## Executable dependencies
$(BUILD_DIR)/release/bin/$(MAPLESAT):	 	$(BUILD_DIR)/release/maplesat/simp/Main.o $(BUILD_DIR)/release/lib/$(MAPLESAT_SLIB)
$(BUILD_DIR)/debug/bin/$(MAPLESAT):	 	$(BUILD_DIR)/debug/maplesat/simp/Main.o $(BUILD_DIR)/debug/lib/$(MAPLESAT_SLIB)
$(BUILD_DIR)/profile/bin/$(MAPLESAT):	 	$(BUILD_DIR)/profile/maplesat/simp/Main.o $(BUILD_DIR)/profile/lib/$(MAPLESAT_SLIB)
# need the main-file be compiled with fpic?
$(BUILD_DIR)/dynamic/bin/$(MAPLESAT):	 	$(BUILD_DIR)/dynamic/maplesat/simp/Main.o $(BUILD_DIR)/dynamic/lib/$(MAPLESAT_DLIB)

## Executable dependencies (core-version)
$(BUILD_DIR)/release/bin/$(MAPLESAT_CORE):	$(BUILD_DIR)/release/maplesat/core/Main.o $(BUILD_DIR)/release/lib/$(MAPLESAT_SLIB)
$(BUILD_DIR)/debug/bin/$(MAPLESAT_CORE):	 	$(BUILD_DIR)/debug/maplesat/core/Main.o $(BUILD_DIR)/debug/lib/$(MAPLESAT_SLIB)
$(BUILD_DIR)/profile/bin/$(MAPLESAT_CORE):	$(BUILD_DIR)/profile/maplesat/core/Main.o $(BUILD_DIR)/profile/lib/$(MAPLESAT_SLIB)
# need the main-file be compiled with fpic?
$(BUILD_DIR)/dynamic/bin/$(MAPLESAT_CORE): 	$(BUILD_DIR)/dynamic/maplesat/core/Main.o $(BUILD_DIR)/dynamic/lib/$(MAPLESAT_DLIB)

## Library dependencies
$(BUILD_DIR)/release/lib/$(MAPLESAT_SLIB):	$(foreach o,$(OBJS),$(BUILD_DIR)/release/$(o))
$(BUILD_DIR)/debug/lib/$(MAPLESAT_SLIB):		$(foreach o,$(OBJS),$(BUILD_DIR)/debug/$(o))
$(BUILD_DIR)/profile/lib/$(MAPLESAT_SLIB):	$(foreach o,$(OBJS),$(BUILD_DIR)/profile/$(o))
$(BUILD_DIR)/dynamic/lib/$(MAPLESAT_DLIB).$(SOMAJOR).$(SOMINOR)$(SORELEASE)\
 $(BUILD_DIR)/dynamic/lib/$(MAPLESAT_DLIB).$(SOMAJOR)\
 $(BUILD_DIR)/dynamic/lib/$(MAPLESAT_DLIB):	$(foreach o,$(OBJS),$(BUILD_DIR)/dynamic/$(o))

## Compile rules (these should be unified, buit I have not yet found a way which works in GNU Make)
$(BUILD_DIR)/release/%.o:	%.cc
	$(ECHO) Compiling: $@
	$(VERB) mkdir -p $(dir $@)
	$(VERB) $(CXX) $(MAPLESAT_CXXFLAGS) $(CXXFLAGS) -c -o $@ $< -MMD -MF $(BUILD_DIR)/release/$*.d

$(BUILD_DIR)/profile/%.o:	%.cc
	$(ECHO) Compiling: $@
	$(VERB) mkdir -p $(dir $@)
	$(VERB) $(CXX) $(MAPLESAT_CXXFLAGS) $(CXXFLAGS) -c -o $@ $< -MMD -MF $(BUILD_DIR)/profile/$*.d

$(BUILD_DIR)/debug/%.o:	%.cc
	$(ECHO) Compiling: $@
	$(VERB) mkdir -p $(dir $@)
	$(VERB) $(CXX) $(MAPLESAT_CXXFLAGS) $(CXXFLAGS) -c -o $@ $< -MMD -MF $(BUILD_DIR)/debug/$*.d

$(BUILD_DIR)/dynamic/%.o:	%.cc
	$(ECHO) Compiling: $@
	$(VERB) mkdir -p $(dir $@)
	$(VERB) $(CXX) $(MAPLESAT_CXXFLAGS) $(CXXFLAGS) -c -o $@ $< -MMD -MF $(BUILD_DIR)/dynamic/$*.d

## Linking rule
$(BUILD_DIR)/release/bin/$(MAPLESAT) $(BUILD_DIR)/debug/bin/$(MAPLESAT) $(BUILD_DIR)/profile/bin/$(MAPLESAT) $(BUILD_DIR)/dynamic/bin/$(MAPLESAT)\
$(BUILD_DIR)/release/bin/$(MAPLESAT_CORE) $(BUILD_DIR)/debug/bin/$(MAPLESAT_CORE) $(BUILD_DIR)/profile/bin/$(MAPLESAT_CORE) $(BUILD_DIR)/dynamic/bin/$(MAPLESAT_CORE):
	$(ECHO) Linking Binary: $@
	$(VERB) mkdir -p $(dir $@)
	$(VERB) $(CXX) $^ $(MAPLESAT_LDFLAGS) $(LDFLAGS) -o $@

## Static Library rule
%/lib/$(MAPLESAT_SLIB):
	$(ECHO) Linking Static Library: $@
	$(VERB) mkdir -p $(dir $@)
	$(VERB) $(AR) -rcs $@ $^

## Shared Library rule
$(BUILD_DIR)/dynamic/lib/$(MAPLESAT_DLIB).$(SOMAJOR).$(SOMINOR)$(SORELEASE)\
 $(BUILD_DIR)/dynamic/lib/$(MAPLESAT_DLIB).$(SOMAJOR)\
 $(BUILD_DIR)/dynamic/lib/$(MAPLESAT_DLIB):
	$(ECHO) Linking Shared Library: $@
	$(VERB) mkdir -p $(dir $@)
	$(VERB) $(CXX) $(MAPLESAT_LDFLAGS) $(LDFLAGS) -o $@ -shared -Wl,-soname,$(MAPLESAT_DLIB).$(SOMAJOR) $^
	$(VERB) ln -sf $(MAPLESAT_DLIB).$(SOMAJOR).$(SOMINOR)$(SORELEASE) $(BUILD_DIR)/dynamic/lib/$(MAPLESAT_DLIB).$(SOMAJOR)
	$(VERB) ln -sf $(MAPLESAT_DLIB).$(SOMAJOR) $(BUILD_DIR)/dynamic/lib/$(MAPLESAT_DLIB)

install:	install-headers install-lib install-bin
install-debug:	install-headers install-lib-debug

install-headers:
#       Create directories
	$(INSTALL) -d $(DESTDIR)$(includedir)/maplesat
	for dir in mtl utils core simp; do \
	  $(INSTALL) -d $(DESTDIR)$(includedir)/maplesat/$$dir ; \
	done
#       Install headers
	for h in $(HDRS) ; do \
	  $(INSTALL) -m 644 $$h $(DESTDIR)$(includedir)/$$h ; \
	done

install-lib-debug: $(BUILD_DIR)/debug/lib/$(MAPLESAT_SLIB)
	$(INSTALL) -d $(DESTDIR)$(libdir)
	$(INSTALL) -m 644 $(BUILD_DIR)/debug/lib/$(MAPLESAT_SLIB) $(DESTDIR)$(libdir)

install-lib: $(BUILD_DIR)/release/lib/$(MAPLESAT_SLIB) $(BUILD_DIR)/dynamic/lib/$(MAPLESAT_DLIB).$(SOMAJOR).$(SOMINOR)$(SORELEASE)
	$(INSTALL) -d $(DESTDIR)$(libdir)
	$(INSTALL) -m 644 $(BUILD_DIR)/dynamic/lib/$(MAPLESAT_DLIB).$(SOMAJOR).$(SOMINOR)$(SORELEASE) $(DESTDIR)$(libdir)
	ln -sf $(MAPLESAT_DLIB).$(SOMAJOR).$(SOMINOR)$(SORELEASE) $(DESTDIR)$(libdir)/$(MAPLESAT_DLIB).$(SOMAJOR)
	ln -sf $(MAPLESAT_DLIB).$(SOMAJOR) $(DESTDIR)$(libdir)/$(MAPLESAT_DLIB)
	$(INSTALL) -m 644 $(BUILD_DIR)/release/lib/$(MAPLESAT_SLIB) $(DESTDIR)$(libdir)

install-bin: $(BUILD_DIR)/dynamic/bin/$(MAPLESAT)
	$(INSTALL) -d $(DESTDIR)$(bindir)
	$(INSTALL) -m 755 $(BUILD_DIR)/dynamic/bin/$(MAPLESAT) $(DESTDIR)$(bindir)

clean:
	rm -rf $(foreach t, release debug profile dynamic, $(foreach o, $(SRCS:.cc=.o), $(BUILD_DIR)/$t/$o)) \
          $(foreach t, release debug profile dynamic, $(foreach d, $(SRCS:.cc=.d), $(BUILD_DIR)/$t/$d)) \
	  $(foreach t, release debug profile dynamic, $(BUILD_DIR)/$t/bin/$(MAPLESAT_CORE) $(BUILD_DIR)/$t/bin/$(MAPLESAT)) \
	  $(foreach t, release debug profile, $(BUILD_DIR)/$t/lib/$(MAPLESAT_SLIB)) \
	  $(BUILD_DIR)/dynamic/lib/$(MAPLESAT_DLIB).$(SOMAJOR).$(SOMINOR)$(SORELEASE)\
	  $(BUILD_DIR)/dynamic/lib/$(MAPLESAT_DLIB).$(SOMAJOR)\
	  $(BUILD_DIR)/dynamic/lib/$(MAPLESAT_DLIB) \
     $(BUILD_DIR)

distclean:	clean
	rm -f config.mk

## Include generated dependencies
-include $(foreach s, $(SRCS:.cc=.d), $(BUILD_DIR)/release/$s)
-include $(foreach s, $(SRCS:.cc=.d), $(BUILD_DIR)/debug/$s)
-include $(foreach s, $(SRCS:.cc=.d), $(BUILD_DIR)/profile/$s)
-include $(foreach s, $(SRCS:.cc=.d), $(BUILD_DIR)/dynamic/$s)
