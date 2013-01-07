
# Platform-specific settings
ifeq ($(MAKE),mingw32-make)
  # assume Windows
  PLATFORM = WINDOWS
  RM = del
  LIBEXT = dll
else
  # assume Linux
  PLATFORM = LINUX
  RM = rm -f
  LIBEXT = so
endif

all: $(PROGS) $(addsuffix .dir, $(SUBDIRS)) $(addprefix lib, $(addsuffix .$(LIBEXT), $(LIBS)))

# GCC
CC = gcc
FC = gfortran
CFLAGS += -std=c99 -fPIC -g -fbounds-check -Wall -D$(PLATFORM)
FFLAGS += -fPIC -g -fbounds-check -Wall
LDLIBS += -lm
SHARED = -shared -Wl,-soname,$(basename $<)

%: %.c
	$(CC) $(CFLAGS) $(LDLIBS) -o $@ $^

lib%.$(LIBEXT): %.c $(addsuffix .o, $(COMMON))
	$(CC) $(CFLAGS) $(LDLIBS) $(SHARED) -o $@ $^

lib%.$(LIBEXT): %.f $(addsuffix .o, $(COMMON))
	$(FC) $(FFLAGS) $(LDLIBS) $(SHARED) -o $@ $^


%.dir:
	@(cd $(basename $@) && $(MAKE))

%.clean:
	@(cd $(basename $@) && $(MAKE) clean)

clean: $(addsuffix .clean, $(SUBDIRS))
	@( $(RM) $(PROGS) *.o *.$(LIBEXT) *.pyc && exit 0 )

