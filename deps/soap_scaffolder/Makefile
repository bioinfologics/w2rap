####################################
#
#
###################################

GCCVERSIONMAJOR := $(shell expr `$(CC) -dumpversion | cut -f1 -d.` \>= 4)
GCCVERSIONMINOR := $(shell expr `$(CC) -dumpversion | cut -f2 -d.` \>= 4)
ifdef debug
CFLAGS=         -O0 -g -fomit-frame-pointer #-static #-mcrc32 -march=core2 -msse4.1 -msse4.2
else
CFLAGS=         -O4 -fomit-frame-pointer #-static #-mcrc32 -march=core2 -msse4.1 -msse4.2
endif
DFLAGS=         
PREPARE_OBJS= prepare.o check.o
MAP_OBJS=   map.o newhash.o kmer.o prlHashCtg.o seq.o \
        check.o hashFunction.o prlRead2Ctg.o lib.o attachPEinfo.o \
        connect.o mem_manager.o stack.o readseq1by1.o
SCAFF_OBJS= scaffold.o orderContig.o fibHeap.o fib.o attachPEinfo.o check.o \
        lib.o connect.o mem_manager.o stack.o output_contig.o \
        loadPreGraph.o kmer.o arc.o darray.o loadGraph.o seq.o prlReadFillGap.o localAsm.o newhash.o searchPath.o
INCLUDES=   -Iinc
SUBDIRS=    .
LIBPATH=    -L/lib64 -L/usr/lib64 -L./inc -L./inc -L/tgac/software/testing/libraries/zlib/1.2.8/x86_64/lib/
LIBS=       -pthread -lz -lm 
EXTRA_FLAGS=

BIT_ERR = 0
ifeq (,$(findstring $(shell uname -m), x86_64 ppc64 ia64))
BIT_ERR = 1
endif

CFLAGS += -DMER127

ifneq (,$(findstring Linux,$(shell uname)))
EXTRA_FLAGS += -Wl,--hash-style=both
LIBS += -lrt 
endif

ifneq (,$(findstring Unix,$(shell uname)))
EXTRA_FLAGS += -Wl,--hash-style=both
LIBS += -lrt
endif

ifneq (,$(findstring $(shell uname -m), x86_64))
CFLAGS += -m64
endif

ifneq (,$(findstring $(shell uname -m), ia64))
CFLAGS += 
endif

ifneq (,$(findstring $(shell uname -m), ppc64))
CFLAGS += -mpowerpc64
endif

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< 

all:		clean s_prepare s_map s_scaff

#SOAPdenovo

.PHONY:all clean install

envTest:
		@test $(BIT_ERR) != 1 || sh -c 'echo "Fatal: 64bit CPU and Operating System required!";false;'
		@test $(GCCVERSIONMAJOR) == 1 || sh -c 'echo "GCC version lower than 4.4.0";false;'
		@test $(GCCVERSIONMINOR) == 1 || sh -c 'echo "GCC version lower than 4.4.0";false;'

s_prepare:	$(PREPARE_OBJS)
		$(CC) $(CFLAGS) -o s_prepare $(PREPARE_OBJS) $(LIBPATH) $(LIBS) $(ENTRAFLAGS)

s_map:	$(MAP_OBJS)
		$(CC) $(CFLAGS) -o s_map $(MAP_OBJS) $(LIBPATH) $(LIBS) $(ENTRAFLAGS)

s_scaff:	$(SCAFF_OBJS)
		$(CC) $(CFLAGS) -o s_scaff $(SCAFF_OBJS) $(LIBPATH) $(LIBS) $(ENTRAFLAGS)

clean:
		@rm -fr gmon.out *.o a.out *.exe *.dSYM s_prepare s_scaff s_map *~ *.a *.so.* *.so *.dylib
		@printf "cleaning done.\n";
