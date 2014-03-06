C=		gcc
CXX=		g++
CFLAGS=		-g -Wall -O3 -std=c99 -fopenmp
CXXFLAGS=	$(CFLAGS)
DFLAGS=		-DHAVE_PTHREAD
DBGFLAGS= 	#-DDEBUG_ENABLED
OBJS=		common.o dbs.o tol.o precalc.o query.o step.o occ.o selectmarker.o candidate.o prepDB.o filt.o main.o
PROG=		speedb
INCLUDES=	
LIBS=		-lm -lz -lpthread

.SUFFIXES:.c .o .cc

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $(DBGFLAGS) $< -o $@
.cc.o:
		$(CXX) -c $(CXXFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

speedb:$(OBJS)
		$(CC) $(CFLAGS) $(DFLAGS) $(OBJS) -o $@ $(LIBS)

clean:
		rm -f gmon.out *.o a.out $(PROG) *~ *.a
