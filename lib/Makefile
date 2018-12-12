
CC = gcc
#
# For intel compilers
#
#CC = icc

CFLAGS = -c -g -Wall

CFLAGS += -O3

AR = ar
ARFLAGS = -r

OBJS = bbox2d.o \
	delaunay2d.o \
	rjmcmc_debug.o \
	solver2x2.o \
	tracking.o \
	slog.o

TARGETS = libdelaunay.a

all : $(TARGETS)

libdelaunay.a : $(OBJS)
	$(AR) $(ARFLAGS) $@ $(OBJS)

%.o : %.c
	$(CC) $(CFLAGS) -o $*.o $*.c

clean :
	rm -f $(TARGETS) *.o
