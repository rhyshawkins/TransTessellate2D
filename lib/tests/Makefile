
CC = gcc
CFLAGS = -c -g -Wall -I..

LIBS=-L.. -ldelaunay -lm

TARGETS = test_voronoi_delaunay \
	test_linear_delaunay \
	test_ct_delaunay \
	test_max_edges

all : $(TARGETS)

test_voronoi_delaunay : test_voronoi_delaunay.o
	$(CC) -o test_voronoi_delaunay test_voronoi_delaunay.o $(LIBS)

test_linear_delaunay : test_linear_delaunay.o
	$(CC) -o test_linear_delaunay test_linear_delaunay.o $(LIBS)

test_ct_delaunay : test_ct_delaunay.o
	$(CC) -o test_ct_delaunay test_ct_delaunay.o $(LIBS)

test_max_edges : test_max_edges.o
	$(CC) -o test_max_edges test_max_edges.o $(LIBS)

%.o : %.c
	$(CC) $(CFLAGS) -o $*.o $*.c

clean :
	rm -f $(TARGETS) *.o


