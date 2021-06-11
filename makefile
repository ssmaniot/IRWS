CC := gcc 
CFLAGS := -std=gnu89 -Wall -pedantic -O3
LIBS := -lm

all: pagerank hits 

pagerank: pagerank.o
	$(CC) -o pagerank pagerank.o $(CFLAGS) $(LIBS)

hits: hits.o
	$(CC) -o hits hits.o $(CFLAGS) $(LIBS)

pagerank.o: src/pagerank.c 
	$(CC) -c src/pagerank.c $(CFLAGS)

hits.o: src/hits.c 
	$(CC) -c src/hits.c $(CFLAGS)

clean:
	rm -f *.o