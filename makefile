CC = gcc 
CFLAGS = -std=gnu89 -Wall -pedantic -O3
LIBS = -lm

pagerank: pagerank.o
	$(CC) -o pagerank pagerank.o $(CFLAGS) $(LIBS)

pagerank.o: src/pagerank.c 
	$(CC) -c src/pagerank.c $(CFLAGS)

clean:
	rm -f *.o