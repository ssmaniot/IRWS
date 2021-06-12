CC := gcc 
override CFLAGS += -std=gnu89 -Wall -pedantic -O3
LDFLAGS := -lm
EXEC := pagerank hits

all: $(EXEC)

pagerank: pagerank.o
	$(CC) -o pagerank pagerank.o $(CFLAGS) $(LDFLAGS)
	
hits: hits.o
	$(CC) -o hits hits.o $(CFLAGS) $(LDFLAGS)

pagerank.o: src/pagerank.c 
	$(CC) -c src/pagerank.c $(CFLAGS)

hits.o: src/hits.c 
	$(CC) -c src/hits.c $(CFLAGS)

clean:
	rm -f *.o $(EXEC)