Target=Main
QDWH=qdwh.c
Main=Main.c
OBJS=$(Main:.c=.o) $(QDWH:.c=.o)
CC=gcc
CFLAGS=-g -O0

all: $(Target)
$(Target):$(Main) $(QDWH)
	$(CC) $(CFLAGS) -c $(Main) $(QDWH)
	$(CC) -o $(Target) $(OBJS) -lblas -llapack -lm
	rm -f $(OBJS)

clean:
	rm -f $(Target)
