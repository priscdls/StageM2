CC=gcc
CFLAGS=-Wall -g
LDFLAGS=-lR -lm
MEMFLAGS=-fsanitize=address -static-libasan
INCPATH=-I/usr/share/R/include
EXEC=test
SRC= $(wildcard *.c)
OBJ= $(SRC:.c=.o)

all: $(EXEC)

test: $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) $(INCPATH) -o $@ -c $< $(CFLAGS)

.PHONY: clean mrproper all mem

clean:
	rm -rf *.o

mrproper: clean
	rm -rf $(EXEC) Results