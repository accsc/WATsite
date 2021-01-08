CC = gcc
CFLAGS = -std=c99

all:
	$(CC) $(CFLAGS) src/main.c src/clustering.c src/input.c src/util.c src/entropy.c src/enthalpy.c -lm -o bin/watsite
clean:
	rm *.o 

