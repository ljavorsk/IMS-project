# IMS project
# file: Makefile
# 
# (C) Lukas Javorsky (xjavor20)
# (C) Patrik Ondriga (xondri08)

CC=gcc
CFLAGS=-lm -g
DEPS = 
OBJ = epidemic.o

.PHONY: clean, run

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

epidemic: $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^

run: epidemic
	./epidemic

clean:
	rm -f *.o epidemic
