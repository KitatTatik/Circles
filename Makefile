CC=gcc
CFLAGS=-c -g -O0 -Wall -I.
LDLIBS=-lSDL2 -lSDL2_gfx -lm -lSDL2_image

SRC=meatball.c
OBJ=$(SRC:.c=.o)
EXE=mb

all: $(EXE)

$(EXE): $(OBJ)
	$(CC) $^ $(LDLIBS) -o $@

%.o : %.c $(DEP)
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f $(EXE) $(OBJ)
