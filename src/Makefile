BUILDDIR=../build
INCLDIR=../include

LIBS= -lm -pthread

DEPS = main.c satellite.c tsunameter.c

CC=mpicc
CFLAGS=-Wall

main: $(BUILDDIR) $(OBJ)
	$(CC) -o $(BUILDDIR)/$@ $(DEPS) $(LIBS)

.PHONY: main

clean:
	rm -f ./*.o