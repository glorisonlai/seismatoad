BUILDDIR=../build
INCLDIR=../include

LIBS= -lm -pthread

DEPS = main.c tsunameter.c

CC=mpicc
CFLAGS=-Wall

main: $(BUILDDIR) $(OBJ)
	$(CC) -o $(BUILDDIR)/$@ $(DEPS) $(LIBS)

$(BUILDDIR):
	mkdir -p $@

.PHONY: main