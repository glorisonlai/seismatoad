CC=mpicc
CFLAGS=-Wall

main: main.o tsunameter.o
	mkdir -p build
	$(CC) -o build/main main.o tsunameter.o
	rm -f ./*.o