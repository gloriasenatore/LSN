CC = g++
CFLAGS = -Wall -O3 --std=c++11

esercizio10.1.exe : esercizio10.1.o random.o funzioni.o
	$(CC) random.o funzioni.o esercizio10.1.o -o esercizio10.1.exe
esercizio10.1.o : esercizio10.1.cpp
	$(CC) -c esercizio10.1.cpp -o esercizio10.1.o $(CFLAGS)

main.exe : main.o funzioni.o random.o
	$(CC) random.o funzioni.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
funzioni.o : funzioni.cpp funzioni.h
	$(CC) -c funzioni.cpp -o funzioni.o $(CFLAGS)

clean :
	rm *.o main.exe seed.out
	rm *.o esercizio10.1.exe
