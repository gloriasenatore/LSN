CC = g++
CFLAGS = -Wall -O3 --std=c++11

esercizio10.1.exe : esercizio10.1.o random.o funzioni.o
	$(CC) random.o funzioni.o esercizio10.1.o -o esercizio10.1.exe
esercizio10.1.o : esercizio10.1.cpp
	$(CC) -c esercizio10.1.cpp -o esercizio10.1.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
funzioni.o : funzioni.cpp funzioni.h
	$(CC) -c funzioni.cpp -o funzioni.o $(CFLAGS)

clean :
	rm *.o esercizio10.1.exe
