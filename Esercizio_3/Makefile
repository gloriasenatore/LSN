CC = g++
CFLAGS = -Wall -O3 --std=c++11

esercizio3.exe : esercizio3.o random.o funzioni.o
	$(CC) random.o funzioni.o esercizio3.o -o esercizio3.exe
esercizio3.o : esercizio3.cpp
	$(CC) -c esercizio3.cpp -o esercizio3.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
funzioni.o : funzioni.cpp funzioni.h
	$(CC) -c funzioni.cpp -o funzioni.o $(CFLAGS)

clean :
	rm *.o esercizio3.exe
