CC = g++
CFLAGS = -Wall -O3 --std=c++11

esercizio9.exe : esercizio9.o random.o funzioni.o
	$(CC) random.o funzioni.o esercizio9.o -o esercizio9.exe
esercizio9.o : esercizio9.cpp
	$(CC) -c esercizio9.cpp -o esercizio9.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
funzioni.o : funzioni.cpp funzioni.h
	$(CC) -c funzioni.cpp -o funzioni.o $(CFLAGS)

clean :
	rm *.o esercizio9.exe
